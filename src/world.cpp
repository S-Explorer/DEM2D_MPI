#include "world.h"

#include <cmath>
#include <limits>
#include <mpi.h>

// 单进程 DEM 时间步进：半隐式欧拉积分 + Hertz-Mindlin 接触
// （粒子-粒子、粒子-墙）。每步三阶段：
//   1) grid.build      重建邻域网格
//   2) compute_force   清零力 -> 计算接触力 -> 修剪分离接触
//   3) integrate       半隐式欧拉推进速度 / 位置

/**
 * @brief 粒子-粒子接触历史的查表键，打包为 64 位
 * @param i 粒子索引（较小者）
 * @param j 粒子索引（较大者）
 * @return 64 位查表键
 */
static inline uint64_t pp_contact_key(id_t i, id_t j) {
    return (static_cast<uint64_t>(i) << 32) | static_cast<uint32_t>(j);
}

/**
 * @brief 构造世界
 * @param box_min 模拟域左下角
 * @param box_max 模拟域右上角
 * @param cell_size 临近网格尺寸，需不小于最大粒子直径以保证 3x3 邻域搜索完整
 * @param mats 材料表，构造后不再变动
 * @param grav 均匀体力加速度
 */
World::World(Vec2 box_min, Vec2 box_max, real_t cell_size,
            const std::vector<Material>& mats,
            Vec2 grav,
            const int rank, const int px, const int py, const int me_x, const int me_y,
            MPI_Comm cart_comm)
    : materials(mats),
      g_min(box_min), g_max(box_max),
      gravity(grav),
      contact_model(materials.data(), static_cast<int>(materials.size())),
      wall_model(materials.data(), static_cast<int>(materials.size())),
      cart_rank(rank),
      cart_comm(cart_comm),
      px(px),
      py(py),
      me_x(me_x),
      me_y(me_y),
      grid(
          {box_min.x + (box_max.x - box_min.x) * me_x / px,
           box_min.y + (box_max.y - box_min.y) * me_y / py},
          {box_min.x + (box_max.x - box_min.x) * (me_x + 1) / px,
           box_min.y + (box_max.y - box_min.y) * (me_y + 1) / py},
          cell_size
      )
      {
          MPI_Cart_shift(cart_comm, 1, 1, &left, &right);
          MPI_Cart_shift(cart_comm, 0, 1, &down, &up);
      }

/**
 * @brief 添加粒子
 * @param p 粒子
 * @return 粒子本地索引
 */
int World::add_particle(const Particle& p) {

    if (!grid.owns(p.pos)){return -1;}

    int idx = static_cast<int>(particles.size());
    particles.push_back(p);
    n_local = idx + 1;
    return idx;
}

/**
 * @brief 添加墙
 * @param w 墙
 * @return 墙索引
 */
int World::add_wall(const Wall& w) {
    int idx = static_cast<int>(walls.size());
    walls.push_back(w);
    return idx;
}

int World::get_nlocal() const {
    return n_local;
}

/**
 * @brief Rayleigh 时间步，取所有粒子中最小者作为最保守的稳定步长上限
 * @return 时间步长
 */
real_t World::critical_timestep() const {
    // Rayleigh 时间步：t_R = π R sqrt(ρ/G) / (0.1631 ν + 0.8766)
    // 取所有粒子中最小者，作为最保守的稳定步长上限。
    // if (particles.empty()) return real_t(0);
    real_t t_min = std::numeric_limits<real_t>::max();
    for (int i = 0; i < n_local; i++) {
        auto& p = particles[i];
        const Material& m = materials[p.mid];
        real_t t_R = pi * p.radius * std::sqrt(m.density / m.G)
                   / (real_t(0.1631) * m.nu + real_t(0.8766));
        if (t_R < t_min) t_min = t_R;
    }
    double t_global = t_min;
    MPI_Allreduce(&t_min, &t_global, 1, MPI_DOUBLE, MPI_MIN, cart_comm);

    return t_global;
}

/**
 * @brief 按给定步数运行模拟
 * @param num_steps 步数
 * @param dt 时间步长，小于等于 0 时按 critical_timestep 的 0.2 倍自动选取
 */
void World::run(int num_steps, real_t dt) {
    if (dt <= 0) dt = critical_timestep() * real_t(0.2);    // 安全系数 0.2
    if (dt <= 0) return;                                    // 空世界
    for (int s = 0; s < num_steps; ++s)
        step(dt);
}

/**
 * @brief 单步推进：重建邻域网格、计算接触力、时间积分
 * @param dt 时间步长
 */
void World::step(real_t dt) {
    current_dt = dt;
    migrate();
    exchange_ghosts();
    grid.build(particles);   // 阶段 1：邻域搜索
    compute_force();         // 阶段 2：接触力
    integrate(dt);           // 阶段 3：时间积分
}

/**
 * @brief 计算所有接触力并累加到粒子，随后修剪已分离的接触
 */
void World::compute_force() {
    for (int i = 0; i < n_local; i++ ) {
        auto& p = particles[i];
        p.reset_force();
    }

    // 先把现有接触标记为“本步待确认”，便于后续修剪已分离者
    for (auto& kv : contacts)      kv.second.active = false;
    for (auto& kv : wall_contacts) kv.second.active = false;

    for (int i = 0; i < n_local; ++i) {
        // --- 粒子-粒子 ---
        grid.for_each_neighbor(i, particles, [&](int j) {
            if (j <= i) return;                       // 每对只算一次 (i < j)
            Vec2 d = particles[j].pos - particles[i].pos;
            real_t rsum = particles[i].radius + particles[j].radius;
            if (d.length_squared() >= rsum * rsum) return;   // 预筛：不接触则不查表
            id_t id_i = particles[i].id, id_j = particles[j].id;
            if (id_i > id_j) std::swap(id_i, id_j);
            Contact& c = contacts[pp_contact_key(id_i, id_j)];
            c.i = id_i; c.j = id_j;
            contact_model.compute_force(particles[i], particles[j], c, current_dt);
        });
        // --- 粒子-墙 ---
        for (int wi = 0; wi < static_cast<int>(walls.size()); ++wi) {
            if (wall_overlap(particles[i], walls[wi]) <= 0) continue;  // 预筛
            WallContact& wc = wall_contacts[wall_contact_key(particles[i].id, wi)];
            wall_model.compute_force(particles[i], walls[wi], wc, current_dt);
        }
    }

    prune_contacts();   // 修剪本步未确认（已分离）的接触
}

/**
 * @brief 半隐式欧拉积分：先用当前力与重力更新速度，再用新速度推进位置
 * @param dt 时间步长
 */
void World::integrate(real_t dt) {
    for (int i = 0; i < n_local; i++) {
        auto& p = particles[i];
        real_t inv_m = real_t(1) / p.mass;
        real_t inv_I = real_t(1) / p.inertia;
        p.vel   += p.force * (inv_m * dt) + gravity * dt;
        p.omega += p.torque * (inv_I * dt);
        p.pos   += p.vel * dt;
    }
}

/**
 * @brief 修剪本步未确认（已分离）的粒子-粒子与粒子-墙接触
 */
void World::prune_contacts() {
    for (auto it = contacts.begin(); it != contacts.end();) {
        if (it->second.active) ++it;
        else it = contacts.erase(it);
    }
    for (auto it = wall_contacts.begin(); it != wall_contacts.end();) {
        if (it->second.active) ++it;
        else it = wall_contacts.erase(it);
    }
}

void World::exchange_ghosts(){
    const int TAG_X_CNT = 0, TAG_X_DAT = 1;
    const int TAG_Y_CNT = 2, TAG_Y_DAT = 3;

    // 打包函数
    auto pack = [&](std::vector<Particle>& buf, int end, auto pred) {
        buf.clear();
        for ( int k = 0; k < end; k++) {
            if (pred(particles[k])) buf.push_back(particles[k]);
        }
    };

    // 邻居数据交换
    auto exchange_axis = [&](int lo, int hi, int tag_cnt, int tag_dat,
                             int pack_end, auto pred_lo, auto pred_hi){
        // lo and hi pack
        if ( lo != MPI_PROC_NULL) pack(send_buf[0], pack_end, pred_lo);
        else                      send_buf[0].clear();
        if ( hi != MPI_PROC_NULL) pack(send_buf[1], pack_end, pred_hi);
        else                      send_buf[1].clear();

        //
        int cnt_lo = (int)send_buf[0].size();
        int cnt_hi = (int)send_buf[1].size();
        int recv_lo = 0;
        int recv_hi = 0;

        // 缓冲区大小确定
        MPI_Request r[4]; int nr = 0;
        if (lo != MPI_PROC_NULL) {
            MPI_Isend(&cnt_lo,  1, MPI_INT, lo, tag_cnt, cart_comm, &r[nr++]);
            MPI_Irecv(&recv_lo, 1, MPI_INT, lo, tag_cnt, cart_comm, &r[nr++]);
        }
        if (hi != MPI_PROC_NULL) {
            MPI_Isend(&cnt_hi,  1, MPI_INT, hi, tag_cnt, cart_comm, &r[nr++]);
            MPI_Irecv(&recv_hi, 1, MPI_INT, hi, tag_cnt, cart_comm, &r[nr++]);
        }

        if (nr) MPI_Waitall(nr, r, MPI_STATUSES_IGNORE);

        // 交换缓冲区数据
        recv_buf[0].resize(recv_lo);
        recv_buf[1].resize(recv_hi);
        nr = 0;
        if (lo != MPI_PROC_NULL) {
            MPI_Isend(send_buf[0].data(), (int)(cnt_lo * sizeof(Particle)),
                MPI_BYTE, lo, tag_dat, cart_comm, &r[nr++]);
            MPI_Irecv(recv_buf[0].data(), (int)(recv_lo * sizeof(Particle)),
                MPI_BYTE, lo, tag_dat, cart_comm, &r[nr++]);
        }
        if (hi != MPI_PROC_NULL) {
            MPI_Isend(send_buf[1].data(), (int)(cnt_hi * sizeof(Particle)),
                MPI_BYTE, hi, tag_dat, cart_comm, &r[nr++]);
            MPI_Irecv(recv_buf[1].data(), (int)(recv_hi * sizeof(Particle)),
                MPI_BYTE, hi, tag_dat, cart_comm, &r[nr++]);
        }

        if (nr) MPI_Waitall(nr, r, MPI_STATUSES_IGNORE);
        // 追加 ghost 数据到particles中
        for (const Particle& p: recv_buf[0]) particles.push_back(p);
        for (const Particle& p: recv_buf[1]) particles.push_back(p);

    };

    // left right — 只用本地粒子（索引 0..n_local-1）
    exchange_axis(left, right, TAG_X_CNT, TAG_X_DAT, n_local,
        [&](const Particle& p){return p.pos.x < grid.getmin().x + grid.getsize();},
        [&](const Particle& p){return p.pos.x > grid.getmax().x - grid.getsize();}
    );
    // down up — 包含本地 + X方向 ghost，确保角落 ghost 完整
    exchange_axis(down, up,    TAG_Y_CNT, TAG_Y_DAT, (int)particles.size(),
        [&](const Particle& p){return p.pos.y < grid.getmin().y + grid.getsize();},
        [&](const Particle& p){return p.pos.y > grid.getmax().y - grid.getsize();}
    );
}


void World::migrate() {
    const int TAG_X_CNT = 4, TAG_X_DAT = 5;
    const int TAG_Y_CNT = 6, TAG_Y_DAT = 7;

    auto migrate_axis = [&](int lo, int hi, int tag_cnt, int tag_dat,
                            int pred_end, auto pred_lo, auto pred_hi){
        send_buf[0].clear();
        send_buf[1].clear();
        int k = 0;
        for ( int i = 0; i < pred_end; i++){
            Particle& p = particles[i];
            if (lo != MPI_PROC_NULL && pred_lo(p))      send_buf[0].push_back(std::move(p));
            else if (hi != MPI_PROC_NULL && pred_hi(p)) send_buf[1].push_back((std::move(p)));
            else {
                if (k != i) particles[k] = std::move(p);
                k++;
            }
        }

        int cnt_lo = (int)send_buf[0].size(), cnt_hi = (int)send_buf[1].size();
        int recv_lo = 0, recv_hi = 0;

        // buf size
        MPI_Request r[4]; int nr = 0;
        if (lo != MPI_PROC_NULL) {
            MPI_Isend(&cnt_lo,  1, MPI_INT, lo, tag_cnt, cart_comm, &r[nr++]);
            MPI_Irecv(&recv_lo, 1, MPI_INT, lo, tag_cnt, cart_comm, &r[nr++]);
        }
        if (hi != MPI_PROC_NULL) {
            MPI_Isend(&cnt_hi,  1, MPI_INT, hi, tag_cnt, cart_comm, &r[nr++]);
            MPI_Irecv(&recv_hi, 1, MPI_INT, hi, tag_cnt, cart_comm, &r[nr++]);
        }
        if (nr) MPI_Waitall(nr, r, MPI_STATUSES_IGNORE);

        // 交换buf
        recv_buf[0].resize(recv_lo);
        recv_buf[1].resize(recv_hi);
        nr = 0;
        if (lo != MPI_PROC_NULL) {
            MPI_Isend(send_buf[0].data(), (int)(cnt_lo * sizeof(Particle)),
                MPI_BYTE, lo, tag_dat, cart_comm, &r[nr++]);
            MPI_Irecv(recv_buf[0].data(), (int)(recv_lo * sizeof(Particle)),
                MPI_BYTE, lo, tag_dat, cart_comm, &r[nr++]);
        }
        if (hi != MPI_PROC_NULL) {
            MPI_Isend(send_buf[1].data(), (int)(cnt_hi * sizeof(Particle)),
                MPI_BYTE, hi, tag_dat, cart_comm, &r[nr++]);
            MPI_Irecv(recv_buf[1].data(), (int)(recv_hi * sizeof(Particle)),
                MPI_BYTE, hi, tag_dat, cart_comm, &r[nr++]);
        }
        if (nr) MPI_Waitall(nr, r, MPI_STATUSES_IGNORE);

        // local + recv
        particles.resize(k);
        for (auto& p: recv_buf[0]) particles.push_back(std::move(p));
        for (auto& p: recv_buf[1]) particles.push_back(std::move(p));
    };

    // left right
    migrate_axis(left, right, TAG_X_CNT, TAG_X_DAT, n_local,
        [&](const Particle& p){return p.pos.x < grid.getmin().x;},
        [&](const Particle& p){return p.pos.x >= grid.getmax().x;});
    // down up
    migrate_axis(down, up, TAG_Y_CNT, TAG_Y_DAT, (int)(particles.size()),
        [&](const Particle& p){return p.pos.y < grid.getmin().y;},
        [&](const Particle& p){return p.pos.y >= grid.getmax().y;});

    n_local = (int)particles.size();

}
