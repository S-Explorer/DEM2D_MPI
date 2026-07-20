#include "world.h"

#include <cmath>
#include <limits>

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
static inline uint64_t pp_contact_key(int i, int j) {
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
             const std::vector<Material>& mats, Vec2 grav)
    : materials(mats),
      grid(box_min, box_max, cell_size),
      gravity(grav),
      contact_model(materials.data(), static_cast<int>(materials.size())),
      wall_model(materials.data(), static_cast<int>(materials.size())) {}

/**
 * @brief 添加粒子
 * @param p 粒子
 * @return 粒子本地索引
 */
int World::add_particle(const Particle& p) {
    int idx = static_cast<int>(particles.size());
    particles.push_back(p);
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

/**
 * @brief Rayleigh 时间步，取所有粒子中最小者作为最保守的稳定步长上限
 * @return 时间步长
 */
real_t World::critical_timestep() const {
    // Rayleigh 时间步：t_R = π R sqrt(ρ/G) / (0.1631 ν + 0.8766)
    // 取所有粒子中最小者，作为最保守的稳定步长上限。
    if (particles.empty()) return real_t(0);
    real_t t_min = std::numeric_limits<real_t>::max();
    for (const Particle& p : particles) {
        const Material& m = materials[p.mid];
        real_t t_R = pi * p.radius * std::sqrt(m.density / m.G)
                   / (real_t(0.1631) * m.nu + real_t(0.8766));
        if (t_R < t_min) t_min = t_R;
    }
    return t_min;
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
    grid.build(particles);   // 阶段 1：邻域搜索
    compute_force();         // 阶段 2：接触力
    integrate(dt);           // 阶段 3：时间积分
}

/**
 * @brief 计算所有接触力并累加到粒子，随后修剪已分离的接触
 */
void World::compute_force() {
    for (Particle& p : particles) p.reset_force();

    // 先把现有接触标记为“本步待确认”，便于后续修剪已分离者
    for (auto& kv : contacts)      kv.second.active = false;
    for (auto& kv : wall_contacts) kv.second.active = false;

    const int n = static_cast<int>(particles.size());
    for (int i = 0; i < n; ++i) {
        // --- 粒子-粒子 ---
        grid.for_each_neighbor(i, particles, [&](int j) {
            if (j <= i) return;                       // 每对只算一次 (i < j)
            Vec2 d = particles[j].pos - particles[i].pos;
            real_t rsum = particles[i].radius + particles[j].radius;
            if (d.length_squared() >= rsum * rsum) return;   // 预筛：不接触则不查表
            Contact& c = contacts[pp_contact_key(i, j)];
            c.i = i; c.j = j;
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
    for (Particle& p : particles) {
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
