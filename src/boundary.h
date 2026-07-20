#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "types.h"
#include "particle.h"

/**
 * @brief 轴对齐刚体墙（无限大质量），粒子在墙的单侧与之碰撞
 *        axis=0 为竖直墙 x=pos，axis=1 为水平墙 y=pos；
 *        dir 取 ±1 为外法线方向（指向粒子合法所在一侧），
 *        axis=0 时 n=(dir,0)，axis=1 时 n=(0,dir)
 */
struct Wall {
    real_t   pos;
    uint16_t axis;
    int16_t  dir;
};

/**
 * @brief 单个 粒子-墙 接触的历史状态
 *        墙只有一侧，仅切向位移需要跨步累积
 */
struct WallContact {
    Vec2 t_disp = {0, 0};   ///< 累积切向位移（接触局部坐标）
    bool active = false;    ///< 当前是否处于接触
};

/**
 * @brief 粒子-墙接触力模型，复用 Hertz-Mindlin + Tsuji 阻尼（见 contact.cpp）
 *        墙为刚体（E_wall、G_wall、曲率半径均趋于无穷），等效参数退化为只含粒子一侧。
 *        仅持有材料指针，无内部状态；切向历史由调用方持有（如 World 的墙接触表），
 *        热路径零分配，内存友好。
 */
class WallContactModel {
public:
    WallContactModel(const Material* materials, int num_materials);

    void compute_force(Particle& p, const Wall& w, WallContact& wc, real_t dt);

private:
    const Material* materials;
    int num_materials;
};

/**
 * @brief 粒子-墙接触历史的查表键，将 (粒子 id, 墙索引) 打包为 64 位
 *        供 World 用 unordered_map 持久化切向历史
 * @param particle_id 粒子 id
 * @param wall_index 墙索引
 * @return 64 位查表键
 */
inline uint64_t wall_contact_key(id_t particle_id, int wall_index) {
    return (static_cast<uint64_t>(particle_id) << 32) |
           static_cast<uint64_t>(static_cast<uint32_t>(wall_index));
}

/**
 * @brief 粒子与墙的法向重叠量，大于 0 表示接触
 *        供调用方预筛，避免为不接触的 (粒子,墙) 对查/建历史表
 * @param p 粒子
 * @param w 墙
 * @return 法向重叠量
 */
inline real_t wall_overlap(const Particle& p, const Wall& w) {
    real_t s = (w.axis == 0) ? (p.pos.x - w.pos) * static_cast<real_t>(w.dir)
                             : (p.pos.y - w.pos) * static_cast<real_t>(w.dir);
    return p.radius - s;
}

#endif // BOUNDARY_H
