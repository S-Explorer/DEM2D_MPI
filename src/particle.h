#ifndef PARTICLE_H
#define PARTICLE_H

#include "types.h"

/**
 * @brief 材料：弹性与接触参数
 */
struct Material {
    real_t E;       ///< 杨氏模量
    real_t nu;      ///< 泊松比
    real_t mu;      ///< 滑动摩擦系数
    real_t e;       ///< 恢复系数
    real_t density; ///< 密度
    real_t G;       ///< 剪切模量

    /**
     * @brief 由 E 与 nu 推导剪切模量 G
     */
    void init() { G = E / (2.0 * (1.0 + nu));}
};

/**
 * @brief 颗粒：位置、速度、受力及几何与惯性参数
 */
struct Particle {
    Vec2   pos = {0, 0};
    Vec2   vel = {0, 0};
    real_t omega = 0;
    Vec2   force = {0, 0};
    real_t torque = 0;
    real_t radius = 0;
    real_t mass = 0;
    real_t inertia = 0;
    id_t   id = 0;
    mid_t  mid = 0;

    void reset_force();
    void init_from_density(real_t density);
};




#endif
