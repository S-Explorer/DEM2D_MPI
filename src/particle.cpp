#include "particle.h"

/**
 * @brief 清零力与力矩
 */
void Particle::reset_force() {
    force = {0, 0};
    torque = 0;
}

/**
 * @brief 由密度与半径计算质量与转动惯量
 * @param density 密度
 */
void Particle::init_from_density(real_t density) {
    mass = density * pi * radius * radius;
    inertia = 0.5 * mass * radius * radius;
}
