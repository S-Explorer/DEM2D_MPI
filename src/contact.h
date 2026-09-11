#ifndef CONTACT_H
#define CONTACT_H

#include "particle.h"

/**
 * @brief 粒子-粒子接触的历史状态
 */
struct Contact {
    id_t i, j;                 ///< 接触的两个粒子 ID，约定 i < j
    Vec2 t_disp = {0, 0};      ///< 累积切向位移（接触局部坐标）
    bool active = false;       ///< 当前是否处于接触
};

/**
 * @brief 粒子-粒子接触力模型
 *        Hertz 法向 + Mindlin 切向 + Tsuji 阻尼
 */
class ContactModel {
public:
    ContactModel(const Material* materials, int num_materials);

    void compute_force(Particle& a, Particle& b, Contact& c, real_t dt);

private:
    const Material* materials;
    int num_materials;
};

#endif
