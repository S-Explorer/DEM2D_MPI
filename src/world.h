#ifndef WORLD_H
#define WORLD_H

#include "boundary.h"
#include "contact.h"
#include "grid.h"

#include <vector>
#include <unordered_map>

/**
 * @brief DEM 模拟世界：持有粒子、墙、材料与接触历史，负责时间步进
 */
class World {
public:
    World(Vec2 box_min, Vec2 box_max, real_t cell_size,
          const std::vector<Material>& materials,
          Vec2 gravity = {0, 0});

    int  add_particle(const Particle& p);
    int  add_wall(const Wall& w);

    void step(real_t dt);
    void run(int num_steps, real_t dt = 0);
    real_t critical_timestep() const;

    const std::vector<Particle>& get_particles() const { return particles; }
    const std::vector<Wall>&     get_walls()     const { return walls; }
    const std::vector<Material>& get_materials() const { return materials; }

private:
    void compute_force();
    void integrate(real_t dt);
    void prune_contacts();

    std::vector<Particle> particles;
    std::vector<Material> materials;
    Grid                  grid;
    Vec2                  gravity = {0, 0};
    real_t                current_dt = 0;   ///< 供 compute_force 使用的时间步长

    std::unordered_map<uint64_t, Contact>     contacts;       ///< 粒子-粒子接触历史
    std::unordered_map<uint64_t, WallContact> wall_contacts;  ///< 粒子-墙接触历史

    std::vector<Wall>     walls;

    ContactModel     contact_model;
    WallContactModel wall_model;
};

#endif
