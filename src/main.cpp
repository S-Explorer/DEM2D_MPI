#include "world.h"
#include "output.h"

#include <cstdio>
#include <random>
#include <vector>

/**
 * @brief 串行 DEM 闭环测试
 *        计算域 [0,0]x[10,10]，四壁容器，80 个颗粒以抖动网格随机布置于上方，
 *        自由下落并垮塌成堆；每 200 步写一帧粒子 VTK（墙体不写入 VTK），
 *        并打印动能与 y 范围
 */
int main() {
    // --- 材料（玻璃样）---
    Material mat;
    mat.E       = 1e9;
    mat.nu      = 0.3;
    mat.mu      = 0.3;
    mat.e       = 0.8;
    mat.density = 2500.0;
    mat.init();                            // G = E / (2(1+nu))
    std::vector<Material> materials = {mat};

    // --- 世界 ---
    const real_t r = 0.2;
    World world({0, 0}, {10, 10}, /*cell_size*/ 0.5, materials,
                /*gravity*/ {0, -9.81});

    // --- 容器四壁 ---
    world.add_wall({0.0,  0,  1});         // 左  x=0,  法线 +x
    world.add_wall({10.0, 0, -1});         // 右  x=10, 法线 -x
    world.add_wall({0.0,  1,  1});         // 底  y=0,  法线 +y
    world.add_wall({10.0, 1, -1});         // 顶  y=10, 法线 -y

    // --- 80 个颗粒 ---
    std::mt19937 rng(12345u);
    std::uniform_real_distribution<real_t> jitter(-1.0, 1.0);
    const int    cols = 10, rows = 8;
    const real_t x0 = 2.0, x1 = 8.0;       // x 布置范围
    const real_t y0 = 4.0, y1 = 9.5;       // y 布置范围
    const real_t sx  = (x1 - x0) / cols;   // 0.6
    const real_t sy  = (y1 - y0) / rows;   // 0.6875
    const real_t jx  = (sx - 2.0 * r) * 0.5 * 0.9;
    const real_t jy  = (sy - 2.0 * r) * 0.5 * 0.9;

    id_t id = 0;
    for (int j = 0; j < rows; ++j)
        for (int i = 0; i < cols; ++i) {
            Particle p;
            p.radius = r;
            p.mid    = 0;
            p.init_from_density(materials[0].density);
            p.pos    = { x0 + (i + 0.5) * sx + jx * jitter(rng),
                         y0 + (j + 0.5) * sy + jy * jitter(rng) };
            p.vel    = {0, 0};
            p.id     = id++;
            world.add_particle(p);
        }

    // --- 运行 ---
    const real_t crit  = world.critical_timestep();
    const real_t dt    = crit * 0.2;
    const int    total = 12000;
    const int    every = 200;
    VTKWriter    writer;

    std::printf("DEM2D smoke test: particles=%d  crit_dt=%.3e  dt=%.3e  steps=%d\n",
                static_cast<int>(world.get_particles().size()), crit, dt, total);

    for (int s = 0; s < total; ++s) {
        world.step(dt);
        if (s % every == 0) {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "frame_%05d.vtk", s);
            writer.write_frame(world, buf);

            real_t ke = 0, ymin = 1e9, ymax = -1e9;
            for (const Particle& p : world.get_particles()) {
                ke += real_t(0.5) * p.mass * p.vel.length_squared()
                    + real_t(0.5) * p.inertia * p.omega * p.omega;
                if (p.pos.y < ymin) ymin = p.pos.y;
                if (p.pos.y > ymax) ymax = p.pos.y;
            }
            std::printf("step %5d  t=%6.3f  KE=%.3e  y=[%.3f,%.3f]\n",
                        s, s * dt, ke, ymin, ymax);
        }
    }
    std::printf("done.\n");
    return 0;
}
