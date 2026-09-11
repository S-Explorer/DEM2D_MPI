#include "world.h"
#include "output.h"

#include <cstdio>
#include <fstream>
#include <random>
#include <vector>
// #include <print>
#include <mpi.h>
#include <string>

/**
 * @brief 串行 DEM 闭环测试
 *        计算域 [0,0]x[10,10]，四壁容器，80 个颗粒以抖动网格随机布置于上方，
 *        自由下落并垮塌成堆；每 200 步写一帧粒子 VTK（墙体不写入 VTK），
 *        并打印动能与 y 范围
 */
int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::printf("command example : exe pcols prows\n");
        return 1;
    }

    int pcols, prows;
    pcols = std::stoi(std::string(argv[1]));
    prows = std::stoi(std::string(argv[2]));

    MPI_Init(&argc, &argv);

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if ( pcols * prows == 0 || pcols * prows != size ) {
        if (rank == 0) std::printf("error prows and pcols\n");
        return 1;
    }

    // MPI 拓扑关系
    int dims[2] = {prows, pcols};
    int periods[2] = {0,0};
    int cart_rank, coords[2];
    int reorder = 1;
    MPI_Dims_create(size, 2, dims);
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);
    // 获取新的坐标和rank
    MPI_Comm_rank(cart_comm, &cart_rank);
    MPI_Cart_coords(cart_comm, cart_rank, 2, coords);

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
    const real_t r = 0.1;
    World world({0, 0}, {100, 10}, /*cell_size*/ 5*r, materials,
                /*gravity*/ {0, -9.81},
                cart_rank, /*px*/ dims[1], /*py*/ dims[0],
                /*me_x*/ coords[1], /*me_y*/ coords[0],
                cart_comm);

    // --- 容器四壁 ---
    world.add_wall({0.0,  0,  1});         // 左  x=0,  法线 +x
    world.add_wall({100.0, 0, -1});        // 右  x=100, 法线 -x
    world.add_wall({0.0,  1,  1});         // 底  y=0,  法线 +y
    world.add_wall({10.0, 1, -1});         // 顶  y=10, 法线 -y

    // --- 80 个颗粒 ---
    std::mt19937 rng(12345u);
    std::uniform_real_distribution<real_t> jitter(-1.0, 1.0);
    const int    cols = 240, rows = 30;     // 100:10 域，沿 x 多沿 y 少；sx/sy 须 >= 2r=0.2
    const real_t x0 = 1.0, x1 = 99.0;       // x 布置范围
    const real_t y0 = 2.0, y1 = 9.0;       // y 布置范围
    const real_t sx  = (x1 - x0) / cols;   // ≈0.408
    const real_t sy  = (y1 - y0) / rows;   // ≈0.233
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
    const int    total = 15000;
    const int    every = 500;
    VTKWriter    writer;

    std::printf("DEM2D smoke test: particles=%d  crit_dt=%.3e  dt=%.3e  steps=%d\n",
                static_cast<int>(world.get_particles().size()), crit, dt, total);

    for (int s = 0; s < total; ++s) {
        world.step(dt);
        if (s % every == 0) {
            // 规约颗粒数
            auto n_local = world.get_nlocal();
            auto n_global = n_local;
            MPI_Allreduce(&n_local, &n_global, 1, MPI_INT, MPI_SUM, cart_comm);
            if (rank == 0) {
                std::printf("all the particles number %05d\n", n_global);
            }
            // 写存储文件
            char buf[64];
            std::snprintf(buf, sizeof(buf), "frame_%05d_%02d.vtu", s, cart_rank);
            writer.write_frame(world, buf);
        }
    }
    // rank 0 写 .pvd 索引，聚合所有 rank 的 .vtu（文件名按约定 frame_XXXXX_rr.vtu）
    if (rank == 0) {
        std::ofstream pvd("dem2d.pvd");
        pvd << "<?xml version=\"1.0\"?>\n"
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
            << "  <Collection>\n";
        for (int fs = 0; fs < total; fs += every) {
            double t = fs * dt;
            for (int r = 0; r < size; ++r) {
                char fbuf[64];
                std::snprintf(fbuf, sizeof(fbuf), "frame_%05d_%02d.vtu", fs, r);
                pvd << "    <DataSet timestep=\"" << t << "\" group=\"\" part=\""
                    << r << "\" file=\"" << fbuf << "\"/>\n";
            }
        }
        pvd << "  </Collection>\n</VTKFile>\n";
    }

    if (cart_rank == 0) std::printf("done.\n");

    MPI_Finalize();
    return 0;
}
