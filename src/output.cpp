#include "output.h"

#include <fstream>

/**
 * @brief 写一帧粒子数据到文件
 * @param world 世界
 * @param filename 输出文件名
 * @return 写入成功返回 true
 */
bool VTKWriter::write_frame(const World& world, const std::string& filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    out.precision(10);   // 10 位有效数字，坐标/速度够用

    const std::vector<Particle>& particles = world.get_particles();
    const int np = static_cast<int>(particles.size());

    // ---- VTK legacy (ASCII)，UNSTRUCTURED_GRID（便于后续混入墙体线单元）----
    out << "# vtk DataFile Version 3.0\n";
    out << "DEM2D frame\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";

    // 每个粒子一个点（中心，z=0），单元类型 VTK_VERTEX。
    // 在 ParaView 中用 Glyph（Sphere，Scale Array=radius）生成真实颗粒几何。
    out << "POINTS " << np << " double\n";
    for (const Particle& p : particles)
        out << p.pos.x << ' ' << p.pos.y << ' ' << real_t(0) << '\n';

    out << "CELLS " << np << ' ' << (np * 2) << '\n';
    for (int i = 0; i < np; ++i)
        out << "1 " << i << '\n';

    out << "CELL_TYPES " << np << '\n';
    for (int i = 0; i < np; ++i) out << "1\n";   // VTK_VERTEX = 1

    // ---- 粒子属性作为 point data（每点一值），便于 Glyph 按 radius 缩放、按速度等上色 ----
    out << "POINT_DATA " << np << '\n';

    auto scalar_int = [&](const char* name, auto getter) {
        out << "SCALARS " << name << " int 1\nLOOKUP_TABLE default\n";
        for (const Particle& p : particles) out << getter(p) << '\n';
    };
    auto scalar_real = [&](const char* name, auto getter) {
        out << "SCALARS " << name << " double 1\nLOOKUP_TABLE default\n";
        for (const Particle& p : particles) out << getter(p) << '\n';
    };

    scalar_int ("id",        [](const Particle& p){ return p.id; });
    scalar_int ("mid",       [](const Particle& p){ return p.mid; });
    scalar_real("radius",    [](const Particle& p){ return p.radius; });
    scalar_real("speed",     [](const Particle& p){ return p.vel.length(); });
    scalar_real("omega",     [](const Particle& p){ return p.omega; });
    scalar_real("force_mag", [](const Particle& p){ return p.force.length(); });

    return static_cast<bool>(out);
}
