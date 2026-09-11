#include "output.h"

#include <fstream>

/**
 * @brief 写一帧本地点数据到 VTU 文件
 * @param world 世界
 * @param filename 输出文件名（每 rank 各自的 frame_XXXXX_rr.vtu）
 * @return 写入成功返回 true
 */
bool VTKWriter::write_frame(const World& world, const std::string& filename) const {
    std::ofstream out(filename);
    if (!out) return false;
    out.precision(10);

    const std::vector<Particle>& particles = world.get_particles();
    const int np = world.get_nlocal();          // 只写本地点，不含 ghost

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
        << "  <UnstructuredGrid>\n"
        << "    <Piece NumberOfPoints=\"" << np << "\" NumberOfCells=\"" << np << "\">\n";

    // ---- Points：3 分量，z=0 ----
    out << "      <Points>\n"
        << "        <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int i = 0; i < np; ++i)
        out << "          " << particles[i].pos.x << " " << particles[i].pos.y << " 0\n";
    out << "        </DataArray>\n      </Points>\n";

    // ---- Cells：每个粒子一个 VTK_VERTEX (type=1) ----
    out << "      <Cells>\n";
    out << "        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n          ";
    for (int i = 0; i < np; ++i) out << i << " ";
    out << "\n        </DataArray>\n";
    out << "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ";
    for (int i = 0; i < np; ++i) out << (i + 1) << " ";
    out << "\n        </DataArray>\n";
    out << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ";
    for (int i = 0; i < np; ++i) out << "1 ";
    out << "\n        </DataArray>\n      </Cells>\n";

    // ---- PointData：一个块里多个属性 ----
    out << "      <PointData>\n";
    auto scalar_real = [&](const char* name, auto getter) {
        out << "        <DataArray type=\"Float64\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (int i = 0; i < np; ++i) out << "          " << getter(particles[i]) << "\n";
        out << "        </DataArray>\n";
    };
    auto scalar_int = [&](const char* name, auto getter) {
        out << "        <DataArray type=\"Int32\" Name=\"" << name << "\" format=\"ascii\">\n";
        for (int i = 0; i < np; ++i) out << "          " << getter(particles[i]) << "\n";
        out << "        </DataArray>\n";
    };
    scalar_int ("id",        [](const Particle& p){ return p.id; });
    scalar_int ("mid",       [](const Particle& p){ return p.mid; });
    scalar_real("radius",    [](const Particle& p){ return p.radius; });
    scalar_real("speed",     [](const Particle& p){ return p.vel.length(); });
    scalar_real("omega",     [](const Particle& p){ return p.omega; });
    scalar_real("force_mag", [](const Particle& p){ return p.force.length(); });
    out << "      </PointData>\n";

    out << "    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n";
    return static_cast<bool>(out);
}
