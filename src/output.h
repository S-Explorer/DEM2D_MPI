#ifndef OUTPUT_H
#define OUTPUT_H

#include "world.h"

#include <string>

/**
 * @brief legacy VTK (ASCII) 帧写出器
 *        每个粒子写成单个 VTK_VERTEX 点（中心，z=0），并附 point data
 *        （id / mid / radius / speed / omega / |force|）。在 ParaView 中用 Glyph
 *        过滤器（Glyph Type=Sphere，Scale Array=radius）按 radius 生成真实颗粒几何，
 *        并按 speed 等量上色。时间序列由调用方按步编号生成文件名（如 frame_0000.vtk），
 *        ParaView 用 File > Open 选中首帧即可按文件名模式自动聚合为动画。
 *        串行版本：单进程直接写单个文件；MPI 版本后续可改为每 rank 写子域文件，
 *        或用 MPI-IO 并行写同一文件。
 */
class VTKWriter {
public:
    bool write_frame(const World& world, const std::string& filename) const;
};

#endif
