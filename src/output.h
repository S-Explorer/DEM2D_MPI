#ifndef OUTPUT_H
#define OUTPUT_H

#include "world.h"

#include <string>

/**
 * @brief VTU（XML UnstructuredGrid）帧写出器
 *        每个粒子写成单个 VTK_VERTEX 单元（中心，z=0），并附 point data
 *        （id / mid / radius / speed / omega / |force|）。在 ParaView 中用 Glyph
 *        过滤器（Glyph Type=Sphere，Scale Array=radius）按 radius 生成真实颗粒几何，
 *        并按 speed 等量上色。
 *        并行版：每个 rank 只写本地点 [0,n_local) 到自己的 frame_XXXXX_rr.vtu，
 *        由 rank 0 另写一份 .pvd 索引聚合各 rank 的 piece，ParaView 打开 .pvd 看全局动画。
 */
class VTKWriter {
public:
    bool write_frame(const World& world, const std::string& filename) const;
};

#endif
