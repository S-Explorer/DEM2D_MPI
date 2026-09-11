#include "grid.h"
#include <algorithm>

/**
 * @brief Grid的构造函数，主要是初始化网格的左下角右上角以及临近网格的大小
 *        计算临近网格数量
 * @param box_min 网格左下角
 * @param box_max 网格右上角
 * @param cell_size 临近网格尺寸
 */
Grid::Grid(Vec2 box_min, Vec2 box_max, real_t cell_size)
    : box_min(box_min), box_max(box_max), cell_size(cell_size) {
    nx = static_cast<int>((box_max.x - box_min.x) / cell_size) + 1;
    ny = static_cast<int>((box_max.y - box_min.y) / cell_size) + 1;
}

/**
 * @brief 计算当前颗粒所属的临近网格的索引
 * @param pos 颗粒的位置
 * @return 临近网格的索引
 */
int Grid::cell_index(const Vec2& pos) const {
    int cx = static_cast<int>((pos.x - box_min.x) / cell_size);
    int cy = static_cast<int>((pos.y - box_min.y) / cell_size);
    cx = std::max(0, std::min(nx - 1, cx));
    cy = std::max(0, std::min(ny - 1, cy));
    return cy * nx + cx;
}

/**
 * @brief 传入颗粒的数组来编排完整的索引链表
 * @param particles 颗粒的原始数组
 * @return void
 */
void Grid::build(const std::vector<Particle>& particles) {
    head.assign(nx * ny, -1);
    next.assign(particles.size(), -1);

    for (int i = 0; i < static_cast<int>(particles.size()); ++i) {
        int cell = cell_index(particles[i].pos);
        next[i] = head[cell];  // 头插法
        head[cell] = i;
    }
}

/**
 * @brief 判断当前pos是不是在目前的grid所在的region
 * @param pos 传入的位置
 * @return bool 是否在当前区域
 */
bool Grid::owns(const Vec2& pos) const {
    bool x_own{false}, y_own{false};
    if (box_min.x <= pos.x && box_max.x > pos.x) x_own = true;
    if (box_min.y <= pos.y && box_max.y > pos.y) y_own = true;
    return x_own && y_own;
}
