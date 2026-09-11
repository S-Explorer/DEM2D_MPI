#ifndef GRID_H
#define GRID_H

#include "types.h"
#include "particle.h"
#include <vector>

/**
 * @brief 均匀网格，基于索引链表实现粒子的邻域搜索
 */
class Grid {
public:
    Grid(Vec2 box_min, Vec2 box_max, real_t cell_size);

    void build(const std::vector<Particle>& particles);
    bool owns(const Vec2& pos) const;
    Vec2 getmin() const {return box_min;}
    Vec2 getmax() const {return box_max;}
    auto getsize() const {return cell_size;}

    /**
     * @brief 遍历粒子 i 的 3x3 邻域，对每个邻居索引 j 调用回调
     * @tparam F 邻居回调类型，签名为 void(int)
     * @param i 中心粒子索引
     * @param particles 粒子数组
     * @param callback 接收邻居索引 j 的回调
     */
    template<typename F>
    void for_each_neighbor(int i, const std::vector<Particle>& particles, F&& callback) {
        const Vec2& p = particles[i].pos;
        int cx = static_cast<int>((p.x - box_min.x) / cell_size);
        int cy = static_cast<int>((p.y - box_min.y) / cell_size);
        cx = (cx < 0) ? 0 : (cx >= nx ? nx - 1 : cx);
        cy = (cy < 0) ? 0 : (cy >= ny ? ny - 1 : cy);

        for (int dy = -1; dy <= 1; ++dy) {
            int yy = cy + dy;
            if (yy < 0 || yy >= ny) continue;
            for (int dx = -1; dx <= 1; ++dx) {
                int xx = cx + dx;
                if (xx < 0 || xx >= nx) continue;
                int cell = yy * nx + xx;
                for (int j = head[cell]; j != -1; j = next[j]) {
                    if (j == i) continue;
                    callback(j);
                }
            }
        }
    }

private:
    Vec2 box_min, box_max;
    real_t cell_size;
    int nx, ny;

    std::vector<int> head;  ///< head[cell] 指向首个粒子，-1 表示空
    std::vector<int> next;  ///< next[particle] 指向同 cell 下一粒子，-1 表示链尾

    int cell_index(const Vec2& pos) const;
};

#endif
