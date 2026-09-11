# DEM2D_MPI

基于 C++17 与 MPI 的二维离散元方法（DEM, Discrete Element Method）并行求解器，用于模拟大量圆形颗粒在重力与容器约束下的运动、碰撞与堆积过程。

## 功能特性

### 物理模型
- **颗粒-颗粒接触力模型**：Hertz 法向接触 + Mindlin 切向接触 + Tsuji 阻尼，支持切向位移的跨步累积历史（`src/contact.h`）
- **颗粒-墙接触力模型**：复用 Hertz-Mindlin + Tsuji 阻尼，墙体为刚体（无限大质量），支持任意轴对齐单侧墙（`src/boundary.h`）
- **材料参数**：支持多材料定义（杨氏模量、泊松比、摩擦系数、恢复系数、密度）
- **时间步长**：自动计算临界时间步长（`World::critical_timestep`）

### 空间离散与邻域搜索
- **均匀网格**：基于索引链表的 cell-list 邻域搜索，3×3 邻域遍历（`src/grid.h`）

### MPI 并行（`src/world.h` / `src/world.cpp`）
- **二维笛卡尔拓扑**：通过 `MPI_Cart_create` 将进程组织为 pcols×prows 的网格，计算域按进程网格划分
- **Ghost 粒子交换**：进程边界处与相邻进程交换 ghost 层粒子，保证跨进程接触力计算正确
- **粒子迁移**：粒子跨过子域边界时自动迁移到所属进程
- **全局归约**：颗粒数等全局量通过 `MPI_Allreduce` 汇总

### 结果输出（`src/output.h`）
- **VTK/ParaView 可视化**：每个 rank 写出本地点为 `frame_XXXXX_rr.vtu`，rank 0 聚合生成 `.pvd` 索引文件；粒子数据附带 id / 材料号 / 半径 / 速度 / 角速度 / 力模长
- ParaView 中打开 `.pvd` 后用 Glyph 过滤器（Glyph Type=Sphere, Scale Array=radius）即可按真实半径渲染颗粒并按速度着色

## 编译与运行

依赖：C++17 编译器、CMake ≥ 3.20、MPI 实现（如 OpenMPI / MPICH）。

```bash
# 编译
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build

# 运行：第一个参数为进程网格列数 pcols，第二个为行数 prows，二者乘积须等于进程数
mpirun -np 8 ./build/dem2d 4 2
```

默认场景：[0,0]×[100,10] 计算域、四壁容器，约 7200 个等半径颗粒随机布置后自由下落堆积，共 15000 步，每 500 步输出一帧。

## 代码结构

```
src/
├── main.cpp       # 入口：MPI 初始化、场景搭建、主循环、结果输出
├── types.h        # 基础类型（Vec2、Material、Particle 等）
├── particle.cpp/h # 颗粒与材料定义
├── grid.cpp/h     # 均匀网格邻域搜索
├── contact.cpp/h  # 颗粒-颗粒 Hertz-Mindlin 接触模型
├── boundary.cpp/h # 墙体与颗粒-墙接触模型
├── world.cpp/h    # 模拟世界：时间步进、ghost 交换、粒子迁移
└── output.cpp/h   # VTU 帧写出器
```
