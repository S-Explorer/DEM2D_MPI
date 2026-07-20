#include "contact.h"
#include <cmath>

// 粒子-粒子接触力：Hertz 法向 + Mindlin 切向 + Tsuji 阻尼。

/**
 * @brief 构造接触力模型，绑定材料表
 * @param materials 材料数组
 * @param num_materials 材料数量
 */
ContactModel::ContactModel(const Material* materials, int num_materials)
    : materials(materials), num_materials(num_materials) {}

/**
 * @brief 计算粒子 a 与 b 的接触力，累加到各自的力与力矩
 * @param a 粒子 a
 * @param b 粒子 b
 * @param c 该粒子对的接触历史
 * @param dt 时间步长
 */
void ContactModel::compute_force(Particle& a, Particle& b, Contact& c, real_t dt) {
    const Material& ma = materials[a.mid];
    const Material& mb = materials[b.mid];

    Vec2 dx = {b.pos.x - a.pos.x, b.pos.y - a.pos.y};
    real_t dist = dx.length();

    real_t delta_n = a.radius + b.radius - dist;
    if (delta_n <= real_t(0)) {
        c.active = false;
        return;
    }
    c.active = true;

    // 接触局部坐标：n 由 a 指向 b，t 为 n 逆时针 90°
    real_t inv_dist = real_t(1) / dist;
    Vec2 n = {dx.x * inv_dist, dx.y * inv_dist};
    Vec2 t = {-n.y, n.x};

    // 等效参数 R*、m*
    real_t R_star = (a.radius * b.radius) / (a.radius + b.radius);
    real_t m_star = (a.mass * b.mass) / (a.mass + b.mass);

    // E* = 1 / [ (1-nu1^2)/E1 + (1-nu2^2)/E2 ]
    real_t e1 = (real_t(1) - ma.nu * ma.nu) / ma.E;
    real_t e2 = (real_t(1) - mb.nu * mb.nu) / mb.E;
    real_t E_star = real_t(1) / (e1 + e2);

    // G* = 1 / [ (2-nu1)/G1 + (2-nu2)/G2 ]
    real_t G_star = real_t(1) / ((real_t(2) - ma.nu) / ma.G + (real_t(2) - mb.nu) / mb.G);

    // 接触点相对速度（2D）：v_rel = (v_b - v_a) - (ω_a·R_a + ω_b·R_b)·t
    Vec2 v_rel_center = {b.vel.x - a.vel.x, b.vel.y - a.vel.y};
    real_t omega_term = a.omega * a.radius + b.omega * b.radius;
    real_t v_n   = v_rel_center.dot(n);
    real_t v_t   = v_rel_center.dot(t) - omega_term;

    // 恢复系数与摩擦系数取两者较小者，偏保守
    real_t e   = (ma.e < mb.e) ? ma.e : mb.e;
    real_t mu  = (ma.mu < mb.mu) ? ma.mu : mb.mu;

    // --- 法向力（Hertz）---
    // Fn_elastic = (4/3)·E*·sqrt(R*)·δn^(3/2)
    real_t sqrt_R_star = std::sqrt(R_star);
    real_t kn = (real_t(4) / real_t(3)) * E_star * sqrt_R_star;
    real_t sqrt_delta = std::sqrt(delta_n);
    real_t Fn_elastic = kn * delta_n * sqrt_delta;

    // Tsuji 阻尼：cn = -2·ln(e)·sqrt(5/6·m*·kn)·δn^(1/4)
    real_t cn = -real_t(2) * std::log(e) * std::sqrt(real_t(5) / real_t(6) * m_star * kn)
                * std::pow(delta_n, real_t(0.25));
    // δ̇n = -v_n（接近时 v_n < 0），故阻尼项写作 -cn·v_n
    real_t Fn = Fn_elastic - cn * v_n;

    // --- 切向位移更新 ---
    // dδt = v_t·dt，并投影到当前接触面
    c.t_disp.x += v_t * t.x * dt;
    c.t_disp.y += v_t * t.y * dt;

    real_t proj = c.t_disp.dot(n);
    c.t_disp.x -= proj * n.x;
    c.t_disp.y -= proj * n.y;

    // --- 切向力（Mindlin）---
    // kt = 8·G*·sqrt(R*·δn)，Ft = -kt·δt - ct·v_t
    real_t kt = real_t(8) * G_star * std::sqrt(R_star * delta_n);
    real_t ct = -real_t(2) * std::log(e) * std::sqrt(real_t(5) / real_t(6) * m_star * kt);

    Vec2 Ft = {c.t_disp.x * (-kt) - ct * v_t * t.x,
               c.t_disp.y * (-kt) - ct * v_t * t.y};

    // Coulomb 限幅：|Ft| <= μ·|Fn|
    real_t Ft_mag   = Ft.length();
    real_t Ft_max   = mu * std::abs(Fn);
    if (Ft_max > 0 && Ft_mag > Ft_max) {
        real_t ratio = Ft_max / Ft_mag;
        Ft.x     *= ratio;
        Ft.y     *= ratio;
        c.t_disp.x = -Ft.x / kt;   // 截断后回写 δt，保持与 Ft 一致
        c.t_disp.y = -Ft.y / kt;
    }

    // --- 施加力与力矩 ---
    // 作用力与反作用力：F = Fn·n + Ft
    a.force.x += -Fn * n.x - Ft.x;
    a.force.y += -Fn * n.y - Ft.y;
    b.force.x +=  Fn * n.x + Ft.x;
    b.force.y +=  Fn * n.y + Ft.y;

    // 力矩：τ = r × F，r_a = a.radius·n，r_b = -b.radius·n
    real_t cross_n_Ft = n.x * Ft.y - n.y * Ft.x;
    a.torque -= a.radius * cross_n_Ft;
    b.torque -= b.radius * cross_n_Ft;
}
