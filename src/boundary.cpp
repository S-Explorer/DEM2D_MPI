#include "boundary.h"
#include <cmath>

// 粒子-墙接触力：Hertz 法向 + Mindlin 切向 + Tsuji 阻尼。
// 与 contact.cpp 的粒子-粒子模型同一套公式；墙视作无限大质量刚体
// （E_wall, G_wall -> ∞，曲率半径 -> ∞），故等效参数只含粒子一侧：
//   R* = R_p,   m* = m_p,
//   E* = E_p / (1 - nu_p^2),   G* = G_p / (2 - nu_p)
//
// 约定（与粒子-粒子保持一致，便于对照）：
//   n  墙外法线（墙 -> 粒子），单位向量
//   t  = (-n.y, n.x)        （n 逆时针 90°）
//   δn = R_p - s,  s 为粒子中心到墙的有符号距离（沿 n，合法侧为正）
//   v_n = v·n              （粒子朝墙运动时 < 0）
//   v_t = v·t - ω·R        （接触点相对墙的切向速度，墙固定）
// 法向 Tsuji 阻尼 Fn = kn·δn^(3/2) + cn·δ̇n，其中 δ̇n = -v_n（重叠量随接近而增大），
// 因此代码中写作  Fn = Fn_elastic - cn·v_n。

/**
 * @brief 构造墙接触力模型，绑定材料表
 * @param materials 材料数组
 * @param num_materials 材料数量
 */
WallContactModel::WallContactModel(const Material* materials, int num_materials)
    : materials(materials), num_materials(num_materials) {}

/**
 * @brief 计算墙对粒子的接触力，累加到粒子的力与力矩
 * @param p 粒子
 * @param w 墙
 * @param wc 该 (粒子,墙) 对的切向历史，跨时间步复用
 * @param dt 时间步长
 */
void WallContactModel::compute_force(Particle& p, const Wall& w, WallContact& wc, real_t dt) {
    const Material& m = materials[p.mid];

    // 外法线 n 与粒子中心到墙的有符号距离 s
    Vec2 n;
    real_t s;
    if (w.axis == 0) {
        n = {static_cast<real_t>(w.dir), 0};
        s = (p.pos.x - w.pos) * static_cast<real_t>(w.dir);
    } else {
        n = {0, static_cast<real_t>(w.dir)};
        s = (p.pos.y - w.pos) * static_cast<real_t>(w.dir);
    }

    real_t delta_n = p.radius - s;
    if (delta_n <= real_t(0)) {
        wc.active = false;
        wc.t_disp = {0, 0};   // 脱离接触：清切向历史，避免再次接触时残留旧位移
        return;
    }
    wc.active = true;

    Vec2 t = {-n.y, n.x};

    // 等效参数（墙为刚体，只含粒子一侧）
    real_t R_star = p.radius;
    real_t m_star = p.mass;
    real_t E_star = m.E / (real_t(1) - m.nu * m.nu);
    real_t G_star = m.G / (real_t(2) - m.nu);
    real_t e      = m.e;
    real_t mu     = m.mu;

    // 接触点速度（墙固定）：法向只含平动，切向含转动贡献
    real_t v_n = p.vel.dot(n);
    real_t v_t = p.vel.dot(t) - p.omega * p.radius;

    // --- 法向力（Hertz + Tsuji 阻尼）---
    real_t sqrt_R   = std::sqrt(R_star);
    real_t kn       = (real_t(4) / real_t(3)) * E_star * sqrt_R;
    real_t Fn_elast = kn * delta_n * std::sqrt(delta_n);
    real_t cn       = -real_t(2) * std::log(e)
                    * std::sqrt(real_t(5) / real_t(6) * m_star * kn)
                    * std::pow(delta_n, real_t(0.25));
    real_t Fn = Fn_elast - cn * v_n;     // δ̇n = -v_n

    // --- 切向位移更新并投影到当前接触面 ---
    // 墙法线固定，投影理论上为零；保留以与粒子-粒子代码结构一致，便于教学对照。
    wc.t_disp.x += v_t * t.x * dt;
    wc.t_disp.y += v_t * t.y * dt;
    real_t proj = wc.t_disp.dot(n);
    wc.t_disp.x -= proj * n.x;
    wc.t_disp.y -= proj * n.y;

    // --- 切向力（Mindlin + 阻尼）---
    real_t kt = real_t(8) * G_star * std::sqrt(R_star * delta_n);
    real_t ct = -real_t(2) * std::log(e)
              * std::sqrt(real_t(5) / real_t(6) * m_star * kt);
    Vec2 Ft = { wc.t_disp.x * (-kt) - ct * v_t * t.x,
                wc.t_disp.y * (-kt) - ct * v_t * t.y };

    // Coulomb 限幅：|Ft| <= mu·|Fn|
    real_t Ft_mag = Ft.length();
    real_t Ft_max = mu * std::abs(Fn);
    if (Ft_max > 0 && Ft_mag > Ft_max) {
        real_t ratio = Ft_max / Ft_mag;
        Ft.x *= ratio;
        Ft.y *= ratio;
        wc.t_disp.x = -Ft.x / kt;   // 与截断后的 Ft 保持一致
        wc.t_disp.y = -Ft.y / kt;
    }

    // --- 施加力与力矩 ---
    // 接触点 r_cp = -R·n（靠墙一侧），F = Fn·n + Ft
    // τ = r_cp × F = -R·(n × Ft)
    p.force.x += Fn * n.x + Ft.x;
    p.force.y += Fn * n.y + Ft.y;
    real_t cross = n.x * Ft.y - n.y * Ft.x;
    p.torque   -= p.radius * cross;
}
