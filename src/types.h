#ifndef TYPES_H
#define TYPES_H

#include <cstdint>
#include <math.h>

using real_t = double;
using id_t = uint32_t;
using mid_t = uint32_t;

const real_t pi = 3.14159265358979323846;

/**
 * @brief 二维向量，提供基本向量运算
 */
struct Vec2{
    real_t x, y;

    Vec2 operator+(const Vec2& v) const { return Vec2{x + v.x, y + v.y}; }
    Vec2 operator-(const Vec2& v) const { return Vec2{x - v.x, y - v.y}; }
    Vec2 operator*(real_t s) const { return Vec2{x * s, y * s}; }

    Vec2& operator+=(const Vec2& v) { x += v.x; y += v.y; return *this; }
    Vec2& operator-=(const Vec2& v) { x -= v.x; y -= v.y; return *this; }
    Vec2& operator*=(real_t s) { x *= s; y *= s; return *this; }

    Vec2 operator-() const { return Vec2{-x, -y}; }

    real_t& operator[](size_t i) { return (i == 0) ? x : y;}
    const real_t& operator[](size_t i) const { return (i == 0) ? x : y;}

    /**
     * @brief 与另一向量的点积
     * @param v 另一个向量
     * @return 点积结果
     */
    real_t dot(const Vec2& v) const { return x * v.x + y * v.y; }

    /**
     * @brief 向量长度（欧几里得范数）
     * @return 长度
     */
    real_t length() const { return std::sqrt(x * x + y * y); }

    /**
     * @brief 向量长度的平方，避免开方
     * @return 长度平方
     */
    real_t length_squared() const { return x * x + y * y; }

    /**
     * @brief 归一化为单位向量，零向量返回零向量
     * @return 单位向量
     */
    Vec2 normalize() const {
        auto len = length();
        return (len > real_t(0)) ? Vec2{x / len, y / len} : Vec2{ 0, 0};
    }
};

inline Vec2 operator*(real_t s, const Vec2& v) { return Vec2{v.x * s, v.y * s}; }

#endif
