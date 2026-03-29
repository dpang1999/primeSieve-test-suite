#pragma once
#include "IField.h"
#include <iostream>
#include <cstdint>
#include <functional>

inline uint64_t mod_inverse(uint64_t a, uint64_t m) {
    if (m == 1) return 0;
    int64_t a_i64 = a;
    int64_t m_i64 = m;
    int64_t m0 = m_i64;
    int64_t y = 0, x = 1;
    while (a_i64 > 1) {
        int64_t q = a_i64 / m_i64;
        int64_t t = m_i64;
        m_i64 = a_i64 % m_i64;
        a_i64 = t;
        t = y;
        y = x - q * y;
        x = t;
    }
    if (x < 0) x += m0;
    return (uint64_t)(x % m0);
}

struct IntModP {
    uint64_t val;
    static uint64_t modulus;
    IntModP(uint64_t v = 0) : val(v % modulus) {}
    IntModP a(const IntModP& o) const { return IntModP((val + o.val) % modulus); }
    IntModP s(const IntModP& o) const { return IntModP((val + modulus - (o.val % modulus)) % modulus); }
    IntModP m(const IntModP& o) const { return IntModP((val * o.val) % modulus); }
    IntModP d(const IntModP& o) const { return IntModP((val * mod_inverse(o.val, modulus)) % modulus); }
    IntModP zero() const { return IntModP(0); }
    double coerce_to_f64() const { return (double)val; }
    bool operator==(const IntModP& o) const { return val == o.val; }
    struct Hash { size_t operator()(const IntModP& f) const { return std::hash<uint64_t>()(f.val); } };
    friend std::ostream& operator<<(std::ostream& os, const IntModP& coeff) { return os << coeff.val; }
};
inline uint64_t IntModP::modulus = 7;
