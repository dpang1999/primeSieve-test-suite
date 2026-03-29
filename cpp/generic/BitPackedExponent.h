#pragma once
#include "IExponent.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

struct BitPackedExponent {
    uint64_t packed;

    BitPackedExponent() : packed(0) {}
    BitPackedExponent(uint64_t p) : packed(p) {}
    BitPackedExponent(std::vector<int> exps) {
        packed = 0;
        while(exps.size() < 6) exps.push_back(0);
        for (int i = 0; i < 6; i++) {
            uint64_t shift = 40 - 8 * i;
            packed |= ((uint64_t)exps[i] & 0xFF) << shift;
        }
        uint64_t deg = 0;
        for (int i = 0; i < 6; i++) deg += exps[i];
        packed |= (deg & 0xFFFF) << 48;
    }

    int lex_compare(const BitPackedExponent& other) const {
        uint64_t e1 = packed & 0x0000FFFFFFFFFFFFULL;
        uint64_t e2 = other.packed & 0x0000FFFFFFFFFFFFULL;
        if (e1 < e2) return -1;
        if (e1 > e2) return 1;
        return 0;
    }

    uint64_t degree() const { return (packed >> 48) & 0xFFFF; }

    bool can_reduce(const BitPackedExponent& other) const {
        for (int i = 0; i <= 40; i += 8) {
            uint64_t self_exp = (packed >> i) & 0xFF;
            uint64_t divisor_exp = (other.packed >> i) & 0xFF;
            if (self_exp < divisor_exp) return false;
        }
        return true;
    }

    BitPackedExponent lcm(const BitPackedExponent& other) const {
        uint64_t self_exponents = packed & 0x0000FFFFFFFFFFFFULL;
        uint64_t other_exponents = other.packed & 0x0000FFFFFFFFFFFFULL;
        uint64_t lcm_exponents = 0;
        uint64_t deg = 0;
        for (int i = 0; i <= 40; i += 8) {
            uint64_t self_exp = (self_exponents >> i) & 0xFF;
            uint64_t other_exp = (other_exponents >> i) & 0xFF;
            uint64_t lcm_exp = std::max(self_exp, other_exp);
            lcm_exponents |= (lcm_exp << i);
            deg += lcm_exp;
        }
        lcm_exponents |= (deg & 0xFFFF) << 48;
        return BitPackedExponent(lcm_exponents);
    }

    BitPackedExponent add(const BitPackedExponent& other) const { return BitPackedExponent(packed + other.packed); }
    BitPackedExponent sub(const BitPackedExponent& other) const { return BitPackedExponent(packed - other.packed); }

    bool operator==(const BitPackedExponent& other) const { return packed == other.packed; }

    struct Hash { size_t operator()(const BitPackedExponent& e) const { return std::hash<uint64_t>()(e.packed); } };

    friend std::ostream& operator<<(std::ostream& os, const BitPackedExponent& e) {
        return os << std::hex << "0x" << e.packed << std::dec;
    }
};
