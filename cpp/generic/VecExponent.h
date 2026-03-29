#pragma once
#include "IExponent.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

struct VecExponent {
    std::vector<int> exps;
    VecExponent() {}
    VecExponent(std::vector<int> v) : exps(v) {}

    int lex_compare(const VecExponent& other) const {
        for (size_t i = 0; i < exps.size() && i < other.exps.size(); ++i) {
            if (exps[i] < other.exps[i]) return -1;
            if (exps[i] > other.exps[i]) return 1;
        }
        if (exps.size() < other.exps.size()) return -1;
        if (exps.size() > other.exps.size()) return 1;
        return 0;
    }

    uint64_t degree() const {
        uint64_t deg = 0;
        for (int e : exps) deg += e;
        return deg;
    }

    bool can_reduce(const VecExponent& other) const {
        for (size_t i = 0; i < exps.size(); ++i) {
            int other_val = i < other.exps.size() ? other.exps[i] : 0;
            if (exps[i] < other_val) return false;
        }
        return true;
    }

    VecExponent lcm(const VecExponent& other) const {
        size_t len = std::max(exps.size(), other.exps.size());
        std::vector<int> res(len);
        for (size_t i = 0; i < len; ++i) {
            int e1 = i < exps.size() ? exps[i] : 0;
            int e2 = i < other.exps.size() ? other.exps[i] : 0;
            res[i] = std::max(e1, e2);
        }
        return VecExponent(res);
    }

    VecExponent add(const VecExponent& other) const {
        size_t len = std::max(exps.size(), other.exps.size());
        std::vector<int> res(len);
        for (size_t i = 0; i < len; ++i) {
            int e1 = i < exps.size() ? exps[i] : 0;
            int e2 = i < other.exps.size() ? other.exps[i] : 0;
            res[i] = e1 + e2;
        }
        return VecExponent(res);
    }

    VecExponent sub(const VecExponent& other) const {
        size_t len = std::max(exps.size(), other.exps.size());
        std::vector<int> res(len);
        for (size_t i = 0; i < len; ++i) {
            int e1 = i < exps.size() ? exps[i] : 0;
            int e2 = i < other.exps.size() ? other.exps[i] : 0;
            res[i] = e1 - e2;
        }
        return VecExponent(res);
    }

    bool operator==(const VecExponent& other) const {
        if (exps.size() != other.exps.size()) return false;
        for (size_t i = 0; i < exps.size(); ++i) {
            if (exps[i] != other.exps[i]) return false;
        }
        return true;
    }

    struct Hash { 
        size_t operator()(const VecExponent& e) const { 
            size_t hash = 0;
            for (int v : e.exps) hash ^= std::hash<int>()(v) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            return hash;
        } 
    };

    friend std::ostream& operator<<(std::ostream& os, const VecExponent& e) {
        os << "[";
        for (size_t i=0; i<e.exps.size(); i++) {
            os << e.exps[i] << (i+1==e.exps.size()?"":",");
        }
        return os << "]";
    }
};
