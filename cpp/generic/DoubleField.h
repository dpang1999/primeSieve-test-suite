#pragma once
#include "IField.h"
#include <iostream>
#include <functional>

struct DoubleField {
    double val;
    DoubleField(double v = 0.0) : val(v) {}
    DoubleField a(const DoubleField& o) const { return DoubleField(val + o.val); }
    DoubleField s(const DoubleField& o) const { return DoubleField(val - o.val); }
    DoubleField m(const DoubleField& o) const { return DoubleField(val * o.val); }
    DoubleField d(const DoubleField& o) const { return DoubleField(val / o.val); }
    DoubleField zero() const { return DoubleField(0.0); }
    double coerce_to_f64() const { return val; }
    bool operator==(const DoubleField& o) const { return val == o.val; }
    struct Hash { size_t operator()(const DoubleField& f) const { return std::hash<double>()(f.val); } };
    friend std::ostream& operator<<(std::ostream& os, const DoubleField& coeff) { return os << coeff.val; }
};
