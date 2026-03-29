#pragma once
#include "IField.h"
#include <iostream>
#include <functional>

struct SingleField {
    float val;
    SingleField(float v = 0.0f) : val(v) {}
    SingleField a(const SingleField& o) const { return SingleField(val + o.val); }
    SingleField s(const SingleField& o) const { return SingleField(val - o.val); }
    SingleField m(const SingleField& o) const { return SingleField(val * o.val); }
    SingleField d(const SingleField& o) const { return SingleField(val / o.val); }
    SingleField zero() const { return SingleField(0.0f); }
    double coerce_to_f64() const { return (double)val; }
    bool operator==(const SingleField& o) const { return val == o.val; }
    struct Hash { size_t operator()(const SingleField& f) const { return std::hash<float>()(f.val); } };
    friend std::ostream& operator<<(std::ostream& os, const SingleField& coeff) { return os << coeff.val; }
};
