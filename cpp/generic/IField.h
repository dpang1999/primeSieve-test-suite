#pragma once
#include <cstddef>
// In C++17, templates use duck typing. This serves as documentation for the IField concept.
// Any type C used as a coefficient must implement:
// C a(const C& o) const;
// C s(const C& o) const;
// C m(const C& o) const;
// C d(const C& o) const;
// C zero() const;
// double coerce_to_f64() const;
// bool operator==(const C& o) const;
// struct Hash { size_t operator()(const C& f) const; };
