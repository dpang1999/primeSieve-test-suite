#pragma once
#include <cstdint>
#include <cstddef>
// In C++17, templates use duck typing. This serves as documentation for the IExponent concept.
// Any type E used as an exponent must implement:
// int lex_compare(const E& other) const;
// uint64_t degree() const;
// bool can_reduce(const E& other) const;
// E lcm(const E& other) const;
// E add(const E& other) const;
// E sub(const E& other) const;
// bool operator==(const E& other) const;
// struct Hash { size_t operator()(const E& e) const; };
