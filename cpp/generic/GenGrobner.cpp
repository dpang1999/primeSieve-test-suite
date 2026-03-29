#include <iostream>
#include <vector>
#include <chrono>
#include <unordered_set>
#include "IntModP.h"
#include "VecExponent.h"
#include "BitPackedExponent.h"
#include "DoubleField.h"
#include "SingleField.h"

enum TermOrder { Lex, GrLex, RevLex };
extern TermOrder TERM_ORDER;


template<typename C, typename E>
struct Term {
    C coefficient;
    E exponents;

    Term(C c, E e) : coefficient(c), exponents(e) {}

    int compare(const Term& other) const {
        if (TERM_ORDER == Lex) {
            return exponents.lex_compare(other.exponents);
        } else if (TERM_ORDER == GrLex) {
            uint64_t d1 = exponents.degree();
            uint64_t d2 = other.exponents.degree();
            if (d1 != d2) return d1 < d2 ? -1 : 1;
            return exponents.lex_compare(other.exponents);
        } else {
            uint64_t d1 = exponents.degree();
            uint64_t d2 = other.exponents.degree();
            if (d1 != d2) return d1 < d2 ? -1 : 1;
            return -exponents.lex_compare(other.exponents); 
        }
    }

    bool can_reduce(const Term& divisor_leading) const {
        return exponents.can_reduce(divisor_leading.exponents);
    }

    E lcm(const Term& other) const {
        return exponents.lcm(other.exponents);
    }

    bool operator==(const Term& o) const { return coefficient == o.coefficient && exponents == o.exponents; }
};

#include <vector>
#include <cmath>
#include <algorithm>

template<typename C, typename E>
struct Polynomial {
    std::vector<Term<C, E>> terms;

    Polynomial() {}
    Polynomial(std::vector<Term<C, E>> _terms) {
        for (auto& t : _terms) {
            if (std::abs(t.coefficient.coerce_to_f64()) > 0.0) {
                terms.push_back(t);
            }
        }
        std::sort(terms.begin(), terms.end(), [](const Term<C, E>& a, const Term<C, E>& b) {
            return a.compare(b) > 0;
        });
    }

    Polynomial addPoly(const Polynomial& other) const {
        std::vector<Term<C, E>> result = terms;
        for (const auto& term : other.terms) {
            bool found = false;
            for (auto& res_term : result) {
                if (res_term.exponents == term.exponents) {
                    res_term.coefficient = res_term.coefficient.a(term.coefficient);
                    found = true;
                    break;
                }
            }
            if (!found) result.push_back(term);
        }
        return Polynomial(result);
    }

    Polynomial add_poly(const Polynomial& other) const {
        std::vector<Term<C, E>> result = terms;
        for (const auto& term2 : other.terms) {
            bool found = false;
            for (auto& res_term : result) {
                if (res_term.exponents == term2.exponents) {
                    res_term.coefficient = res_term.coefficient.a(term2.coefficient);
                    found = true;
                    break;
                }
            }
            if (!found) {
                result.push_back(term2);
            }
        }
        return Polynomial(result);
    }

    Polynomial subtract(const Polynomial& other) const {
        std::vector<Term<C, E>> result = terms;
        for (const auto& term : other.terms) {
            bool found = false;
            for (auto& res_term : result) {
                if (res_term.exponents == term.exponents) {
                    res_term.coefficient = res_term.coefficient.s(term.coefficient);
                    found = true;
                    break;
                }
            }
            if (!found) {
                result.push_back(Term<C, E>(term.coefficient.zero().s(term.coefficient), term.exponents));
            }
        }
        return Polynomial(result);
    }
    Polynomial make_monic() const {
        if (terms.empty()) { return *this; }
        C lead_coeff = terms[0].coefficient;
        std::vector<Term<C, E>> new_terms;
        new_terms.reserve(terms.size());
        for (const auto& t : terms) {
            new_terms.push_back(Term<C, E>{t.coefficient.d(lead_coeff), t.exponents});
        }
        return Polynomial(new_terms);
    }


    Polynomial reduce_poly(const std::vector<Polynomial>& divisors) const {
        Polynomial result = *this;
        std::vector<Term<C, E>> remainder;
        
        while (!result.terms.empty()) {
            bool reduced = false;
            Term<C, E> leading_term = result.terms[0];
            
            for (const auto& divisor : divisors) {
                if (divisor.terms.empty()) continue;
                Term<C, E> divisor_leading_term = divisor.terms[0];
                
                if (leading_term.can_reduce(divisor_leading_term)) {
                    C coeff = leading_term.coefficient.d(divisor_leading_term.coefficient);
                    E exps = leading_term.exponents.sub(divisor_leading_term.exponents);
                    Term<C, E> reduction_term(coeff, exps);
                    Polynomial scaled_divisor = divisor.multiply_by_term(reduction_term);
                    result = result.subtract(scaled_divisor);
                    reduced = true;
                    break;
                }
            }
            
            if (!reduced) {
                remainder.push_back(result.terms[0]);
                result.terms.erase(result.terms.begin());
            }
        }
        
        for (const auto& t : remainder) {
            result.terms.push_back(t);
        }
        return Polynomial(result.terms);
    }
    Polynomial multiply_by_term(const Term<C, E>& term) const {
        std::vector<Term<C, E>> new_terms;
        new_terms.reserve(terms.size());
        for (const auto& t : terms) {
            new_terms.push_back(Term<C, E>{t.coefficient.m(term.coefficient), t.exponents.add(term.exponents)});
        }
        return Polynomial(new_terms);
    }

    static Polynomial s_polynomial(const Polynomial& p1, const Polynomial& p2) {
        Term<C, E> leading_term_p1 = p1.terms[0];
        Term<C, E> leading_term_p2 = p2.terms[0];

        E lcm_exponents = leading_term_p1.lcm(leading_term_p2);
        E scale_factor_p1 = lcm_exponents.sub(leading_term_p1.exponents);
        E scale_factor_p2 = lcm_exponents.sub(leading_term_p2.exponents);

        Polynomial scaled_p1 = p1.multiply_by_term(Term<C, E>(leading_term_p2.coefficient, scale_factor_p1));
        Polynomial scaled_p2 = p2.multiply_by_term(Term<C, E>(leading_term_p1.coefficient, scale_factor_p2));
        
        return scaled_p1.subtract(scaled_p2);
    }

    bool operator==(const Polynomial& o) const {
        if (terms.size() != o.terms.size()) return false;
        for (size_t i = 0; i < terms.size(); i++) {
            if (!(terms[i] == o.terms[i])) return false;
        }
        return true;
    }
};

template<typename C, typename E>
struct PolyHash {
    size_t operator()(const Polynomial<C, E>& p) const {
        size_t hash = 0;
        for (const auto& t : p.terms) {
            hash ^= typename C::Hash()(t.coefficient) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            hash ^= typename E::Hash()(t.exponents) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

TermOrder TERM_ORDER = Lex;

template<typename C, typename E>
std::vector<Polynomial<C, E>> naive_grobner_basis(std::vector<Polynomial<C, E>> basis) {
    std::vector<std::pair<size_t, size_t>> pairs;
    for (size_t i = 0; i < basis.size(); ++i) {
        for (size_t j = i + 1; j < basis.size(); ++j) {
            pairs.push_back({i, j});
        }
    }

    std::unordered_set<Polynomial<C, E>, PolyHash<C, E>> basis_set(basis.begin(), basis.end());

    while (!pairs.empty()) {
        auto [i, j] = pairs.front();
        pairs.erase(pairs.begin());

        Polynomial<C, E> s_poly = Polynomial<C, E>::s_polynomial(basis[i], basis[j]);
        Polynomial<C, E> remainder = s_poly.reduce_poly(basis);

        if (!remainder.terms.empty() && basis_set.insert(remainder).second) {
            for (size_t k = 0; k < basis.size(); ++k) {
                pairs.push_back({k, basis.size()});
            }
            basis.push_back(remainder);
        }
    }

    std::vector<Polynomial<C, E>> reduced_basis;
    for (const auto& poly : basis) {
        std::vector<Polynomial<C, E>> basis_excluding_self;
        for (const auto& p : basis) {
            if (!(p == poly)) {
                basis_excluding_self.push_back(p);
            }
        }
        Polynomial<C, E> reduced = poly.reduce_poly(basis_excluding_self);
        if (!reduced.terms.empty()) {
            reduced = reduced.make_monic();
            bool contains = false;
            for (const auto& r : reduced_basis) {
                if (r == reduced) { contains = true; break; }
            }
            if (!contains) {
                reduced_basis.push_back(reduced);
            }
        }
    }

    return reduced_basis;
}


template <typename C, typename E>
void run_benchmark(int n) {
    std::vector<Polynomial<C, E>> sys;
    if constexpr (std::is_same_v<E, VecExponent>) {
        if (n == 4) {
                // Cyclic 4 Vec
                Polynomial<C, E> p1({
                    Term<C, E>(C(1), E({ 1, 0, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 0, 1 })),
                });
                sys.push_back(p1);
                Polynomial<C, E> p2({
                    Term<C, E>(C(1), E({ 1, 1, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 0, 0, 1 })),
                });
                sys.push_back(p2);
                Polynomial<C, E> p3({
                    Term<C, E>(C(1), E({ 1, 1, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 0, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 1, 0, 1 })),
                });
                sys.push_back(p3);
                Polynomial<C, E> p4({
                    Term<C, E>(C(1), E({ 1, 1, 1, 1 })),
                    Term<C, E>(C(C::modulus - 1), E({ 0, 0, 0, 0 }))
                });
                sys.push_back(p4);

        } else if (n == 5) {
                // Cyclic 5 Vec
                Polynomial<C, E> p1({
                    Term<C, E>(C(1), E({ 1, 0, 0, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 0, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 1, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 0, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 0, 0, 1 })),
                });
                sys.push_back(p1);
                Polynomial<C, E> p2({
                    Term<C, E>(C(1), E({ 1, 1, 0, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 1, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 1, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 0, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 0, 0, 0, 1 })),
                });
                sys.push_back(p2);
                Polynomial<C, E> p3({
                    Term<C, E>(C(1), E({ 1, 1, 1, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 1, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 1, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 0, 0, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 1, 0, 0, 1 })),
                });
                sys.push_back(p3);
                Polynomial<C, E> p4({
                    Term<C, E>(C(1), E({ 1, 1, 1, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 1, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 0, 1, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 1, 0, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 1, 1, 0, 1 })),
                });
                sys.push_back(p4);
                Polynomial<C, E> p5({
                    Term<C, E>(C(1), E({ 1, 1, 1, 1, 1 })),
                    Term<C, E>(C(C::modulus - 1), E({ 0, 0, 0, 0, 0 }))
                });
                sys.push_back(p5);

        } else if (n == 6) {
                // Cyclic 6 Vec
                Polynomial<C, E> p1({
                    Term<C, E>(C(1), E({ 1, 0, 0, 0, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 0, 0, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 1, 0, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 0, 1, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 0, 0, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 0, 0, 0, 1 })),
                });
                sys.push_back(p1);
                Polynomial<C, E> p2({
                    Term<C, E>(C(1), E({ 1, 1, 0, 0, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 1, 0, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 1, 1, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 0, 1, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 0, 0, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 0, 0, 0, 0, 1 })),
                });
                sys.push_back(p2);
                Polynomial<C, E> p3({
                    Term<C, E>(C(1), E({ 1, 1, 1, 0, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 1, 1, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 1, 1, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 0, 1, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 0, 0, 0, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 1, 0, 0, 0, 1 })),
                });
                sys.push_back(p3);
                Polynomial<C, E> p4({
                    Term<C, E>(C(1), E({ 1, 1, 1, 1, 0, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 1, 1, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 0, 1, 1, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 0, 0, 1, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 1, 0, 0, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 1, 1, 0, 0, 1 })),
                });
                sys.push_back(p4);
                Polynomial<C, E> p5({
                    Term<C, E>(C(1), E({ 1, 1, 1, 1, 1, 0 })),
                    Term<C, E>(C(1), E({ 0, 1, 1, 1, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 0, 1, 1, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 1, 0, 1, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 1, 1, 0, 1, 1 })),
                    Term<C, E>(C(1), E({ 1, 1, 1, 1, 0, 1 })),
                });
                sys.push_back(p5);
                Polynomial<C, E> p6({
                    Term<C, E>(C(1), E({ 1, 1, 1, 1, 1, 1 })),
                    Term<C, E>(C(C::modulus - 1), E({ 0, 0, 0, 0, 0, 0 }))
                });
                sys.push_back(p6);

        } else {
            std::cout << "Invalid n for hardcoded benchmark" << std::endl; return;
        }
    } else {
        if (n == 4) {
                // Cyclic 4 BitPacked
                Polynomial<C, E> p1({
                    Term<C, E>(C(1), E(0x0001010000000000ULL)),
                    Term<C, E>(C(1), E(0x0001000100000000ULL)),
                    Term<C, E>(C(1), E(0x0001000001000000ULL)),
                    Term<C, E>(C(1), E(0x0001000000010000ULL)),
                });
                sys.push_back(p1);
                Polynomial<C, E> p2({
                    Term<C, E>(C(1), E(0x0002010100000000ULL)),
                    Term<C, E>(C(1), E(0x0002000101000000ULL)),
                    Term<C, E>(C(1), E(0x0002000001010000ULL)),
                    Term<C, E>(C(1), E(0x0002010000010000ULL)),
                });
                sys.push_back(p2);
                Polynomial<C, E> p3({
                    Term<C, E>(C(1), E(0x0003010101000000ULL)),
                    Term<C, E>(C(1), E(0x0003000101010000ULL)),
                    Term<C, E>(C(1), E(0x0003010001010000ULL)),
                    Term<C, E>(C(1), E(0x0003010100010000ULL)),
                });
                sys.push_back(p3);
                Polynomial<C, E> p4({
                    Term<C, E>(C(1), E(0x0004010101010000ULL)),
                    Term<C, E>(C(C::modulus - 1), E(0x0ULL))
                });
                sys.push_back(p4);

        } else if (n == 5) {
                // Cyclic 5 BitPacked
                Polynomial<C, E> p1({
                    Term<C, E>(C(1), E(0x0001010000000000ULL)),
                    Term<C, E>(C(1), E(0x0001000100000000ULL)),
                    Term<C, E>(C(1), E(0x0001000001000000ULL)),
                    Term<C, E>(C(1), E(0x0001000000010000ULL)),
                    Term<C, E>(C(1), E(0x0001000000000100ULL)),
                });
                sys.push_back(p1);
                Polynomial<C, E> p2({
                    Term<C, E>(C(1), E(0x0002010100000000ULL)),
                    Term<C, E>(C(1), E(0x0002000101000000ULL)),
                    Term<C, E>(C(1), E(0x0002000001010000ULL)),
                    Term<C, E>(C(1), E(0x0002000000010100ULL)),
                    Term<C, E>(C(1), E(0x0002010000000100ULL)),
                });
                sys.push_back(p2);
                Polynomial<C, E> p3({
                    Term<C, E>(C(1), E(0x0003010101000000ULL)),
                    Term<C, E>(C(1), E(0x0003000101010000ULL)),
                    Term<C, E>(C(1), E(0x0003000001010100ULL)),
                    Term<C, E>(C(1), E(0x0003010000010100ULL)),
                    Term<C, E>(C(1), E(0x0003010100000100ULL)),
                });
                sys.push_back(p3);
                Polynomial<C, E> p4({
                    Term<C, E>(C(1), E(0x0004010101010000ULL)),
                    Term<C, E>(C(1), E(0x0004000101010100ULL)),
                    Term<C, E>(C(1), E(0x0004010001010100ULL)),
                    Term<C, E>(C(1), E(0x0004010100010100ULL)),
                    Term<C, E>(C(1), E(0x0004010101000100ULL)),
                });
                sys.push_back(p4);
                Polynomial<C, E> p5({
                    Term<C, E>(C(1), E(0x0005010101010100ULL)),
                    Term<C, E>(C(C::modulus - 1), E(0x0ULL))
                });
                sys.push_back(p5);

        } else if (n == 6) {
                // Cyclic 6 BitPacked
                Polynomial<C, E> p1({
                    Term<C, E>(C(1), E(0x0001010000000000ULL)),
                    Term<C, E>(C(1), E(0x0001000100000000ULL)),
                    Term<C, E>(C(1), E(0x0001000001000000ULL)),
                    Term<C, E>(C(1), E(0x0001000000010000ULL)),
                    Term<C, E>(C(1), E(0x0001000000000100ULL)),
                    Term<C, E>(C(1), E(0x0001000000000001ULL)),
                });
                sys.push_back(p1);
                Polynomial<C, E> p2({
                    Term<C, E>(C(1), E(0x0002010100000000ULL)),
                    Term<C, E>(C(1), E(0x0002000101000000ULL)),
                    Term<C, E>(C(1), E(0x0002000001010000ULL)),
                    Term<C, E>(C(1), E(0x0002000000010100ULL)),
                    Term<C, E>(C(1), E(0x0002000000000101ULL)),
                    Term<C, E>(C(1), E(0x0002010000000001ULL)),
                });
                sys.push_back(p2);
                Polynomial<C, E> p3({
                    Term<C, E>(C(1), E(0x0003010101000000ULL)),
                    Term<C, E>(C(1), E(0x0003000101010000ULL)),
                    Term<C, E>(C(1), E(0x0003000001010100ULL)),
                    Term<C, E>(C(1), E(0x0003000000010101ULL)),
                    Term<C, E>(C(1), E(0x0003010000000101ULL)),
                    Term<C, E>(C(1), E(0x0003010100000001ULL)),
                });
                sys.push_back(p3);
                Polynomial<C, E> p4({
                    Term<C, E>(C(1), E(0x0004010101010000ULL)),
                    Term<C, E>(C(1), E(0x0004000101010100ULL)),
                    Term<C, E>(C(1), E(0x0004000001010101ULL)),
                    Term<C, E>(C(1), E(0x0004010000010101ULL)),
                    Term<C, E>(C(1), E(0x0004010100000101ULL)),
                    Term<C, E>(C(1), E(0x0004010101000001ULL)),
                });
                sys.push_back(p4);
                Polynomial<C, E> p5({
                    Term<C, E>(C(1), E(0x0005010101010100ULL)),
                    Term<C, E>(C(1), E(0x0005000101010101ULL)),
                    Term<C, E>(C(1), E(0x0005010001010101ULL)),
                    Term<C, E>(C(1), E(0x0005010100010101ULL)),
                    Term<C, E>(C(1), E(0x0005010101000101ULL)),
                    Term<C, E>(C(1), E(0x0005010101010001ULL)),
                });
                sys.push_back(p5);
                Polynomial<C, E> p6({
                    Term<C, E>(C(1), E(0x0006010101010101ULL)),
                    Term<C, E>(C(C::modulus - 1), E(0x0ULL))
                });
                sys.push_back(p6);

        } else {
            std::cout << "Invalid n for hardcoded benchmark" << std::endl; return;
        }
    }

    std::vector<Polynomial<C, E>> gb;
    for (int i = 0; i < 10; ++i) {
        gb = naive_grobner_basis(sys);
        std::cout << "Iteration " << i << " complete" << std::endl;
    }
    //std::cout << "Basis size: " << gb.size() << std::endl;
}

int main(int argc, char* argv[]) {
    int n = 4;
    int mode = 0;
    
    if (argc > 1) {
        n = std::atoi(argv[1]);
    }
    if (argc > 2) {
        mode = std::atoi(argv[2]);
    }
    
    IntModP::modulus = 7;
    TERM_ORDER = Lex;

    if (mode == 0) {
        std::cout << "C++ generic vec exponent cyclic " << n << std::endl;
        run_benchmark<IntModP, VecExponent>(n);
    } else {
        std::cout << "C++ generic bit packed exponent cyclic " << n << std::endl;
        run_benchmark<IntModP, BitPackedExponent>(n);
    }
    return 0;
}

