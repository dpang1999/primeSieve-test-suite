#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <chrono>
#include <unordered_set>

using namespace std;

enum TermOrder { Lex, GrLex, RevLex };
TermOrder TERM_ORDER = Lex;
uint64_t MODULUS = 7;

uint64_t mod_inverse(uint64_t a, uint64_t m) {
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

struct Term {
    uint64_t coefficient;
    uint64_t exponents;
};

int compare_terms(const Term& t1, const Term& t2) {
    uint64_t e1 = t1.exponents & 0x0000FFFFFFFFFFFFULL;
    uint64_t e2 = t2.exponents & 0x0000FFFFFFFFFFFFULL;
    if (e1 < e2) return -1;
    if (e1 > e2) return 1;
    return 0;
}

class Polynomial {
public:
    vector<Term> terms;
    Polynomial() {}
    Polynomial(vector<Term> _terms) {
        for (auto& t : _terms) {
            if (t.coefficient != 0) {
                terms.push_back(t);
            }
        }
        sort(terms.begin(), terms.end(), [](const Term& a, const Term& b) {
            return compare_terms(a, b) > 0;
        });
    }
};

Polynomial add_poly(const Polynomial& p1, const Polynomial& p2) {
    vector<Term> result = p1.terms;
    for (const auto& term2 : p2.terms) {
        bool found = false;
        for (auto& res_term : result) {
            if (res_term.exponents == term2.exponents) {
                res_term.coefficient = (res_term.coefficient + term2.coefficient) % MODULUS;
                found = true;
                break;
            }
        }
        if (!found) result.push_back(term2);
    }
    return Polynomial(result);
}

Polynomial sub_poly(const Polynomial& p1, const Polynomial& p2) {
    vector<Term> result = p1.terms;
    for (const auto& term : p2.terms) {
        bool found = false;
        for (auto& res_term : result) {
            if (res_term.exponents == term.exponents) {
                res_term.coefficient = (MODULUS + res_term.coefficient - term.coefficient) % MODULUS;
                found = true;
                break;
            }
        }
        if (!found) {
            result.push_back({(MODULUS - (term.coefficient % MODULUS)) % MODULUS, term.exponents});
        }
    }
    return Polynomial(result);
}


Polynomial multiply_by_term(const Polynomial& p, const Term& term) {
    vector<Term> new_terms;
    new_terms.reserve(p.terms.size());
    for (const auto& t : p.terms) {
        new_terms.push_back({(t.coefficient * term.coefficient) % MODULUS, t.exponents + term.exponents});
    }
    return Polynomial(new_terms);
}


Term leading_term(const Polynomial& p) {
    if (p.terms.empty()) return {0, 0};
    return p.terms[0];
}

uint64_t get_lcm(uint64_t a, uint64_t b) {
    uint64_t self_exponents = a & 0x0000FFFFFFFFFFFFULL;
    uint64_t other_exponents = b & 0x0000FFFFFFFFFFFFULL;
    uint64_t lcm_exponents = 0;
    uint64_t deg = 0;
    for (int i = 0; i <= 40; i += 8) {
        uint64_t self_exp = (self_exponents >> i) & 0xFF;
        uint64_t other_exp = (other_exponents >> i) & 0xFF;
        uint64_t lcm_exp = max(self_exp, other_exp);
        lcm_exponents |= (lcm_exp << i);
        deg += lcm_exp;
    }
    lcm_exponents |= (deg & 0xFFFF) << 48;
    return lcm_exponents;
}

bool can_reduce(uint64_t a, uint64_t b) {
    for (int i = 0; i <= 40; i += 8) {
        uint64_t self_exp = (a >> i) & 0xFF;
        uint64_t divisor_exp = (b >> i) & 0xFF;
        if (self_exp < divisor_exp) return false;
    }
    return true;
}

Polynomial s_polynomial(const Polynomial& p1, const Polynomial& p2) {
    Term lt1 = leading_term(p1);
    Term lt2 = leading_term(p2);
    if (p1.terms.empty() || p2.terms.empty()) return Polynomial();

    uint64_t lcm_exps = get_lcm(lt1.exponents, lt2.exponents);
    uint64_t s1_exp = lcm_exps - lt1.exponents;
    uint64_t s2_exp = lcm_exps - lt2.exponents;

    Term s1_term = {lt2.coefficient, s1_exp};
    Term s2_term = {lt1.coefficient, s2_exp};

    return sub_poly(multiply_by_term(p1, s1_term), multiply_by_term(p2, s2_term));
}

Polynomial make_monic(const Polynomial& p) {
    if (p.terms.empty()) return p;
    uint64_t lt_coeff = leading_term(p).coefficient;
    uint64_t inv_coeff = mod_inverse(lt_coeff, MODULUS);
    vector<Term> new_terms;
    for (const auto& t : p.terms) {
        new_terms.push_back({(t.coefficient * inv_coeff) % MODULUS, t.exponents});
    }
    return Polynomial(new_terms);
}

Polynomial reduce_polynomial(const Polynomial& p, const vector<Polynomial>& basis) {
    Polynomial poly = p;
    vector<Term> remainder;
    while (!poly.terms.empty()) {
        bool reduced = false;
        Term lt_poly = leading_term(poly);
        for (const auto& b : basis) {
            if (b.terms.empty()) continue;
            Term lt_b = leading_term(b);
            if (can_reduce(lt_poly.exponents, lt_b.exponents)) {
                uint64_t inv_b = mod_inverse(lt_b.coefficient, MODULUS);
                uint64_t scale_coeff = (lt_poly.coefficient * inv_b) % MODULUS;
                uint64_t scale_exps = lt_poly.exponents - lt_b.exponents;
                
                Polynomial scaled_b = multiply_by_term(b, {scale_coeff, scale_exps});
                poly = sub_poly(poly, scaled_b);
                reduced = true;
                break;
            }
        }
        if (!reduced) {
            remainder.push_back(leading_term(poly));
            poly.terms.erase(poly.terms.begin());
        }
    }
    return Polynomial(remainder);
}


bool polynomials_equal(const Polynomial& p1, const Polynomial& p2) {
    if (p1.terms.size() != p2.terms.size()) return false;
    for (size_t i = 0; i < p1.terms.size(); i++) {
        if (p1.terms[i].coefficient != p2.terms[i].coefficient || p1.terms[i].exponents != p2.terms[i].exponents) {
            return false;
        }
    }
    return true;
}

inline bool operator==(const Polynomial& p1, const Polynomial& p2) {
    return polynomials_equal(p1, p2);
}

struct PolyHash {
    size_t operator()(const Polynomial& p) const {
        size_t hash = 0;
        for (const auto& t : p.terms) {
            hash ^= t.coefficient + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            hash ^= t.exponents + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
        return hash;
    }
};

vector<Polynomial> grobner_basis(vector<Polynomial> basis) {
    unordered_set<Polynomial, PolyHash> basis_set;
    for (const auto& poly : basis) {
        basis_set.insert(poly);
    }

    vector<pair<int, int>> pair_set;
    for (size_t i = 0; i < basis.size(); i++) {
        for (size_t j = i + 1; j < basis.size(); j++) {
            pair_set.push_back({(int)i, (int)j});
        }
    }
    while (!pair_set.empty()) {
        auto p = pair_set.front();
        pair_set.erase(pair_set.begin());
        
        Polynomial s_poly = s_polynomial(basis[p.first], basis[p.second]);
        Polynomial h = reduce_polynomial(s_poly, basis);
        
        if (!h.terms.empty() && basis_set.insert(h).second) {
            int n = basis.size();
            basis.push_back(h);
            for (int k = 0; k < n; k++) {
                pair_set.push_back({k, n});
            }
        }
    }
    vector<Polynomial> reduced_basis;
    for (size_t i = 0; i < basis.size(); i++) {
        vector<Polynomial> basis_excluding_self;
        for (size_t j = 0; j < basis.size(); j++) {
            if (i != j) basis_excluding_self.push_back(basis[j]);
        }
        Polynomial reduced = reduce_polynomial(basis[i], basis_excluding_self);
        if (!reduced.terms.empty()) {
            reduced_basis.push_back(make_monic(reduced));
        }
    }
    return reduced_basis;
}

uint64_t pack_exponents(vector<int> exps) {
    uint64_t packed = 0;
    while(exps.size() < 6) exps.push_back(0);
    for (int i = 0; i < 6; i++) {
        uint64_t shift = 40 - 8 * i;
        packed |= ((uint64_t)exps[i] & 0xFF) << shift;
    }
    uint64_t deg = 0;
    for (int i = 0; i < 6; i++) deg += exps[i];
    packed |= (deg & 0xFFFF) << 48;
    return packed;
}

int main(int argc, char** argv) {
    int n = 4;
    if (argc > 1) {
        n = stoi(argv[1]);
    }
    
    MODULUS = 7;
    TERM_ORDER = Lex;
    vector<Polynomial> polynomials;

    if (n == 4) {
        cout << "C++ specialized bit-packed cyclic 4\n";
        polynomials = {
            Polynomial({ {1, pack_exponents({1, 0, 0, 0})}, {1, pack_exponents({0, 1, 0, 0})}, {1, pack_exponents({0, 0, 1, 0})}, {1, pack_exponents({0, 0, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 0, 0})}, {1, pack_exponents({0, 1, 1, 0})}, {1, pack_exponents({0, 0, 1, 1})}, {1, pack_exponents({1, 0, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 1, 0})}, {1, pack_exponents({0, 1, 1, 1})}, {1, pack_exponents({1, 0, 1, 1})}, {1, pack_exponents({1, 1, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 1, 1})}, {6, pack_exponents({0, 0, 0, 0})} })
        };
    } else if (n == 5) {
        cout << "C++ specialized bit-packed cyclic 5\n";
        polynomials = {
            Polynomial({ {1, pack_exponents({1, 0, 0, 0, 0})}, {1, pack_exponents({0, 1, 0, 0, 0})}, {1, pack_exponents({0, 0, 1, 0, 0})}, {1, pack_exponents({0, 0, 0, 1, 0})}, {1, pack_exponents({0, 0, 0, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 0, 0, 0})}, {1, pack_exponents({0, 1, 1, 0, 0})}, {1, pack_exponents({0, 0, 1, 1, 0})}, {1, pack_exponents({0, 0, 0, 1, 1})}, {1, pack_exponents({1, 0, 0, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 1, 0, 0})}, {1, pack_exponents({0, 1, 1, 1, 0})}, {1, pack_exponents({0, 0, 1, 1, 1})}, {1, pack_exponents({1, 0, 0, 1, 1})}, {1, pack_exponents({1, 1, 0, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 1, 1, 0})}, {1, pack_exponents({0, 1, 1, 1, 1})}, {1, pack_exponents({1, 0, 1, 1, 1})}, {1, pack_exponents({1, 1, 0, 1, 1})}, {1, pack_exponents({1, 1, 1, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 1, 1, 1})}, {6, pack_exponents({0, 0, 0, 0, 0})} })
        };
    } else if (n == 6) {
        cout << "C++ specialized bit-packed cyclic 6\n";
        polynomials = {
            Polynomial({ {1, pack_exponents({1, 0, 0, 0, 0, 0})}, {1, pack_exponents({0, 1, 0, 0, 0, 0})}, {1, pack_exponents({0, 0, 1, 0, 0, 0})}, {1, pack_exponents({0, 0, 0, 1, 0, 0})}, {1, pack_exponents({0, 0, 0, 0, 1, 0})}, {1, pack_exponents({0, 0, 0, 0, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 0, 0, 0, 0})}, {1, pack_exponents({0, 1, 1, 0, 0, 0})}, {1, pack_exponents({0, 0, 1, 1, 0, 0})}, {1, pack_exponents({0, 0, 0, 1, 1, 0})}, {1, pack_exponents({0, 0, 0, 0, 1, 1})}, {1, pack_exponents({1, 0, 0, 0, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 1, 0, 0, 0})}, {1, pack_exponents({0, 1, 1, 1, 0, 0})}, {1, pack_exponents({0, 0, 1, 1, 1, 0})}, {1, pack_exponents({0, 0, 0, 1, 1, 1})}, {1, pack_exponents({1, 0, 0, 0, 1, 1})}, {1, pack_exponents({1, 1, 0, 0, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 1, 1, 0, 0})}, {1, pack_exponents({0, 1, 1, 1, 1, 0})}, {1, pack_exponents({0, 0, 1, 1, 1, 1})}, {1, pack_exponents({1, 0, 0, 1, 1, 1})}, {1, pack_exponents({1, 1, 0, 0, 1, 1})}, {1, pack_exponents({1, 1, 1, 0, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 1, 1, 1, 0})}, {1, pack_exponents({0, 1, 1, 1, 1, 1})}, {1, pack_exponents({1, 0, 1, 1, 1, 1})}, {1, pack_exponents({1, 1, 0, 1, 1, 1})}, {1, pack_exponents({1, 1, 1, 0, 1, 1})}, {1, pack_exponents({1, 1, 1, 1, 0, 1})} }),
            Polynomial({ {1, pack_exponents({1, 1, 1, 1, 1, 1})}, {6, pack_exponents({0, 0, 0, 0, 0, 0})} })
        };
    }

    for (int i = 0; i < 10; ++i) {
        auto basis = grobner_basis(polynomials);
        cout << "Iteration " << i << " complete\n";
    }

    return 0;
}
