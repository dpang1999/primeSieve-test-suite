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
    vector<int> exponents;
    Term() : coefficient(0) {}
    Term(uint64_t c, vector<int> e) : coefficient(c), exponents(std::move(e)) {}
};

inline int compare_terms(const Term& t1, const Term& t2) {
    // Only Lex implemented as required
    for (size_t i = 0; i < t1.exponents.size() && i < t2.exponents.size(); ++i) {
        if (t1.exponents[i] < t2.exponents[i]) return -1;
        if (t1.exponents[i] > t2.exponents[i]) return 1;
    }
    if (t1.exponents.size() < t2.exponents.size()) return -1;
    if (t1.exponents.size() > t2.exponents.size()) return 1;
    return 0;
}

class Polynomial {
public:
    vector<Term> terms;
    Polynomial() {}
    Polynomial(vector<Term> _terms) {
        terms.reserve(_terms.size());
        for (auto& t : _terms) {
            if (t.coefficient != 0) {
                terms.push_back(std::move(t));
            }
        }
        sort(terms.begin(), terms.end(), [](const Term& a, const Term& b) {
            return compare_terms(a, b) > 0;
        });
    }
};

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
            for (int e : t.exponents) {
                hash ^= e + 0x9e3779b9 + (hash << 6) + (hash >> 2);
            }
        }
        return hash;
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
    return Polynomial(std::move(result));
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
            result.push_back(Term((MODULUS - (term.coefficient % MODULUS)) % MODULUS, term.exponents));
        }
    }
    return Polynomial(std::move(result));
}


Polynomial multiply_by_term(const Polynomial& p, const Term& term) {
    vector<Term> new_terms;
    new_terms.reserve(p.terms.size());
    for (const auto& t : p.terms) {
        vector<int> exps(t.exponents.size());
        for (size_t i = 0; i < exps.size(); ++i) {
            exps[i] = t.exponents[i] + term.exponents[i];
        }
        new_terms.push_back(Term((t.coefficient * term.coefficient) % MODULUS, std::move(exps)));
    }
    return Polynomial(std::move(new_terms));
}


Term leading_term(const Polynomial& p) {
    if (p.terms.empty()) return {0, vector<int>()};
    return p.terms[0];
}

Polynomial s_polynomial(const Polynomial& p1, const Polynomial& p2) {
    Term lt1 = leading_term(p1);
    Term lt2 = leading_term(p2);
    if (p1.terms.empty() || p2.terms.empty()) return Polynomial();

    vector<int> lcm_exps(lt1.exponents.size());
    vector<int> s1_exp(lt1.exponents.size());
    vector<int> s2_exp(lt2.exponents.size());
    for (size_t i = 0; i < lcm_exps.size(); ++i) {
        lcm_exps[i] = max(lt1.exponents[i], lt2.exponents[i]);
        s1_exp[i] = lcm_exps[i] - lt1.exponents[i];
        s2_exp[i] = lcm_exps[i] - lt2.exponents[i];
    }

    Term s1_term(lt2.coefficient, std::move(s1_exp));
    Term s2_term(lt1.coefficient, std::move(s2_exp));

    return sub_poly(multiply_by_term(p1, s1_term), multiply_by_term(p2, s2_term));
}

Polynomial make_monic(const Polynomial& p) {
    if (p.terms.empty()) return p;
    uint64_t lt_coeff = leading_term(p).coefficient;
    uint64_t inv_coeff = mod_inverse(lt_coeff, MODULUS);
    vector<Term> new_terms;
    for (const auto& t : p.terms) {
        new_terms.push_back(Term((t.coefficient * inv_coeff) % MODULUS, t.exponents));
    }
    return Polynomial(std::move(new_terms));
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
            bool can_reduce = true;
            for (size_t i = 0; i < lt_poly.exponents.size(); ++i) {
                if (lt_poly.exponents[i] < lt_b.exponents[i]) {
                    can_reduce = false;
                    break;
                }
            }
            if (can_reduce) {
                uint64_t inv_b = mod_inverse(lt_b.coefficient, MODULUS);
                uint64_t scale_coeff = (lt_poly.coefficient * inv_b) % MODULUS;
                vector<int> scale_exps(lt_poly.exponents.size());
                for (size_t i = 0; i < scale_exps.size(); ++i) {
                    scale_exps[i] = lt_poly.exponents[i] - lt_b.exponents[i];
                }
                Polynomial scaled_b = multiply_by_term(b, {scale_coeff, scale_exps});
                poly = sub_poly(poly, scaled_b);
                reduced = true;
                break;
            }
        }
        if (!reduced) {
            remainder.push_back(std::move(poly.terms[0]));
            poly.terms.erase(poly.terms.begin());
        }
    }
    return Polynomial(std::move(remainder));
}

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
        // h = make_monic(h);
        
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

int main(int argc, char** argv) {
    int n = 4;
    if (argc > 1) {
        n = stoi(argv[1]);
    }
    
    MODULUS = 7;
    TERM_ORDER = Lex;
    vector<Polynomial> polynomials;

    if (n == 4) {
        cout << "C++ specialized vec exponent cyclic 4\n";
        polynomials = {
            Polynomial({ {1, {1, 0, 0, 0}}, {1, {0, 1, 0, 0}}, {1, {0, 0, 1, 0}}, {1, {0, 0, 0, 1}} }),
            Polynomial({ {1, {1, 1, 0, 0}}, {1, {0, 1, 1, 0}}, {1, {0, 0, 1, 1}}, {1, {1, 0, 0, 1}} }),
            Polynomial({ {1, {1, 1, 1, 0}}, {1, {0, 1, 1, 1}}, {1, {1, 0, 1, 1}}, {1, {1, 1, 0, 1}} }),
            Polynomial({ {1, {1, 1, 1, 1}}, {6, {0, 0, 0, 0}} })
        };
    } else if (n == 5) {
        cout << "C++ specialized vec exponent cyclic 5\n";
        polynomials = {
            Polynomial({ {1, {1, 0, 0, 0, 0}}, {1, {0, 1, 0, 0, 0}}, {1, {0, 0, 1, 0, 0}}, {1, {0, 0, 0, 1, 0}}, {1, {0, 0, 0, 0, 1}} }),
            Polynomial({ {1, {1, 1, 0, 0, 0}}, {1, {0, 1, 1, 0, 0}}, {1, {0, 0, 1, 1, 0}}, {1, {0, 0, 0, 1, 1}}, {1, {1, 0, 0, 0, 1}} }),
            Polynomial({ {1, {1, 1, 1, 0, 0}}, {1, {0, 1, 1, 1, 0}}, {1, {0, 0, 1, 1, 1}}, {1, {1, 0, 0, 1, 1}}, {1, {1, 1, 0, 0, 1}} }),
            Polynomial({ {1, {1, 1, 1, 1, 0}}, {1, {0, 1, 1, 1, 1}}, {1, {1, 0, 1, 1, 1}}, {1, {1, 1, 0, 1, 1}}, {1, {1, 1, 1, 0, 1}} }),
            Polynomial({ {1, {1, 1, 1, 1, 1}}, {6, {0, 0, 0, 0, 0}} })
        };
    } else if (n == 6) {
        cout << "C++ specialized vec exponent cyclic 6\n";
        polynomials = {
            Polynomial({ {1, {1, 0, 0, 0, 0, 0}}, {1, {0, 1, 0, 0, 0, 0}}, {1, {0, 0, 1, 0, 0, 0}}, {1, {0, 0, 0, 1, 0, 0}}, {1, {0, 0, 0, 0, 1, 0}}, {1, {0, 0, 0, 0, 0, 1}} }),
            Polynomial({ {1, {1, 1, 0, 0, 0, 0}}, {1, {0, 1, 1, 0, 0, 0}}, {1, {0, 0, 1, 1, 0, 0}}, {1, {0, 0, 0, 1, 1, 0}}, {1, {0, 0, 0, 0, 1, 1}}, {1, {1, 0, 0, 0, 0, 1}} }),
            Polynomial({ {1, {1, 1, 1, 0, 0, 0}}, {1, {0, 1, 1, 1, 0, 0}}, {1, {0, 0, 1, 1, 1, 0}}, {1, {0, 0, 0, 1, 1, 1}}, {1, {1, 0, 0, 0, 1, 1}}, {1, {1, 1, 0, 0, 0, 1}} }),
            Polynomial({ {1, {1, 1, 1, 1, 0, 0}}, {1, {0, 1, 1, 1, 1, 0}}, {1, {0, 0, 1, 1, 1, 1}}, {1, {1, 0, 0, 1, 1, 1}}, {1, {1, 1, 0, 0, 1, 1}}, {1, {1, 1, 1, 0, 0, 1}} }),
            Polynomial({ {1, {1, 1, 1, 1, 1, 0}}, {1, {0, 1, 1, 1, 1, 1}}, {1, {1, 0, 1, 1, 1, 1}}, {1, {1, 1, 0, 1, 1, 1}}, {1, {1, 1, 1, 0, 1, 1}}, {1, {1, 1, 1, 1, 0, 1}} }),
            Polynomial({ {1, {1, 1, 1, 1, 1, 1}}, {6, {0, 0, 0, 0, 0, 0}} })
        };
    }

    for (int i = 0; i < 10; ++i) {
        auto basis = grobner_basis(polynomials);
        cout << "Iteration " << i << " complete\n";
    }

    return 0;
}
