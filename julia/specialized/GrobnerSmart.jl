module GrobnerSmart
    import Base: ==, hash
    export Term, Polynomial, TermOrder, set_term_order, set_modulus, grobner_basis, s_polynomial, reduce_polynomial

    mutable struct Term
        coefficient::UInt64
        exponents::UInt64
    end

    ==(t1::Term, t2::Term) = t1.coefficient == t2.coefficient && t1.exponents == t2.exponents
    hash(t::Term, h::UInt) = hash(t.coefficient, hash(t.exponents, h))

    @enum TermOrder Lex GrLex RevLex

    const TERM_ORDER = Ref(Lex)
    const MODULUS = Ref{UInt64}(7)

    function set_term_order(order::TermOrder)
        TERM_ORDER[] = order
    end

    function set_modulus(m::UInt64)
        MODULUS[] = m
    end

    function get_modulus()
        return MODULUS[]
    end

    function mod_inverse(a::UInt64, m::UInt64)
        a_i64 = Int64(a)
        m_i64 = Int64(m)
        m0 = m_i64
        y = 0
        x = 1

        if m == 1
            return UInt64(0)
        end

        while a_i64 > 1
            q = div(a_i64, m_i64)
            t = m_i64

            m_i64 = a_i64 % m_i64
            a_i64 = t
            t = y

            y = x - q * y
            x = t
        end

        if x < 0
            x = x + m0
        end

        return UInt64(x % m0)
    end

    function compare(t1::Term, t2::Term)
        order = TERM_ORDER[]
        if order == Lex
            e1 = t1.exponents & 0x0000FFFFFFFFFFFF
            e2 = t2.exponents & 0x0000FFFFFFFFFFFF
            return cmp(e1, e2)
        elseif order == GrLex
            d1 = (t1.exponents >> 48) & 0xFFFF
            d2 = (t2.exponents >> 48) & 0xFFFF
            if d1 != d2
                return cmp(d1, d2)
            else
                e1 = t1.exponents & 0x0000FFFFFFFFFFFF
                e2 = t2.exponents & 0x0000FFFFFFFFFFFF
                return cmp(e1, e2)
            end
        elseif order == RevLex
            # We don't strictly use RevLex in benchmark natively, but falling back for now
            d1 = (t1.exponents >> 48) & 0xFFFF
            d2 = (t2.exponents >> 48) & 0xFFFF
            if d1 != d2
                return cmp(d1, d2)
            else
                e1 = t1.exponents & 0x0000FFFFFFFFFFFF
                e2 = t2.exponents & 0x0000FFFFFFFFFFFF
                return cmp(e2, e1) # crude revlex approximation
            end
        end
    end

    mutable struct Polynomial
        terms::Vector{Term}
        function Polynomial(terms::Vector{Term})
            sort!(terms, by=t -> t, lt= (a,b) -> compare(a,b) == 1)
            filter!(t -> t.coefficient != 0, terms)
            new(terms)
        end
    end

    ==(p1::Polynomial, p2::Polynomial) = p1.terms == p2.terms
    hash(p::Polynomial, h::UInt) = hash(p.terms, h)

    import Base: +, -, *

    function +(p1::Polynomial, p2::Polynomial)
        result = deepcopy(p1.terms)
        m = MODULUS[]
        for term in p2.terms
            found = false
            for res_term in result
                if res_term.exponents == term.exponents
                    res_term.coefficient = (res_term.coefficient + term.coefficient) % m
                    found = true
                    break
                end
            end
            if !found
                push!(result, deepcopy(term))
            end
        end
        Polynomial(result)
    end

    function -(p1::Polynomial, p2::Polynomial)
        m = MODULUS[]
        neg_p2_terms = [Term((m - (term.coefficient % m)) % m, term.exponents) for term in p2.terms]
        p1 + Polynomial(neg_p2_terms)
    end

    function *(p1::Polynomial, p2::Polynomial)
        new_terms = Term[]
        m = MODULUS[]
        for t1 in p1.terms
            for t2 in p2.terms
                coeff = (t1.coefficient * t2.coefficient) % m
                exps = t1.exponents + t2.exponents
                push!(new_terms, Term(coeff, exps))
            end
        end
        # combine like terms
        combined_terms = Term[]
        for term in new_terms
            found = false
            for c_term in combined_terms
                if c_term.exponents == term.exponents
                    c_term.coefficient = (c_term.coefficient + term.coefficient) % m
                    found = true
                    break
                end
            end
            if !found
                push!(combined_terms, term)
            end
        end
        Polynomial(combined_terms)
    end

    function leading_term(p::Polynomial)
        if isempty(p.terms)
            return Term(UInt64(0), UInt64(0))
        end
        return p.terms[1]
    end

    function get_lcm(a::UInt64, b::UInt64)
        self_exponents = a & 0x0000FFFFFFFFFFFF
        other_exponents = b & 0x0000FFFFFFFFFFFF
        
        lcm_exponents = UInt64(0)
        deg = UInt64(0)
        
        for i in 0:8:40
            self_exp = (self_exponents >> i) & 0xFF
            other_exp = (other_exponents >> i) & 0xFF
            lcm_exp = max(self_exp, other_exp)
            lcm_exponents |= (lcm_exp << i)
            deg += lcm_exp
        end
        
        lcm_exponents |= (deg & 0xFFFF) << 48
        return lcm_exponents
    end

    function can_reduce(a::UInt64, b::UInt64)
        for i in 0:8:40
            self_exp = (a >> i) & 0xFF
            divisor_exp = (b >> i) & 0xFF
            if self_exp < divisor_exp
                return false
            end
        end
        return true
    end

    function s_polynomial(p1::Polynomial, p2::Polynomial)
        lt1 = leading_term(p1)
        lt2 = leading_term(p2)

        if isempty(p1.terms) || isempty(p2.terms)
            return Polynomial(Term[])
        end

        lcm_exponents = get_lcm(lt1.exponents, lt2.exponents)

        s1_exp = lcm_exponents - lt1.exponents
        s2_exp = lcm_exponents - lt2.exponents

        m = MODULUS[]
        inv1 = mod_inverse(lt1.coefficient, m)
        inv2 = mod_inverse(lt2.coefficient, m)

        s1_term = Term(lt2.coefficient, s1_exp)
        s2_term = Term(lt1.coefficient, s2_exp)

        term_poly1 = Polynomial([s1_term])
        term_poly2 = Polynomial([s2_term])

        return (p1 * term_poly1) - (p2 * term_poly2)
    end

    function make_monic(p::Polynomial)
        if isempty(p.terms)
            return p
        end
        m = MODULUS[]
        lt_coeff = leading_term(p).coefficient
        inv_coeff = mod_inverse(lt_coeff, m)
        new_terms = [Term((t.coefficient * inv_coeff) % m, t.exponents) for t in p.terms]
        return Polynomial(new_terms)
    end

    function reduce_polynomial(p::Polynomial, basis::Vector{Polynomial})
        poly = deepcopy(p)
        m = MODULUS[]
        remainder = Term[]

        while !isempty(poly.terms)
            reduced = false
            lt_poly = leading_term(poly)

            for b in basis
                if isempty(b.terms)
                    continue
                end
                lt_b = leading_term(b)
                if can_reduce(lt_poly.exponents, lt_b.exponents)
                    inv_b = mod_inverse(lt_b.coefficient, m)
                    scale_coeff = (lt_poly.coefficient * inv_b) % m
                    scale_exps = lt_poly.exponents - lt_b.exponents

                    scale_term = Term(scale_coeff, scale_exps)
                    scaled_b = b * Polynomial([scale_term])
                    poly = poly - scaled_b
                    reduced = true
                    break
                end
            end

            if !reduced
                push!(remainder, deepcopy(leading_term(poly)))
                popfirst!(poly.terms)
            end
        end
        return Polynomial(remainder)
    end

    function grobner_basis(polynomials::Vector{Polynomial})
        basis = deepcopy(polynomials)
        basis_set = Set{Polynomial}(basis)
        n_orig = length(basis)
        pair_set = [(i, j) for i in 1:length(basis) for j in (i+1):length(basis)]

        while !isempty(pair_set)
            (i, j) = popfirst!(pair_set)
            s_poly = s_polynomial(basis[i], basis[j])
            h = reduce_polynomial(s_poly, basis)

            if !isempty(h.terms) && !in(h, basis_set)
                push!(basis, h)
                push!(basis_set, h)
                n = length(basis)
                for k in 1:(n-1)
                    push!(pair_set, (k, n))
                end
            end
        end

        reduced_basis = Polynomial[]
        for i in 1:length(basis)
            basis_excluding_self = vcat(basis[1:i-1], basis[i+1:end])
            reduced = reduce_polynomial(basis[i], basis_excluding_self)
            if !isempty(reduced.terms)
                push!(reduced_basis, make_monic(reduced))
            end
        end

        return reduced_basis
    end

    function pack_exponents(exponents::Vector{Int})
        packed = UInt64(0)
        padded_exps = copy(exponents)
        while length(padded_exps) < 6
            push!(padded_exps, 0)
        end
        for i in 1:6
            exp = padded_exps[i]
            shift = UInt64(40 - 8 * (i - 1))
            packed |= (UInt64(exp) & 0xFF) << shift
        end
        deg = UInt64(sum(padded_exps[1:6]))
        packed |= (deg & 0xFFFF) << 48
        return packed
    end

    function main(n::Int)
        set_modulus(UInt64(7))
        polynomials = Polynomial[]
        if n == 4
            println("Julia specialized bit-packed cyclic 4")
            p1 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 0, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 0, 1])),
            ])
            p2 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 0, 0, 1])),
            ])
            p3 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 0, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 1, 0, 1])),
            ])
            p4 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 1, 1])),
                Term(UInt64(get_modulus() - 1), pack_exponents([0, 0, 0, 0])),
            ])
            polynomials = [p1, p2, p3, p4]
        end
        if n == 5
            println("Julia specialized bit-packed cyclic 5")
            p1 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 0, 0, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 0, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 1, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 0, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 0, 0, 1])),
            ])
            p2 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 0, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 1, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 1, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 0, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 0, 0, 0, 1])),
            ])
            p3 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 1, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 1, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 1, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 0, 0, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 1, 0, 0, 1])),
            ])
            p4 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 1, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 1, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 0, 1, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 1, 0, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 1, 1, 0, 1])),
            ])
            p5 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 1, 1, 1])),
                Term(UInt64(get_modulus() - 1), pack_exponents([0, 0, 0, 0, 0])),
            ])

            polynomials = [p1, p2, p3, p4, p5]
        end
        if n == 6
            println("Julia specialized bit-packed cyclic 6")
            p1 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 0, 0, 0, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 0, 0, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 1, 0, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 0, 1, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 0, 0, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 0, 0, 0, 1])),
            ])
            p2 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 0, 0, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 1, 0, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 1, 1, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 0, 1, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 0, 0, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 0, 0, 0, 0, 1])),
            ])
            p3 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 1, 0, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 1, 1, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 1, 1, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 0, 1, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 0, 0, 0, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 1, 0, 0, 0, 1])),
            ])
            p4 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 1, 1, 0, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 1, 1, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 0, 1, 1, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 0, 0, 1, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 1, 0, 0, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 1, 1, 0, 0, 1])),
            ])
            p5 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 1, 1, 1, 0])),
                Term(UInt64(1), pack_exponents([0, 1, 1, 1, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 0, 1, 1, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 1, 0, 1, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 1, 1, 0, 1, 1])),
                Term(UInt64(1), pack_exponents([1, 1, 1, 1, 0, 1])),
            ])
            p6 = Polynomial([
                Term(UInt64(1), pack_exponents([1, 1, 1, 1, 1, 1])),
                Term(UInt64(get_modulus() - 1), pack_exponents([0, 0, 0, 0, 0, 0])),
            ])

            polynomials = [p1, p2, p3, p4, p5, p6]
        end

        set_term_order(Lex)
        for i in 0:9 
            basis = grobner_basis(polynomials)
            println("Iteration $i complete")
            #if i == 9
                #println("Grobner basis computed with size: ", length(basis))
            #end
        end
        
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    n = 4 # Default to cyclic 4
    if length(ARGS) > 0
        n = parse(Int, ARGS[1])
    end
    GrobnerSmart.main(n)
end
