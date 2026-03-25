include("generic/IField.jl")
include("generic/IExponent.jl")
include("generic/int_mod_p.jl")
include("generic/double_field.jl")
include("generic/single_field.jl")
include("generic/vec_exponent.jl")
include("generic/bit_packed_exponent.jl")

module GenGrobner
    import ..IField
    import ..IExponent
    import ..IntModP
    import ..DoubleField
    import ..SingleField
    import ..VecExponent
    # import ..BitPackedExponent
    import ..is_zero, ..can_reduce, ..get_lcm, ..lex_compare, ..degree

    export Term, Polynomial, TermOrder, set_term_order, grobner_basis, s_polynomial, reduce_polynomial

    mutable struct Term{C <: IField, E <: IExponent}
        coefficient::C
        exponents::E
    end

    @enum TermOrder Lex GrLex RevLex

    const TERM_ORDER = Ref(Lex)

    function set_term_order(order::TermOrder)
        TERM_ORDER[] = order
    end

    function compare(t1::Term{C, E}, t2::Term{C, E}) where {C <: IField, E <: IExponent}
        order = TERM_ORDER[]
        if order == Lex
            return lex_compare(t1.exponents, t2.exponents)
        elseif order == GrLex
            d1 = degree(t1.exponents)
            d2 = degree(t2.exponents)
            if d1 != d2
                return cmp(d1, d2)
            else
                return lex_compare(t1.exponents, t2.exponents)
            end
        elseif order == RevLex
            d1 = degree(t1.exponents)
            d2 = degree(t2.exponents)
            if d1 != d2
                return cmp(d1, d2)
            else
                return lex_compare(t2.exponents, t1.exponents) # RevLex
            end
        end
    end

    mutable struct Polynomial{C <: IField, E <: IExponent}
        terms::Vector{Term{C, E}}
    end

    function Polynomial(terms::Vector{Term{C, E}}) where {C, E}
        sort!(terms, by=t -> t, lt= (a,b) -> compare(a,b) == 1)
        filter!(t -> !is_zero(t.coefficient), terms)
        Polynomial{C, E}(terms)
    end

    import Base: +, -, *

    function +(p1::Polynomial{C, E}, p2::Polynomial{C, E}) where {C, E}
        result = deepcopy(p1.terms)
        for term in p2.terms
            found = false
            for res_term in result
                if lex_compare(res_term.exponents, term.exponents) == 0
                    res_term.coefficient = res_term.coefficient + term.coefficient
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

    function -(p1::Polynomial{C, E}, p2::Polynomial{C, E}) where {C, E}
        neg_p2_terms = [Term(Base.zero(term.coefficient) - term.coefficient, term.exponents) for term in p2.terms]
        p1 + Polynomial(neg_p2_terms)
    end

    function *(p1::Polynomial{C, E}, p2::Polynomial{C, E}) where {C, E}
        new_terms = Term{C, E}[]
        for t1 in p1.terms
            for t2 in p2.terms
                coeff = t1.coefficient * t2.coefficient
                exps = t1.exponents + t2.exponents
                push!(new_terms, Term(coeff, exps))
            end
        end
        # combine like terms
        combined_terms = Term{C, E}[]
        for term in new_terms
            found = false
            for c_term in combined_terms
                if lex_compare(c_term.exponents, term.exponents) == 0
                    c_term.coefficient = c_term.coefficient + term.coefficient
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
        return p.terms[1]
    end

    function s_polynomial(p1::Polynomial{C, E}, p2::Polynomial{C, E}) where {C, E}
        lt1 = leading_term(p1)
        lt2 = leading_term(p2)

        lcm_exponents = get_lcm(lt1.exponents, lt2.exponents)

        s1_exp = lcm_exponents - lt1.exponents
        s2_exp = lcm_exponents - lt2.exponents

        s1_term = Term(lt2.coefficient, s1_exp)
        s2_term = Term(lt1.coefficient, s2_exp)

        term_poly1 = Polynomial([s1_term])
        term_poly2 = Polynomial([s2_term])

        return (p1 * term_poly1) - (p2 * term_poly2)
    end

    function make_monic(p::Polynomial{C, E}) where {C, E}
        if isempty(p.terms)
            return p
        end
        lt_coeff = leading_term(p).coefficient
        inv_coeff = Base.one(C) / lt_coeff # C needs inverse or division defined
        new_terms = [Term(t.coefficient * inv_coeff, t.exponents) for t in p.terms]
        return Polynomial(new_terms)
    end

    function reduce_polynomial(p::Polynomial{C, E}, basis::Vector{Polynomial{C, E}}) where {C, E}
        poly = deepcopy(p)
        remainder = Term{C, E}[]

        while !isempty(poly.terms)
            reduced = false
            lt_poly = leading_term(poly)

            for b in basis
                if isempty(b.terms)
                    continue
                end
                lt_b = leading_term(b)
                if can_reduce(lt_poly.exponents, lt_b.exponents)
                    scale_coeff = lt_poly.coefficient / lt_b.coefficient
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

    function grobner_basis(polynomials::Vector{Polynomial{C, E}}) where {C, E}
        basis = deepcopy(polynomials)
        pair_set = [(i, j) for i in 1:length(basis) for j in (i+1):length(basis)]

        while !isempty(pair_set)
            (i, j) = popfirst!(pair_set)
            s_poly = s_polynomial(basis[i], basis[j])
            if !isempty(s_poly.terms)
                h = reduce_polynomial(s_poly, basis)
            else
                h = Polynomial(Term{C,E}[])
            end

            if !isempty(h.terms)
                push!(basis, h)
                n = length(basis)
                for k in 1:(n-1)
                    push!(pair_set, (k, n))
                end
            end
        end

        reduced_basis = Polynomial{C, E}[]
        for i in 1:length(basis)
            basis_excluding_self = vcat(basis[1:i-1], basis[i+1:end])
            reduced = reduce_polynomial(basis[i], basis_excluding_self)
            if !isempty(reduced.terms)
                push!(reduced_basis, make_monic(reduced))
            end
        end

        return reduced_basis
    end

    function main()
        println("Generic Grobner basis execution...")
        Main.set_modulus(UInt64(7))
        p1 = Polynomial([Term(IntModP(1), VecExponent([2, 0, 0])), Term(IntModP(6), VecExponent([0, 1, 0]))])
        p2 = Polynomial([Term(IntModP(1), VecExponent([3, 0, 0])), Term(IntModP(6), VecExponent([0, 0, 1]))])
        basis = grobner_basis([p1, p2])
        println("Basis size: ", length(basis))
        for b in basis
            for t in b.terms
                print("(", t.coefficient.value, " * ", t.exponents.exps, ") ")
            end
            println()
        end
        IntModP.Main.set_modulus(UInt64(7))
        p1 = Polynomial([Term(IntModP(1), VecExponent([2, 0, 0])), Term(IntModP(6), VecExponent([0, 1, 0]))])
        p2 = Polynomial([Term(IntModP(1), VecExponent([3, 0, 0])), Term(IntModP(6), VecExponent([0, 0, 1]))])
        basis = grobner_basis([p1, p2])
        println("Basis size: ", length(basis))
        for b in basis
            for t in b.terms
                print("(", t.coefficient.value, " * ", t.exponents.exps, ") ")
            end
            println()
        end
        # Add basic test to make sure it works
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    GenGrobner.main()
end
