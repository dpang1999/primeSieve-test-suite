include("generic/IField.jl")
include("generic/IExponent.jl")
include("generic/int_mod_p.jl")
include("generic/vec_exponent.jl")
include("generic/bit_packed_exponent.jl")

module GenGrobner
    import ..IField
    import ..IExponent
    import ..IntModP
    import ..VecExponent
    import ..BitPackedExponent
    import ..is_zero, ..can_reduce, ..get_lcm, ..lex_compare, ..degree

    import ..set_modulus, ..get_modulus

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
        function Polynomial{C, E}(terms::Vector{Term{C, E}}) where {C <: IField, E <: IExponent}
            sort!(terms, by=t -> t, lt= (a,b) -> compare(a,b) == 1)
            filter!(t -> !is_zero(t.coefficient), terms)
            new{C, E}(terms)
        end
    end

    function Polynomial(terms::Vector{Term{C, E}}) where {C <: IField, E <: IExponent}
        Polynomial{C, E}(terms)
    end

    import Base: ==, hash

    function ==(t1::Term{C, E}, t2::Term{C, E}) where {C <: IField, E <: IExponent}
        t1.coefficient == t2.coefficient && t1.exponents == t2.exponents
    end

    function hash(t::Term{C, E}, h::UInt) where {C <: IField, E <: IExponent}
        hash(t.coefficient, hash(t.exponents, h))
    end

    function ==(p1::Polynomial{C, E}, p2::Polynomial{C, E}) where {C <: IField, E <: IExponent}
        length(p1.terms) == length(p2.terms) && all(t1 == t2 for (t1, t2) in zip(p1.terms, p2.terms))
    end

    function hash(p::Polynomial{C, E}, h::UInt) where {C <: IField, E <: IExponent}
        h2 = hash(length(p.terms), h)
        for t in p.terms
            h2 = hash(t, h2)
        end
        return h2
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
        basis_set = Set{Polynomial{C, E}}(basis)
        pair_set = [(i, j) for i in 1:length(basis) for j in (i+1):length(basis)]

        while !isempty(pair_set)
            (i, j) = popfirst!(pair_set)
            s_poly = s_polynomial(basis[i], basis[j])
            if !isempty(s_poly.terms)
                h = reduce_polynomial(s_poly, basis)
            else
                h = Polynomial(Term{C,E}[])
            end

            if !isempty(h.terms) && !in(h, basis_set)
                push!(basis, h)
                push!(basis_set, h)
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

    function main(n::Int, mode::Int)
        # mode = 0 for vec, 1 for bit packed
        set_modulus(UInt64(7))
        basis = []
        if n == 4
            if mode == 0
                # cyclic 4
                println("Julia generic vec exponent cyclic 4")
                p1 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 0, 0, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 0, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 1, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 0, 1])),
                ])
                p2 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 0, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 1, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 0, 0, 1])),
                ])
                p3 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 1, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 0, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 1, 0, 1])),
                ])
                p4 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 1, 1])), 
                    Term(IntModP(get_modulus() - 1), VecExponent([0, 0, 0, 0])),
                ])
            else
                # cyclic 4 with bit packed exponents
                println("Julia generic bit packed exponent cyclic 4")
                p1 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000000000001)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000000000100)),
                    Term(IntModP(1), BitPackedExponent(0x0000000000010000)),
                    Term(IntModP(1), BitPackedExponent(0x0000000001000000)),
                ])
                p2 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000000000101)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000000010100)),
                    Term(IntModP(1), BitPackedExponent(0x0000000001010000)),
                    Term(IntModP(1), BitPackedExponent(0x0000000001000001)),
                ])
                p3 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000000010101)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000001010100)),
                    Term(IntModP(1), BitPackedExponent(0x0000000001010001)),
                    Term(IntModP(1), BitPackedExponent(0x0000000001000101)),
                ])
                p4 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000001010101)), 
                    Term(IntModP(get_modulus() - 1), BitPackedExponent(0x0000000000000000)),
                ])
            end
            for i in 0:9
                basis = grobner_basis([p1, p2, p3, p4])
                println("Iteration $i complete")
            end
        elseif n == 5
            if mode == 0
                # cyclic 5
                println("Julia generic vec exponent cyclic 5")
                p1 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 0, 0, 0, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 0, 0, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 1, 0, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 0, 1, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 0, 0, 1])),
                ])
                p2 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 0, 0, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 1, 0, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 1, 1, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 0, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 0, 0, 0, 1])),
                ])
                p3 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 1, 0, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 1, 1, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 1, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 0, 0, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 1, 0, 0, 1])),
                ])
                p4 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 1, 1, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 1, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 0, 1, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 1, 0, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 1, 1, 0, 1])),
                ])
                p5 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 1, 1, 1])), 
                    Term(IntModP(get_modulus() - 1), VecExponent([0, 0, 0, 0, 0])),
                ])
            else
                # cyclic 5 with bit packed exponents
                println("Julia generic bit packed exponent cyclic 5")
                p1 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000000000001)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000000000100)),
                    Term(IntModP(1), BitPackedExponent(0x0000000000010000)),
                    Term(IntModP(1), BitPackedExponent(0x0000000001000000)),
                    Term(IntModP(1), BitPackedExponent(0x0000000100000000)),
                ])
                p2 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000000000101)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000000010100)),
                    Term(IntModP(1), BitPackedExponent(0x0000000001010000)),
                    Term(IntModP(1), BitPackedExponent(0x0000000101000000)),
                    Term(IntModP(1), BitPackedExponent(0x0000000100000001)),
                ])
                p3 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000000010101)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000001010100)),
                    Term(IntModP(1), BitPackedExponent(0x0000000101010000)),
                    Term(IntModP(1), BitPackedExponent(0x0000000101000001)),
                    Term(IntModP(1), BitPackedExponent(0x0000000100000101)),
                ])
                p4 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000001010101)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000101010100)),
                    Term(IntModP(1), BitPackedExponent(0x0000000101010001)),
                    Term(IntModP(1), BitPackedExponent(0x0000000101000101)),
                    Term(IntModP(1), BitPackedExponent(0x0000000100010101)),
                ])
                p5 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000101010101)), 
                    Term(IntModP(get_modulus() - 1), BitPackedExponent(0x0000000000000000)),
                ])
            end
            for i in 0:9
                basis = grobner_basis([p1, p2, p3, p4, p5])
                println("Iteration $i complete")
            end
        elseif n == 6
            if mode == 0
                # cyclic 6
                println("Julia generic vec exponent cyclic 6")
                p1 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 0, 0, 0, 0, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 0, 0, 0, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 1, 0, 0, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 0, 1, 0, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 0, 0, 1, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 0, 0, 0, 1])),
                ])
                p2 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 0, 0, 0, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 1, 0, 0, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 1, 1, 0, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 0, 1, 1, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 0, 0, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 0, 0, 0, 0, 1])),
                ])
                p3 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 1, 0, 0, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 1, 1, 0, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 1, 1, 1, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 0, 1, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 0, 0, 0, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 1, 0, 0, 0, 1])),
                ])
                p4 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 1, 1, 0, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 1, 1, 1, 0])),
                    Term(IntModP(1), VecExponent([0, 0, 1, 1, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 0, 0, 1, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 1, 0, 0, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 1, 1, 0, 0, 1])),
                ])
                p5 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 1, 1, 1, 0])), 
                    Term(IntModP(1), VecExponent([0, 1, 1, 1, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 0, 1, 1, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 1, 0, 1, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 1, 1, 0, 1, 1])),
                    Term(IntModP(1), VecExponent([1, 1, 1, 1, 0, 1])),
                ])
                p6 = Polynomial([
                    Term(IntModP(1), VecExponent([1, 1, 1, 1, 1, 1])), 
                    Term(IntModP(get_modulus() - 1), VecExponent([0, 0, 0, 0, 0, 0])),
                ])
            else
                # cyclic 6 with bit packed exponents
                println("Julia generic bit packed exponent cyclic 6")
                p1 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000000000001)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000000000100)),
                    Term(IntModP(1), BitPackedExponent(0x0000000000010000)),
                    Term(IntModP(1), BitPackedExponent(0x0000000001000000)),
                    Term(IntModP(1), BitPackedExponent(0x0000000100000000)),
                    Term(IntModP(1), BitPackedExponent(0x0000010000000000)),
                ])
                p2 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000000000101)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000000010100)),
                    Term(IntModP(1), BitPackedExponent(0x0000000001010000)),
                    Term(IntModP(1), BitPackedExponent(0x0000000101000000)),
                    Term(IntModP(1), BitPackedExponent(0x0000010100000000)),
                    Term(IntModP(1), BitPackedExponent(0x0000010000000001)),
                ])
                p3 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000000010101)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000001010100)),
                    Term(IntModP(1), BitPackedExponent(0x0000000101010000)),
                    Term(IntModP(1), BitPackedExponent(0x0000010101000000)),
                    Term(IntModP(1), BitPackedExponent(0x0000010100000001)),
                    Term(IntModP(1), BitPackedExponent(0x0000010000000101)),
                ])
                p4 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000001010101)), 
                    Term(IntModP(1), BitPackedExponent(0x0000000101010100)),
                    Term(IntModP(1), BitPackedExponent(0x0000010101010000)),
                    Term(IntModP(1), BitPackedExponent(0x0000010101000001)),
                    Term(IntModP(1), BitPackedExponent(0x0000010100000101)),
                    Term(IntModP(1), BitPackedExponent(0x0000010000010101)),
                ])
                p5 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000000101010101)), 
                    Term(IntModP(1), BitPackedExponent(0x0000010101010100)),
                    Term(IntModP(1), BitPackedExponent(0x0000010101010001)),
                    Term(IntModP(1), BitPackedExponent(0x0000010101000101)),
                    Term(IntModP(1), BitPackedExponent(0x0000010100010101)),
                    Term(IntModP(1), BitPackedExponent(0x0000010001010101)),
                ])
                p6 = Polynomial([
                    Term(IntModP(1), BitPackedExponent(0x0000010101010101)), 
                    Term(IntModP(get_modulus() - 1), BitPackedExponent(0x0000000000000000)),
                ])
            end
            for i in 0:9
                basis = grobner_basis([p1, p2, p3, p4, p5, p6])
                println("Iteration $i complete")
            end
        end
        #println("Basis size: ", length(basis))
        #= for b in basis
            for t in b.terms
                print("(", t.coefficient.value, " * ", t.exponents.exps, ") ")
            end
            println()
        end =#
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
     n = 4 # Default to cyclic 4
    if length(ARGS) > 0
        n = parse(Int, ARGS[1])
    end
    mode = 0
    if length(ARGS) > 1
        mode = parse(Int, ARGS[2])
    end
    GenGrobner.main(n, mode)
end
