import .IExponent

struct VecExponent <: IExponent
    exps::Vector{Int}
end

Base.:+(a::VecExponent, b::VecExponent) = VecExponent(a.exps .+ b.exps)
Base.:-(a::VecExponent, b::VecExponent) = VecExponent(a.exps .- b.exps)

function can_reduce(a::VecExponent, b::VecExponent)
    return all(a.exps .>= b.exps)
end

function get_lcm(a::VecExponent, b::VecExponent)
    return VecExponent(max.(a.exps, b.exps))
end

function lex_compare(a::VecExponent, b::VecExponent)
    return cmp(a.exps, b.exps)
end

function degree(a::VecExponent)
    return sum(a.exps)
end
