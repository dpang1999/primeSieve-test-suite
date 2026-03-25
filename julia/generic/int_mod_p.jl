# IntModP implementation
import .IField

const GLOBAL_MODULUS = Ref{UInt64}(7)

function set_modulus(m::UInt64)
    GLOBAL_MODULUS[] = m
end

struct IntModP <: IField
    value::UInt64
end

IntModP(v::Integer) = IntModP(UInt64(v % GLOBAL_MODULUS[]))

Base.:+(a::IntModP, b::IntModP) = IntModP((a.value + b.value) % GLOBAL_MODULUS[])
Base.:-(a::IntModP, b::IntModP) = IntModP((GLOBAL_MODULUS[] + a.value - b.value) % GLOBAL_MODULUS[])
Base.:*(a::IntModP, b::IntModP) = IntModP((a.value * b.value) % GLOBAL_MODULUS[])
Base.:/(a::IntModP, b::IntModP) = a * mod_inverse(b)
Base.zero(::Type{IntModP}) = IntModP(0)
Base.one(::Type{IntModP}) = IntModP(1)
Base.zero(a::IntModP) = IntModP(0)
is_zero(a::IntModP) = a.value == 0
Base.show(io::IO, a::IntModP) = print(io, a.value)

function mod_inverse(a::IntModP)
    v::Int64 = a.value
    m::Int64 = GLOBAL_MODULUS[]
    m0::Int64 = m
    y::Int64 = 0
    x::Int64 = 1

    if m == 1
        return IntModP(0)
    end

    while v > 1
        q = div(v, m)
        t = m

        m = v % m
        v = t
        t = y

        y = x - q * y
        x = t
    end

    if x < 0
        x = x + m0
    end

    return IntModP(UInt64(x % m0))
end
