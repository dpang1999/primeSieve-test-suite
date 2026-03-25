import .IExponent

struct BitPackedExponent <: IExponent
    exps::UInt64
end

# Assuming 8 bits per exponent, up to 8 variables
function degree(a::BitPackedExponent)
    d = 0
    t = a.exps
    while t > 0
        d += t & 0xFF
        t >>= 8
    end
    return d
end

function lex_compare(a::BitPackedExponent, b::BitPackedExponent)
    # Simple uint compare is not lex compare, this is a placeholder
    return cmp(a.exps, b.exps)
end

function can_reduce(a::BitPackedExponent, b::BitPackedExponent)
    # Check if each 8-bit chunk in a >= b
    t_a = a.exps
    t_b = b.exps
    for i in 1:8
        if (t_a & 0xFF) < (t_b & 0xFF)
            return false
        end
        t_a >>= 8
        t_b >>= 8
    end
    return true
end

function get_lcm(a::BitPackedExponent, b::BitPackedExponent)
    # max of each 8-bit chunk
    res = UInt64(0)
    t_a = a.exps
    t_b = b.exps
    for i in 0:7
        va = t_a & 0xFF
        vb = t_b & 0xFF
        m = max(va, vb)
        res |= (m << (i * 8))
        t_a >>= 8
        t_b >>= 8
    end
    return BitPackedExponent(res)
end

Base.:+(a::BitPackedExponent, b::BitPackedExponent) = BitPackedExponent(a.exps + b.exps)
Base.:-(a::BitPackedExponent, b::BitPackedExponent) = BitPackedExponent(a.exps - b.exps)
