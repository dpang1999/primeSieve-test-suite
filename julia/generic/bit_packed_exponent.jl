import .IExponent

struct BitPackedExponent <: IExponent
    exps::UInt64
end

function BitPackedExponent(exponents::Vector{Int})
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
    
    return BitPackedExponent(packed)
end

function degree(a::BitPackedExponent)
    return (a.exps >> 48) & 0xFFFF
end

function lex_compare(a::BitPackedExponent, b::BitPackedExponent)
    self_exp = a.exps & 0x0000FFFFFFFFFFFF
    other_exp = b.exps & 0x0000FFFFFFFFFFFF
    return cmp(self_exp, other_exp)
end

function can_reduce(a::BitPackedExponent, b::BitPackedExponent)
    for i in 0:8:40
        self_exp = (a.exps >> i) & 0xFF
        divisor_exp = (b.exps >> i) & 0xFF
        if self_exp < divisor_exp
            return false
        end
    end
    return true
end

function get_lcm(a::BitPackedExponent, b::BitPackedExponent)
    self_exponents = a.exps & 0x0000FFFFFFFFFFFF
    other_exponents = b.exps & 0x0000FFFFFFFFFFFF
    
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
    return BitPackedExponent(lcm_exponents)
end


function Base.show(io::IO, a::BitPackedExponent)
    deg = (a.exps >> 48) & 0xFFFF
    print(io, "Degree: ", string(deg, base=16, pad=4), ", Exponents (hex): ")
    for i in 40:-8:0
        exp = (a.exps >> i) & 0xFF
        print(io, string(exp, base=16, pad=2), " ")
    end
end

Base.:+(a::BitPackedExponent, b::BitPackedExponent) = BitPackedExponent(a.exps + b.exps)
Base.:-(a::BitPackedExponent, b::BitPackedExponent) = BitPackedExponent(a.exps - b.exps)
Base.:(==)(a::BitPackedExponent, b::BitPackedExponent) = a.exps == b.exps
Base.hash(a::BitPackedExponent, h::UInt) = hash(a.exps, h)
