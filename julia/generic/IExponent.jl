abstract type IExponent end
import Base: +, -, ==, hash

+(a::IExponent, b::IExponent) = error("Not implemented")
-(a::IExponent, b::IExponent) = error("Not implemented")
==(a::IExponent, b::IExponent) = error("Not implemented")
hash(a::IExponent, h::UInt) = error("Not implemented")

# Custom operations for Gröbner Bases
get_lcm(a::IExponent, b::IExponent) = error("Not implemented")
degree(a::IExponent) = error("Not implemented")
lex_compare(a::IExponent, b::IExponent) = error("Not implemented")
can_reduce(a::IExponent, b::IExponent) = error("Not implemented")
