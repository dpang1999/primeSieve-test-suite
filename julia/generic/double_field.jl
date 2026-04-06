import .IField

struct DoubleField <: IField
    value::Float64
end
Base.:+(a::DoubleField, b::DoubleField) = DoubleField(a.value + b.value)
Base.:-(a::DoubleField, b::DoubleField) = DoubleField(a.value - b.value)
Base.:*(a::DoubleField, b::DoubleField) = DoubleField(a.value * b.value)
Base.:/(a::DoubleField, b::DoubleField) = DoubleField(a.value / b.value)

Base.zero(::Type{DoubleField}) = DoubleField(0.0)
Base.one(::Type{DoubleField}) = DoubleField(1.0)
is_zero(a::DoubleField) = abs(a.value - 0.0) < 1e-2
is_one(a::DoubleField) = abs(a.value - 1.0) < 1e-2

coerce_to_f64(a::DoubleField) = a.value
coerce_from_int(::Type{DoubleField}, value::Int32) = DoubleField(Float64(value))
coerce(::Type{DoubleField}, value::Float64) = DoubleField(value)

