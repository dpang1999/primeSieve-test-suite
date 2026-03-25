import .IField

struct DoubleField <: IField
    value::Float64
end
Base.:+(a::DoubleField, b::DoubleField) = DoubleField(a.value + b.value)
Base.:-(a::DoubleField, b::DoubleField) = DoubleField(a.value - b.value)
Base.:*(a::DoubleField, b::DoubleField) = DoubleField(a.value * b.value)
Base.:/(a::DoubleField, b::DoubleField) = DoubleField(a.value / b.value)
is_zero(a::DoubleField) = abs(a.value - 0.0) < 1e-2
