import .IField

struct SingleField <: IField
    value::Float32
end
Base.:+(a::SingleField, b::SingleField) = SingleField(a.value + b.value)
Base.:-(a::SingleField, b::SingleField) = SingleField(a.value - b.value)
Base.:*(a::SingleField, b::SingleField) = SingleField(a.value * b.value)
Base.:/(a::SingleField, b::SingleField) = SingleField(a.value / b.value)
is_zero(a::SingleField) = abs(a.value - 0.0f0) < 1e-2
