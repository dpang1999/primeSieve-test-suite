import .IField

struct SingleField <: IField
    value::Float32
end
Base.:+(a::SingleField, b::SingleField) = SingleField(a.value + b.value)
Base.:-(a::SingleField, b::SingleField) = SingleField(a.value - b.value)
Base.:*(a::SingleField, b::SingleField) = SingleField(a.value * b.value)
Base.:/(a::SingleField, b::SingleField) = SingleField(a.value / b.value)

Base.zero(::Type{SingleField}) = SingleField(0.0f0)
Base.one(::Type{SingleField}) = SingleField(1.0f0)
is_zero(a::SingleField) = abs(a.value - 0.0f0) < 1e-2
is_one(a::SingleField) = abs(a.value - 1.0f0) < 1e-2

coerce_to_f64(a::SingleField) = Float64(a.value)
coerce_from_int(::Type{SingleField}, value::Int32) = SingleField(Float32(value))
coerce(::Type{SingleField}, value::Float64) = SingleField(Float32(value))

