abstract type IField end

import Base: +, -, *, /, zero, one

# Arithmetic operations
+(a::IField, b::IField) = error("Not implemented")
-(a::IField, b::IField) = error("Not implemented")
*(a::IField, b::IField) = error("Not implemented")
/(a::IField, b::IField) = error("Not implemented")

# Identity elements
zero(::Type{<:IField}) = error("Not implemented")
one(::Type{<:IField}) = error("Not implemented")

# Predicates
is_zero(a::IField) = error("Not implemented")
is_one(a::IField) = error("Not implemented")

# Coercions (matching Rust coercions)
coerce_to_f64(a::IField) = error("Not implemented")
coerce_from_int(::Type{<:IField}, value::Int32) = error("Not implemented")
coerce(::Type{<:IField}, value::Float64) = error("Not implemented")

