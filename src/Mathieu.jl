module Mathieu
using LinearAlgebra
using LinearAlgebra: BlasReal
using Base.Iterators: product, filter

const Order = Union{Integer, AbstractVector{T} where T<:Integer}

include("integerorder_mat.jl")
include("integerorder_evals.jl")
include("integerorder_coeffs.jl")
include("integerorder_funs.jl")
include("integerorder_derivs.jl")
include("integerorder_special.jl")
include("utils.jl")

end # module
