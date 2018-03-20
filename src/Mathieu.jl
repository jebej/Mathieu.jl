module Mathieu
using Base: product, Iterators.filter, LinAlg.BlasReal
using Compat

@compat const Order{T<:Integer}  = Union{T,AbstractVector{T}}

include("integerorder_mat.jl")
include("integerorder_evals.jl")
include("integerorder_coeffs.jl")
include("integerorder_funs.jl")
include("integerorder_derivs.jl")
include("integerorder_special.jl")
include("utils.jl")

end # module
