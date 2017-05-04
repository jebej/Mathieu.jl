module Mathieu
import Base: product, LinAlg.BlasReal

typealias Order{T<:Integer} Union{T,AbstractVector{T}}

include("integerorder_mat.jl")
include("integerorder_evals.jl")
include("integerorder_funs.jl")
include("utils.jl")

end # module
