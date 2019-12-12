```@setup mathieu_plots
using Mathieu, PyPlot
!isdir("img") && mkdir("img")
include(joinpath(@__DIR__,"..","..","..","examples","functionplots.jl"))
savefig(joinpath("img","funs_10.svg"),transparent=true); close()
savefig(joinpath("img","funs_1.svg"),transparent=true); close()
```
# Mathieu Function Plots

We plot Mathieu functions for q = 1 and 10 in order to show the periodic and
anti-periodic properties of the sine-elliptic and cosine-elliptic solutions.

The plots were inspired by the NIST Digital Library of Mathematical Functions.

![function_1](img/funs_1.svg)

![function_10](img/funs_10.svg)

The code used to generate these plots is also available in the example folder of
the package.

````@eval
using Markdown: parse
parse("""
```julia
$(read(joinpath(@__DIR__,"..","..","..","examples","functionplots.jl"),String))
```
""")
````
