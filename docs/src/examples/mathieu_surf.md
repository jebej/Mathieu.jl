```@setup mathieu_plots
using Mathieu, PyPlot
!isdir("img") && mkdir("img")
include(joinpath(@__DIR__,"..","..","..","examples","surfaceplots.jl"))
savefig(joinpath("img","surf_5.svg"),transparent=true); close()
savefig(joinpath("img","surf_4.svg"),transparent=true); close()
savefig(joinpath("img","surf_3.svg"),transparent=true); close()
savefig(joinpath("img","surf_2.svg"),transparent=true); close()
savefig(joinpath("img","surf_1.svg"),transparent=true); close()
```
# Mathieu Surface Plots

We plot Mathieu functions for q  between 0 and 10 in order to show the periodic and
anti-periodic properties of the sine-elliptic and cosine-elliptic solutions.

The plots were inspired by the NIST Digital Library of Mathematical Functions.

![surf_1](img/surf_1.svg)
![surf_2](img/surf_2.svg)
![surf_3](img/surf_3.svg)
![surf_4](img/surf_4.svg)
![surf_5](img/surf_5.svg)

The code used to generate these plots is also available in the example folder of
the package.

````@eval
using Markdown: parse
parse("""
```julia
$(read(joinpath(@__DIR__,"..","..","..","examples","surfaceplots.jl"),String))
```
""")
````
