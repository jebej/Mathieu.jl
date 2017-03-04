using Base.Test, Mathieu
import Base.product

read_testfile(filename) = parse.([BigFloat],readcsv(joinpath(dirname(@__FILE__),filename),String))

# Order n and parameter q settings
nn = 0:18
qq = 0:2φ:200
q_0to10    = 1:4
q_10to40   = 5:13
q_40to100  = 14:31
q_100to200 = 32:62

charval_testfiles = ["MathieuCharA_Periodic"]

for file in charval_testfiles
    # Read Mathematica reference file
    arr_ref = read_testfile(string(file,".csv"))
    # Calculate corresponding values
    arr_jul = map(x->Mathieu.a_p(x[2],x[1]),product(qq,nn))
    # Compute difference in multiples of eps
    a = Float64.(abs.(arr_jul-arr_ref)./eps.(arr_jul))
    # Test that this relative difference is below some threshold
    # Maximum differences are pretty bad, between 5 and 500 eps
    @test maximum(a[q_0to10,:]) < 6
    @test maximum(a[q_10to40,:]) < 60
    @test maximum(a[q_40to100,:]) < 80
    @test maximum(a[q_100to200,:]) < 500
    # Mean differences are better, although within 2 eps would be preferable
    @test mean(a[q_0to10,:]) < 2
    @test mean(a[q_10to40,:]) < 3
    @test mean(a[q_40to100,:]) < 4
    @test mean(a[q_100to200,:]) < 5
end

# Comparison with GSL
#import GSL: sf_mathieu_a
#arr_gsl = map(x->sf_mathieu_a(x[2],x[1]).val,product(qq,2nn))
#b = Float64.(abs.(arr_gsl-arr_ref)./eps.(arr_gsl))
#maximum(b[q_0to10,:])
#maximum(b[q_10to40,:])
#maximum(b[q_40to100,:])
#maximum(b[q_100to200,:])
#mean(b[q_0to10,:])
#mean(b[q_10to40,:])
#mean(b[q_40to100,:])
#mean(b[q_100to200,:])
