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

charval_testfiles = ["MathieuCharA_Periodic","MathieuCharA_Antiperiodic"]
charval_funs = [x->Mathieu.a_p(x[2],x[1]),x->Mathieu.a_a(x[2],x[1])]

charval_max1 = [6,21]
charval_max2 = [60,60]
charval_max3 = [80,800]
charval_max4 = [420,2000]
charval_mean1 = [1.8,2]
charval_mean2 = [2.9,3]
charval_mean3 = [3.1,6]
charval_mean4 = [4.2,11]

for i = 1:length(charval_testfiles)
    # Read Mathematica reference file
    arr_ref = read_testfile(string(charval_testfiles[i],".csv"))
    # Calculate corresponding values
    arr_jul = map(charval_funs[i],product(qq,nn))
    # Compute difference in multiples of eps
    a = Float64.(abs.(arr_jul-arr_ref)./eps.(arr_jul))
    # Test that this relative difference is below some threshold
    # Maximum differences are pretty bad
    @test maximum(a[q_0to10,:])    < charval_max1[i]
    @test maximum(a[q_10to40,:])   < charval_max2[i]
    @test maximum(a[q_40to100,:])  < charval_max3[i]
    @test maximum(a[q_100to200,:]) < charval_max4[i]
    # Mean differences are better, although within 2 eps would be preferable
    @test mean(a[q_0to10,:])    < charval_mean1[i]
    @test mean(a[q_10to40,:])   < charval_mean2[i]
    @test mean(a[q_40to100,:])  < charval_mean3[i]
    @test mean(a[q_100to200,:]) < charval_mean4[i]
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
