using Base.Test, Mathieu
import Base.product

PRINTVALS = true

# Order n and parameter q settings
nn = 0:18
qq = 0:2φ:200
q_0to10    = 1:4
q_10to40   = 5:13
q_40to100  = 14:31
q_100to200 = 32:62

charval_testfiles = ["MathieuCharA_Periodic","MathieuCharA_Antiperiodic","MathieuCharB_Antiperiodic","MathieuCharB_Periodic"]
charval_funs = [q->Mathieu.ap(nn,q),q->Mathieu.aa(nn,q),q->Mathieu.ba(nn,q),q->Mathieu.bp(nn,q)]
charval_max1  = [6,   21,   6,   7]
charval_max2  = [53,  21,   37,  16]
charval_max3  = [76,  783,  149, 446]
charval_max4  = [416, 1837, 83,  419]
charval_mean1 = [1.8, 2.0,  1.7, 1.4]
charval_mean2 = [2.9, 2.7,  2.6, 2.4]
charval_mean3 = [3.1, 5.1,  3.2, 3.7]
charval_mean4 = [4.2, 10.2, 2.8, 3.7]

for i = 1:length(charval_testfiles)
    println(string("Testing ",charval_testfiles[i]))
    # Read Mathematica reference values
    arr_ref = parse.([BigFloat],readcsv(joinpath(dirname(@__FILE__),string(charval_testfiles[i],".csv")),String))
    # Calculate corresponding values
    arr_jul = mapreduce(charval_funs[i],hcat,qq).'
    # Compute difference in multiples of eps
    a = Float64.(abs.(arr_jul-arr_ref)./eps.(arr_jul))

    # Test that this relative difference is below some threshold
    # Maximum differences are pretty bad
    @test maximum(a[q_0to10,:])    < charval_max1[i]
    @test maximum(a[q_10to40,:])   < charval_max2[i]
    @test maximum(a[q_40to100,:])  < charval_max3[i]
    @test maximum(a[q_100to200,:]) < charval_max4[i]
    # Mean differences are better, although within 2 eps would be preferable
    @test mean(a[q_0to10,:])       < charval_mean1[i]
    @test mean(a[q_10to40,:])      < charval_mean2[i]
    @test mean(a[q_40to100,:])     < charval_mean3[i]
    @test mean(a[q_100to200,:])    < charval_mean4[i]

    if PRINTVALS # See the numbers!
        println("Maximum relative difference with Mathematica (in eps):")
        println(string("q∈[0,10]:    ",maximum(a[q_0to10,:])))
        println(string("q∈[10,40]:   ",maximum(a[q_10to40,:])))
        println(string("q∈[40,100]:  ",maximum(a[q_40to100,:])))
        println(string("q∈[100,200]: ",maximum(a[q_100to200,:])))
        println("Mean relative difference with Mathematica (in eps):")
        println(string("q∈[0,10]:    ",mean(a[q_0to10,:])))
        println(string("q∈[10,40]:   ",mean(a[q_10to40,:])))
        println(string("q∈[40,100]:  ",mean(a[q_40to100,:])))
        println(string("q∈[100,200]: ",mean(a[q_100to200,:])))
        println()
    end
end
