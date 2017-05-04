using Base.Test, Mathieu
#PRINTVALS = true
import Base.product
import GSL: sf_mathieu_a, sf_mathieu_b, sf_mathieu_ce, sf_mathieu_se
gsl_mathieu_a(n,q) = sf_mathieu_a(n,q).val
gsl_mathieu_b(n,q) = sf_mathieu_b(n,q).val
gsl_mathieu_ce(n::Integer,q,z::Number) = sf_mathieu_ce(n,q,z).val
gsl_mathieu_se(n::Integer,q,z::Number) = sf_mathieu_se(n,q,z).val
gsl_mathieu_ce(n::AbstractVector,q,z::AbstractVector) = mapreduce(n->gsl_mathieu_ce.(n,q,z),hcat,n)
gsl_mathieu_se(n::AbstractVector,q,z::AbstractVector) = mapreduce(n->gsl_mathieu_se.(n,q,z),hcat,n)

# Order n and parameter q settings
nn = 0:20
qq = 0:2φ:30
zz = linspace(0,4φ,10)
cases = ((q->Mathieu.ap(nn,q), q->gsl_mathieu_a.(2nn,q), "MathieuCharA_Periodic"),
         (q->Mathieu.aa(nn,q), q->gsl_mathieu_a.(2nn+1,q), "MathieuCharA_Antiperiodic"),
         (q->Mathieu.ba(nn,q), q->gsl_mathieu_b.(2nn+1,q), "MathieuCharB_Antiperiodic"),
         (q->Mathieu.bp(nn,q), q->gsl_mathieu_b.(2nn+2,q), "MathieuCharB_Periodic"),
         (q->Mathieu.cep(nn,q,zz), q->gsl_mathieu_ce(2nn,q,zz), "MathieuC_Periodic"),
         (q->Mathieu.cea(nn,q,zz), q->gsl_mathieu_ce(2nn+1,q,zz), "MathieuC_Antiperiodic"),
         (q->Mathieu.sea(nn,q,zz), q->gsl_mathieu_se(2nn+1,q,zz), "MathieuS_Antiperiodic"),
         (q->Mathieu.sep(nn,q,zz), q->gsl_mathieu_se(2nn+2,q,zz), "MathieuS_Periodic"),
         (q->Mathieu.cep_prime(nn,q,zz), q->gsl_mathieu_ce(2nn,q,zz), "MathieuCPrime_Periodic"),
         (q->Mathieu.cea_prime(nn,q,zz), q->gsl_mathieu_ce(2nn+1,q,zz), "MathieuCPrime_Antiperiodic"),
         (q->Mathieu.sea_prime(nn,q,zz), q->gsl_mathieu_se(2nn+1,q,zz), "MathieuSPrime_Antiperiodic"),
         (q->Mathieu.sep_prime(nn,q,zz), q->gsl_mathieu_se(2nn+2,q,zz), "MathieuSPrime_Periodic"))

refdir = joinpath(dirname(@__FILE__),"ref")

for (fjl,fgsl,file) in cases
    print(string("Testing ",file,"... "))
    # Read Mathematica reference values
    arr_ref = readcsv(joinpath(refdir,string(file,".csv")),Float64)
    # Calculate with GSL
    #arr_gsl = mapreduce(fgsl,hcat,qq)
    # Calculate corresponding values with Mathieu.jl
    arr_jul = mapreduce(fjl,hcat,qq)
    # Make sure the values are approximately equal
    @test arr_jul ≈ arr_ref
    a = abs.(arr_jul-arr_ref)./eps.(arr_jul)
    #@test arr_jul ≈ arr_gsl
    #a = abs.(arr_jul-arr_gsl)./eps.(arr_jul)
    println(@sprintf("mean error: %d eps, median: %d eps, max: %d eps.",mean(a),median(a),maximum(a)))
end

## Order n and parameter q settings
#nn = 0:18
#qq = 0:2φ:200
#q_0to10    = 1:4
#q_10to40   = 5:13
#q_40to100  = 14:31
#q_100to200 = 32:62
#
#charval_testfiles = #["MathieuCharA_Periodic","MathieuCharA_Antiperiodic","MathieuCharB_Antiperiodic","MathieuCharB_Periodic"]
#charval_funs = [q->Mathieu.ap(nn,q),q->Mathieu.aa(nn,q),q->Mathieu.ba(nn,q),q->Mathieu.bp(nn,q)]
#charval_max1  = [6,   21,   6,   7]
#charval_max2  = [53,  21,   37,  16]
#charval_max3  = [76,  783,  149, 446]
#charval_max4  = [416, 1837, 83,  419]
#charval_mean1 = [1.8, 2.0,  1.7, 1.4]
#charval_mean2 = [2.9, 2.7,  2.6, 2.4]
#charval_mean3 = [3.1, 5.1,  3.2, 3.7]
#charval_mean4 = [4.2, 10.2, 2.8, 3.7]
#
#for i = 1:length(charval_testfiles)
#    println(string("Testing ",charval_testfiles[i]))
#    # Read Mathematica reference values
#    arr_ref = #parse.([BigFloat],readcsv(joinpath(dirname(@__FILE__),string(charval_testfiles[i],".csv")),String))
#    # Calculate corresponding values
#    arr_jul = mapreduce(charval_funs[i],hcat,qq).'
#    # Compute difference in multiples of eps
#    a = Float64.(abs.(arr_jul-arr_ref)./eps.(arr_jul))
#    arr_jul ≈ arr_ref
#    # Test that this relative difference is below some threshold
#    # Maximum differences are pretty bad
#    @test maximum(a[q_0to10,:])    < charval_max1[i]
#    @test maximum(a[q_10to40,:])   < charval_max2[i]
#    @test maximum(a[q_40to100,:])  < charval_max3[i]
#    @test maximum(a[q_100to200,:]) < charval_max4[i]
#    # Mean differences are better, although within 2 eps would be preferable
#    @test mean(a[q_0to10,:])       < charval_mean1[i]
#    @test mean(a[q_10to40,:])      < charval_mean2[i]
#    @test mean(a[q_40to100,:])     < charval_mean3[i]
#    @test mean(a[q_100to200,:])    < charval_mean4[i]
#
#    if PRINTVALS # See the numbers!
#        println("Maximum relative difference with Mathematica (in eps):")
#        println(string("q∈[0,10]:    ",maximum(a[q_0to10,:])))
#        println(string("q∈[10,40]:   ",maximum(a[q_10to40,:])))
#        println(string("q∈[40,100]:  ",maximum(a[q_40to100,:])))
#        println(string("q∈[100,200]: ",maximum(a[q_100to200,:])))
#        println("Mean relative difference with Mathematica (in eps):")
#        println(string("q∈[0,10]:    ",mean(a[q_0to10,:])))
#        println(string("q∈[10,40]:   ",mean(a[q_10to40,:])))
#        println(string("q∈[40,100]:  ",mean(a[q_40to100,:])))
#        println(string("q∈[100,200]: ",mean(a[q_100to200,:])))
#        println()
#    end
#end
