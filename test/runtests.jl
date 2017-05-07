using Base.Test, Mathieu

## Test vs Mathematica
n1 = 0:24
q1 = 0:2φ:80
z1 = linspace(0,4φ,10)
cases = ((q->Mathieu.ap(n1,q),           "MathieuCharA_Periodic"),
         (q->Mathieu.aa(n1,q),           "MathieuCharA_Antiperiodic"),
         (q->Mathieu.ba(n1,q),           "MathieuCharB_Antiperiodic"),
         (q->Mathieu.bp(n1,q),           "MathieuCharB_Periodic"),
         (q->Mathieu.cep(n1,q,z1),       "MathieuC_Periodic"),
         (q->Mathieu.cea(n1,q,z1),       "MathieuC_Antiperiodic"),
         (q->Mathieu.sea(n1,q,z1),       "MathieuS_Antiperiodic"),
         (q->Mathieu.sep(n1,q,z1),       "MathieuS_Periodic"),
         (q->Mathieu.cep_prime(n1,q,z1), "MathieuCPrime_Periodic"),
         (q->Mathieu.cea_prime(n1,q,z1), "MathieuCPrime_Antiperiodic"),
         (q->Mathieu.sea_prime(n1,q,z1), "MathieuSPrime_Antiperiodic"),
         (q->Mathieu.sep_prime(n1,q,z1), "MathieuSPrime_Periodic"),
         )
refdir = joinpath(dirname(@__FILE__),"ref")
for (fjl,file) in cases
    println(string("Comparing ",file," with Mathematica:"))
    # Read Mathematica reference values
    arr_ref = readcsv(joinpath(refdir,string(file,".csv")),Float64)
    # Calculate corresponding values with Mathieu.jl
    arr_jul = mapreduce(fjl,hcat,q1)
    # Make sure the values are approximately equal
    @test arr_jul ≈ arr_ref
    a = abs.(arr_jul-arr_ref)./eps.(arr_jul)
    println(@sprintf("  mean error: %d eps, median: %d eps, max: %d eps.",mean(a),median(a),maximum(a)))
end

println(string("Comparing ","Int_Per_PerPrime"," with Mathematica:"))
# Read Mathematica reference values
arr_ref = readcsv(joinpath(refdir,string("Int_Per_PerPrime",".csv")),Float64)
# Calculate corresponding values with Mathieu.jl
arr_jul = Mathieu.int_per_perprime(n1,n1,223/100,[0.5,3])
# Make sure the values are approximately equal
@test arr_jul ≈ arr_ref[1:25,1:25]
a = abs.(arr_jul-arr_ref[1:25,1:25])./eps.(arr_jul)
println(@sprintf("  mean error: %d eps, median: %d eps, max: %d eps.",mean(a),median(a),maximum(a)))

## Test q change of sign
n2 = [0:20,7:21,10:35,12:46,50:100,100:150,150:200]
q2 = linspace(0,100,202) #0:2φ:30
z2 = linspace(0,4φ,10)
for nn in n2
    println("Testing change of sign for n = $nn:")
    N2 = maximum(nn)+12 # Make N different to hit different eigendecompositions
    print("  ap")
    A1 = mapreduce(q->Mathieu.ap(nn,q),hcat,-q2)
    A2 = mapreduce(q->Mathieu.ap(nn,q,N2),hcat,q2)
    @test A1 ≈ A2
    print(", aa, ba")
    A1 = mapreduce(q->Mathieu.aa(nn,q),hcat,-q2)
    A2 = mapreduce(q->Mathieu.ba(nn,q,N2),hcat,q2)
    @test A1 ≈ A2
    print(", bp")
    A1 = mapreduce(q->Mathieu.bp(nn,q),hcat,-q2)
    A2 = mapreduce(q->Mathieu.bp(nn,q,N2),hcat,q2)
    @test A1 ≈ A2
    print(", cep")
    A1 = mapreduce(q->Mathieu.cep(nn,q,z2),vcat,-q2)
    A2 = mapreduce(q->Mathieu.cep(nn,q,π/2-z2,N2).*(-1).^nn.',vcat,q2)
    @test A1 ≈ A2
    print(", cea")
    A1 = mapreduce(q->Mathieu.cea(nn,q,z2),vcat,-q2)
    A2 = mapreduce(q->Mathieu.sea(nn,q,π/2-z2,N2).*(-1).^nn.',vcat,q2)
    @test A1 ≈ A2
    print(", sea")
    A1 = mapreduce(q->Mathieu.sea(nn,q,z2),vcat,-q2)
    A2 = mapreduce(q->Mathieu.cea(nn,q,π/2-z2,N2).*(-1).^nn.',vcat,q2)
    @test A1 ≈ A2
    println(", sep.")
    A1 = mapreduce(q->Mathieu.sep(nn,q,z2),vcat,-q2)
    A2 = mapreduce(q->Mathieu.sep(nn,q,π/2-z2,N2).*(-1).^nn.',vcat,q2)
    @test A1 ≈ A2
end
