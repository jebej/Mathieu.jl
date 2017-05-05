using Mathieu,Plots
plotlyjs()
z = linspace(0,2π,101)
q = linspace(0,10,101)

# ce_0(z,q)
surface(q,z,mapreduce(q->Mathieu.cep(0,q,z),hcat,q),zlim=(0,1.5))

# ce_1(z,q)
surface(q,z,mapreduce(q->Mathieu.cea(0,q,z),hcat,q),zlim=(-1.2,1.2))

# se_1(z,q)
surface(q,z,mapreduce(q->Mathieu.sea(0,q,z),hcat,q),zlim=(-1.5,1.5))

# se_2(z,q)
surface(q,z,mapreduce(q->Mathieu.sep(0,q,z),hcat,q),zlim=(-1.2,1.2))

# ce_2(z,q)
surface(q,z,mapreduce(q->Mathieu.cep(1,q,z),hcat,q),zlim=(-1,1.1))
