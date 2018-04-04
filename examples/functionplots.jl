using Mathieu, PyPlot

z = linspace(0,π/2,101)
n = 0:3

f,ax = subplots(2,2,figsize=(12,8)); f[:suptitle]("Mathieu Functions for \$q=1\$")
ax[1][:set_title]("Even π-Periodic Solutions")
ax[1][:plot](z,Mathieu.cep(n,1,z)); ax[1][:legend](string.("\$ce_",2n,"\$"));

ax[2][:set_title]("Even π-Antiperiodic Solutions")
ax[2][:plot](z,Mathieu.cea(n,1,z)); ax[2][:legend](string.("\$ce_",2n+1,"\$"));

ax[3][:set_title]("Odd π-Antiperiodic Solutions")
ax[3][:plot](z,Mathieu.sea(n,1,z)); ax[3][:legend](string.("\$se_",2n+1,"\$"));

ax[4][:set_title]("Odd π-Periodic Solutions")
ax[4][:plot](z,Mathieu.sep(n,1,z)); ax[4][:legend](string.("\$se_",2n+2,"\$"));
