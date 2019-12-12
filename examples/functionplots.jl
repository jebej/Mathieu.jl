# Plot of various Mathieu functions, inspired by https://dlmf.nist.gov/28.3
using Mathieu, PyPlot
const MPL = PyPlot.matplotlib

# get the same line colors
PyPlot.PyDict(MPL."rcParams")["axes.prop_cycle"] =
    MPL.cycler(color=["#2ca02c","#d62728","#1f77b4","#ff7f0e"])

# plot parameters
z = LinRange(0,π/2,101)
n = 0:3

for q in [1,10]
    f,ax = subplots(2,2,figsize=(10,8));
    f.suptitle("Mathieu Functions for \$q=$q\$");

    ax[1].set_title("Even π-Periodic Solutions");
    ax[1].plot(z, Mathieu.cep(n,q,z));
    ax[1].legend(string.("\$ce_",2n,"\$"));

    ax[2].set_title("Even π-Antiperiodic Solutions");
    ax[2].plot(z, Mathieu.cea(n,q,z));
    ax[2].legend(string.("\$ce_",2n.+1,"\$"));

    ax[3].set_title("Odd π-Antiperiodic Solutions");
    ax[3].plot(z, Mathieu.sea(n,q,z));
    ax[3].legend(string.("\$se_",2n.+1,"\$"));

    ax[4].set_title("Odd π-Periodic Solutions");
    ax[4].plot(z, Mathieu.sep(n,q,z));
    ax[4].legend(string.("\$se_",2n.+2,"\$"));
end
