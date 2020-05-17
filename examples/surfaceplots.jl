# Plot of various Mathieu functions, inspired by https://dlmf.nist.gov/28.3
using Mathieu, PyPlot

z = LinRange(0,2π,101)
q = LinRange(0,10,101)

# ce₀(z,q)
figure(figsize=(10,8));
plot_surface(q,z,mapreduce(q->Mathieu.cep(0,q,z),hcat,q),cmap="coolwarm")
title("ce₀(z,q)"); xlabel("q"); ylabel("z");
zlim([0,1.5]); gca(projection="3d").view_init(azim=-30)

# ce₁(z,q)
figure(figsize=(10,8));
plot_surface(q,z,mapreduce(q->Mathieu.cea(0,q,z),hcat,q),cmap="coolwarm")
title("ce₁(z,q)"); xlabel("q"); ylabel("z");
zlim([-1.2,1.2]); gca(projection="3d").view_init(azim=-30)

# se₁(z,q)
figure(figsize=(10,8));
plot_surface(q,z,mapreduce(q->Mathieu.sea(0,q,z),hcat,q),cmap="coolwarm")
title("se₁(z,q)"); xlabel("q"); ylabel("z");
zlim([-1.5,1.5]); gca(projection="3d").view_init(azim=-30)

# se₂(z,q)
figure(figsize=(10,8));
plot_surface(q,z,mapreduce(q->Mathieu.sep(0,q,z),hcat,q),cmap="coolwarm")
title("se₂(z,q)"); xlabel("q"); ylabel("z");
zlim([-1.2,1.2]); gca(projection="3d").view_init(azim=-30)

# ce₂(z,q)
figure(figsize=(10,8));
plot_surface(q,z,mapreduce(q->Mathieu.cep(1,q,z),hcat,q),cmap="coolwarm")
title("ce₂(z,q)"); xlabel("q"); ylabel("z");
zlim([-1,1.1]); gca(projection="3d").view_init(azim=-30)
