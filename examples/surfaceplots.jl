# Plot of various Mathieu functions, inspired by https://dlmf.nist.gov/28.3
using Mathieu, PyPlot

z = LinRange(0,2π,101)
q = LinRange(0,10,101)

## ce₀(z,q)
f = figure(num=1, figsize=(10,8)); f.clf();
plot_surface(q,z,mapreduce(q->Mathieu.ce(0,q,z),hcat,q),cmap="coolwarm")
title("ce₀(z,q)"); xlabel("q"); ylabel("z"); zlim(0,1.5);
gca(projection="3d").view_init(azim=-30); display(f);

## ce₁(z,q)
f = figure(num=2, figsize=(10,8)); f.clf();
plot_surface(q,z,mapreduce(q->Mathieu.ce(1,q,z),hcat,q),cmap="coolwarm")
title("ce₁(z,q)"); xlabel("q"); ylabel("z"); zlim(-1.2,1.2);
gca(projection="3d").view_init(azim=-30); display(f);

## se₁(z,q)
f = figure(num=3, figsize=(10,8)); f.clf();
plot_surface(q,z,mapreduce(q->Mathieu.se(1,q,z),hcat,q),cmap="coolwarm")
title("se₁(z,q)"); xlabel("q"); ylabel("z"); zlim(-1.5,1.5);
gca(projection="3d").view_init(azim=-30); display(f);

## se₂(z,q)
f = figure(num=4, figsize=(10,8)); f.clf();
plot_surface(q,z,mapreduce(q->Mathieu.se(2,q,z),hcat,q),cmap="coolwarm")
title("se₂(z,q)"); xlabel("q"); ylabel("z"); zlim(-1.2,1.2);
gca(projection="3d").view_init(azim=-30); display(f);

## ce₂(z,q)
f = figure(num=5, figsize=(10,8)); f.clf();
plot_surface(q,z,mapreduce(q->Mathieu.ce(2,q,z),hcat,q),cmap="coolwarm")
title("ce₂(z,q)"); xlabel("q"); ylabel("z"); zlim(-1,1.1);
gca(projection="3d").view_init(azim=-30); display(f);
