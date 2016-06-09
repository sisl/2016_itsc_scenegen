include("../full_scene.jl")

tic()
dsetexp = pull_exportable_subscene_dataset()
toc()

println("nscenes: ", length(dsetexp.scenes))
println("nvehstates: ", length(dsetexp.states))
println("ave veh / scene: ", length(dsetexp.states)/length(dsetexp.scenes))

println("SAVING")
tic()
JLD.save("../output/dsetexport.jld", "dsetexp", dsetexp, "time", now())
toc()

tic()
dset = reconstruct_dataset(dsetexp, create_core_factors())
toc()

"done"