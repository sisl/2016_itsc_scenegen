run(`clear`)

include("../full_scene.jl")
tic()
dsetexp = JLD.load("../output/dsetexport.jld", "dsetexp")
toc()

println("nscenes: ", length(dsetexp.scenes))
println("nvehstates: ", length(dsetexp.states))
println("ave veh / scene: ", length(dsetexp.states)/length(dsetexp.scenes))

tic()
dset = reconstruct_dataset(dsetexp, create_core_factors(), Scene())
toc()

tic()
roadscenes = get_roadscene_dset(dset)
toc()

# opt_res = cyclic_coordinate_ascent_parallel(JointBNChainSceneGenerator, roadscenes, CV_nfolds = 5, CV_rounds = 10)
# JLD.save("../output/BN_opt_res.jld", "opt_res", opt_res)

# ModelOptimizationResults{JointBNChainSceneGenerator}(Dict{Symbol,Any}(:d_front=>10,:yaw=>30,:d_cl=>15,:speed=>6,:v_front=>25,:v_rear=>25,:d_rear=>30),-228679.04182735403,5990.9839734724665,-228679.04182735403,4,5,10)

# train a final model based on opt_res:
nbins = Dict{Symbol,Int}(:d_front=>10,:yaw=>30,:d_cl=>15,:speed=>25,:v_front=>10,:v_rear=>10,:d_rear=>15)
tic()
BN = train_jointbnchainscenegenerator(roadscenes, nbins)
toc()

JLD.save(joinpath(SAVE_DIR, "BN.jld"), "discmap", BN.discmap, "varindeces", BN.varindeces, "ordering", BN.ordering)
write_file( BN.net, joinpath(SAVE_DIR, "BN_net.xdsl") )

dataset_size = 10000

tic()
dset_BN = generate_dset(BN, dset, dataset_size)
toc()
JLD.save(joinpath(SAVE_DIR, "dset_BN.jld"), "dset", dset_BN, "time", string(now()))

println("DONE")
