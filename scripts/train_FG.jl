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

params = StochasticGradientAscentParams()
params.n_samples_monte_carlo_integration = 100
params.n_samples_monte_carlo_pseudolikelihood = 100
params.n_epochs = 1
params.max_n_batches = 100
params.learning_rate = 0.001
params.learning_rate_decay = 0.97
params.regularization = 0.001

params.batch_size = 20
params.batch_rng = MersenneTwister(1)
params.mc_pseudolikelihood_rng = MersenneTwister(1)

params.print_to_stdout = true
params.print_to_stdout_every = 1
params.print_to_file = joinpath(params.save_dir, "log.txt")
params.print_to_file_every = 1
params.save_every = 1

reset_weights!(dset)

println("weights before: ")
for ϕ in dset.factors
    print(ϕ.name, ": [")
    for w in ϕ.weights
        @printf("%10.6f  ", w)
    end
    println("")
end

stochastic_gradient_ascent!(dset, params)

println("weights after: ")
for ϕ in dset.factors
    print(ϕ.name, ": [")
    for w in ϕ.weights
        @printf("%10.6f  ", w)
    end
    println("")
end

########################################
# Generate a Dataset
########################################

dataset_size = 10000

model_file = get_most_resent_file("../output")
println("model file: ", model_file)
saved_model = JLD.load(model_file)
println("save time: ", saved_model["time"])

factors = create_core_factors()
set_weights!(factors, saved_model["weights"])
println("DONE")

transition = MvNormal([2.0, 0.5, 1.0, 0.01]) # s, t, v, theta standard deviations
nsteps_burnin = 1000
FG = FactorGraphSceneGenerator(factors, transition, nsteps_burnin)

tic()
dset_FG = generate_dset(FG, dset, dataset_size)
toc()
JLD.save(joinpath(SAVE_DIR, "dset_FG.jld"), "dset", dset_FG, "time", string(now()))

println("DONE")