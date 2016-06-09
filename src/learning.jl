function get_variable_bounds(index_in_vehicle::Int, scene::Scene, structure::SceneStructure, vehicle_index::Int)
    if index_in_vehicle == IND_S

        veh = scene[vehicle_index]
        veh_index_fore = structure.index_fore[vehicle_index]
        @assert(veh_index_fore > 0)
        veh_fore = scene[veh_index_fore]
        veh_index_rear = structure.index_rear[vehicle_index]
        @assert(veh_index_rear > 0)
        veh_rear = scene[veh_index_rear]

        d_fore = get_headway_dist_between(veh, veh_fore)
        d_rear = get_headway_dist_between(veh_rear, veh)
        @assert(!isnan(d_fore))
        @assert(!isnan(d_rear))
        @assert(d_fore > 0.0)

        if !(d_rear > 0.0)
            println("veh: ")
            println(veh)
            println("veh rear:")
            println(veh_rear)
            println("d_rear: ", d_rear)

            active_laneid = veh_rear.state.posF.laneid
            curve = roadway.centerlines[active_laneid]
            s1 = curve_at(curve, veh_rear.state.posF.extind).s
            s2 = curve_at(curve, veh.state.posF.extind).s
            println("s1: ", s1)
            println("s2: ", s2)
            println("Δ2: ", s2 - s1)
            println("len: ", veh.length)
            println("dist: ", s2 - s1 - veh.length)
        end

        @assert(d_rear > 0.0)

        s = veh.state.posF.s
        (s - d_rear, s + d_fore)
    elseif index_in_vehicle == IND_T
        (BOUNDS_T_LO, BOUNDS_T_HI)
    elseif index_in_vehicle == IND_V
        (BOUNDS_V_LO, BOUNDS_V_HI)
    elseif index_in_vehicle == IND_ϕ
        (BOUNDS_ϕ_LO, BOUNDS_ϕ_HI)
    else
        error("this should never happen")
    end
end
function _get_logP_tilde_denom(lnP_tilde_denom_arr::Vector{Float64}, volume::Float64, n_samples_monte_carlo_integration::Int)
    P_tilde_denom = 0.0
    for v in lnP_tilde_denom_arr
        P_tilde_denom += exp(v)
    end
    log(P_tilde_denom * volume/n_samples_monte_carlo_integration)
end
function _update_lnP_tilde_denom_arr!(
    lnP_tilde_denom_arr::Vector{Float64},
    ϕ::Factor,
    variable_index::Int,
    lo::Float64,
    hi::Float64,
    nsamples::Int,
    roadway::Roadway
    )

    ###########################
    # v -> lnP(v | ~)
    #   vary v and compute the denominator term, ∫ᵥ p(v | ~) dv
    #   use Monte Carlo integration to get ∫ₐ P(a , x₋ⱼ) ≈ V/N ∑ₐ P(a , x₋ⱼ) where V is domain width

    v_original = get_variable(ϕ, variable_index)

    v = lo
    Δval = (hi - lo) / (nsamples-1)
    for k in 1 : nsamples
        set_variable!(ϕ, variable_index, v)
        extract_features!(ϕ, roadway)
        lnP_tilde_denom_arr[k] += factor_dot(ϕ)
        v += Δval
    end
    set_variable!(ϕ, variable_index, v_original)
end
function calc_pseudolikelihood(dset::SceneDataset;
    n_samples_monte_carlo_integration::Int=10,
    )

    return 1.0 # DEBUG

    retval = 0.0
    lnP_tilde_denom_arr_s = Array(Float64, n_samples_monte_carlo_integration)
    lnP_tilde_denom_arr_t = Array(Float64, n_samples_monte_carlo_integration)
    lnP_tilde_denom_arr_v = Array(Float64, n_samples_monte_carlo_integration)
    lnP_tilde_denom_arr_ϕ = Array(Float64, n_samples_monte_carlo_integration)

    M = length(dset.scenes)

    for m in 1 : M

        scene, structure = get_scene_and_structure(dset, m)
        roadway = get_roadway(scene)

        for vehicle_index in structure.active_vehicles

            lnP_tilde = 0.0
            fill!(lnP_tilde_denom_arr_s, 0.0)
            fill!(lnP_tilde_denom_arr_t, 0.0)
            fill!(lnP_tilde_denom_arr_v, 0.0)
            fill!(lnP_tilde_denom_arr_ϕ, 0.0)

            lo_s, hi_s = get_variable_bounds(IND_S, scene, structure, vehicle_index)
            lo_t, hi_t = get_variable_bounds(IND_T, scene, structure, vehicle_index)
            lo_v, hi_v = get_variable_bounds(IND_V, scene, structure, vehicle_index)
            lo_ϕ, hi_ϕ = get_variable_bounds(IND_ϕ, scene, structure, vehicle_index)

            for (factor_index, vehicle_indeces) in structure.factor_assignments

                ϕ = dset.factors[factor_index]
                pull_vehicle_data!(ϕ, scene, vehicle_indeces)
                extract_features!(ϕ, roadway)

                fd = factor_dot(ϕ)
                lnP_tilde += fd

                if in(vehicle_index, vehicle_indeces) # this factor affects the vehicle

                    base_index = VEHICLE_FEATURE_SIZE*(findfirst(vehicle_indeces, vehicle_index)-1)
                    @assert(base_index ≥ 0)

                    _update_lnP_tilde_denom_arr!(lnP_tilde_denom_arr_s, ϕ, base_index+IND_S, lo_s, hi_s, n_samples_monte_carlo_integration, roadway)
                    _update_lnP_tilde_denom_arr!(lnP_tilde_denom_arr_t, ϕ, base_index+IND_T, lo_t, hi_t, n_samples_monte_carlo_integration, roadway)
                    _update_lnP_tilde_denom_arr!(lnP_tilde_denom_arr_v, ϕ, base_index+IND_V, lo_v, hi_v, n_samples_monte_carlo_integration, roadway)
                    _update_lnP_tilde_denom_arr!(lnP_tilde_denom_arr_ϕ, ϕ, base_index+IND_ϕ, lo_ϕ, hi_ϕ, n_samples_monte_carlo_integration, roadway)

                else # this factor does not affect the vehicle
                    for k in 1 : n_samples_monte_carlo_integration
                        lnP_tilde_denom_arr_s[k] += fd
                        lnP_tilde_denom_arr_t[k] += fd
                        lnP_tilde_denom_arr_v[k] += fd
                        lnP_tilde_denom_arr_ϕ[k] += fd
                    end
                end
            end

            retval += lnP_tilde - _get_logP_tilde_denom(lnP_tilde_denom_arr_s, hi_s - lo_s, n_samples_monte_carlo_integration)
            retval += lnP_tilde - _get_logP_tilde_denom(lnP_tilde_denom_arr_t, hi_t - lo_t, n_samples_monte_carlo_integration)
            retval += lnP_tilde - _get_logP_tilde_denom(lnP_tilde_denom_arr_v, hi_v - lo_v, n_samples_monte_carlo_integration)
            retval += lnP_tilde - _get_logP_tilde_denom(lnP_tilde_denom_arr_ϕ, hi_ϕ - lo_ϕ, n_samples_monte_carlo_integration)
        end
    end
    retval /= M

    retval
end
function calc_weight_gradient(factor_index::Int, feature_index::Int, dset::SceneDataset, n_samples_monte_carlo_integration::Int, regularization::Float64, rng::AbstractRNG=Base.GLOBAL_RNG)

    retval = 0.0
    M = length(dset.scenes)
    ϕ = dset.factors[factor_index]
    target_feature = ϕ.features[feature_index]

    varscope = scope(target_feature)
    for j in 1 : length(varscope)
        if varscope[j] # in scope
            for m in 1 : M

                scene, structure = get_scene_and_structure(dset, m)
                roadway = get_roadway(scene)

                for (subfactor_index, vehicle_indeces) in structure.factor_assignments
                    if subfactor_index == factor_index # is the correct factor

                        index_in_vehicle = rem(j-1, VEHICLE_FEATURE_SIZE)+1
                        index_in_vehicle_indeces = div(j-1, VEHICLE_FEATURE_SIZE)+1
                        vehicle_index = vehicle_indeces[index_in_vehicle_indeces]

                        if in(vehicle_index, structure.active_vehicles)

                            pull_vehicle_data!(ϕ, scene, vehicle_indeces)

                            lo, hi = get_variable_bounds(index_in_vehicle, scene, structure, vehicle_index)
                            volume = hi - lo

                            # first component
                            if isa(target_feature, RoadFeature)
                                retval += evaluate(target_feature::RoadFeature, ϕ.x)
                                @assert(!isnan(retval))
                                @assert(!isinf(retval))
                            elseif isa(target_feature, RoadFeatureCont)
                                retval += evaluate(target_feature::RoadFeatureCont, ϕ.x)
                                @assert(!isnan(retval))
                                @assert(!isinf(retval))
                            elseif isa(target_feature, FollowFeature)
                                retval += evaluate(target_feature::FollowFeature, ϕ.x)
                                @assert(!isnan(retval))
                                @assert(!isinf(retval))
                            elseif isa(target_feature, FollowFeatureCont)
                                retval += evaluate(target_feature::FollowFeatureCont, ϕ.x, ϕ.vehicles)
                                @assert(!isnan(retval))
                                @assert(!isinf(retval))
                            else
                                retval += evaluate(target_feature::NeighborFeature, ϕ.x, ϕ.vehicles, ϕ.orig_states, roadway)
                                @assert(!isnan(retval))
                                @assert(!isinf(retval))
                            end

                            # compute expectation term
                            #   evaluated using Importance Sampling
                            E_numerator = 0.0
                            E_denominator = 0.0
                            for r in 1 : n_samples_monte_carlo_integration

                                ϕ.x[j] = volume*rand(rng) + lo
                                extract_features!(ϕ, roadway)

                                if isa(target_feature, RoadFeature)
                                    f = evaluate(target_feature::RoadFeature, ϕ.x)
                                elseif isa(target_feature, RoadFeatureCont)
                                    f = evaluate(target_feature::RoadFeatureCont, ϕ.x)
                                elseif isa(target_feature, FollowFeature)
                                    f = evaluate(target_feature::FollowFeature, ϕ.x)
                                elseif isa(target_feature, FollowFeatureCont)
                                    f = evaluate(target_feature::FollowFeatureCont, ϕ.x, ϕ.vehicles)
                                else
                                    f = evaluate(target_feature::NeighborFeature, ϕ.x, ϕ.vehicles, ϕ.orig_states, roadway)
                                end

                                p_true = eval_factor(ϕ)
                                p_importance = 1.0/volume
                                W = p_true / p_importance

                                if isnan(W) || isinf(W)
                                    println("feature 3: ", ϕ.features[3])
                                    println("x: ",         ϕ.x)
                                    println("weights: ", ϕ.weights)
                                    println("values:  ", ϕ.feature_values)
                                    println("factor dot: ", factor_dot(ϕ))
                                    println("p_true: ", p_true)
                                    println("p_imp:  ", p_importance)
                                    println("W:      ", W)
                                end

                                @assert(!isnan(p_true))
                                @assert(!isnan(p_importance))
                                @assert(!isnan(W))

                                E_numerator += f*W
                                E_denominator += W
                            end
                            # @assert(!isnan(E_numerator))
                            # @assert(!isnan(E_denominator))
                            # @assert(!isapprox(E_denominator, 0.0))

                            E = E_numerator / E_denominator

                            if isnan(E)
                                println("E:             ", E)
                                println("E_numerator:   ", E_numerator)
                                println("E_denominator: ", E_denominator)
                                sleep(0.1)
                                E = 0.0
                            end

                            @assert(!isnan(E))
                            @assert(!isinf(E))

                            retval -= E
                            @assert(!isnan(retval))
                            @assert(!isinf(retval))
                        end
                    end
                end
            end
        end
    end

    retval /= M

    # add regularization
    # retval -= 2*regularization*ϕ.weights[feature_index] # TEMP_REMOVE

    retval
end

type StochasticGradientAscentParams

    batch_size::Int
    batch_rng::AbstractRNG
    mc_pseudolikelihood_rng::AbstractRNG

    n_epochs::Int
    n_samples_monte_carlo_integration::Int
    n_samples_monte_carlo_pseudolikelihood::Int
    max_n_batches::Int

    learning_rate::Float64
    learning_rate_decay::Float64
    momentum_param::Float64      # ∈ [0,1), v ← γv + α∇ₜ(batch), γ typically starts at 0.5 and then is set to 0.9 later
    regularization::Float64
    factor_weight_min::Float64
    factor_weight_max::Float64

    save_every::Int            # [batches], set to -1 to never save
    save_dir::AbstractString   # directory in which to store checkpointed models
    same_name::AbstractString  # name for saved models -> ex 'model' -> model_001.jld

    print_to_stdout::Bool      # whether to print to STDOUT
    print_to_stdout_every::Int # [batches]
    print_to_file::AbstractString
    print_to_file_every::Int   # [batches]

    fout::Union{Void, IO}      # TODO: delete me

    function StochasticGradientAscentParams()
        retval = new()

        retval.batch_size = 100
        retval.batch_rng = Base.GLOBAL_RNG
        retval.mc_pseudolikelihood_rng = Base.GLOBAL_RNG

        retval.n_epochs = 5
        retval.n_samples_monte_carlo_integration = 10
        retval.n_samples_monte_carlo_pseudolikelihood = 10
        retval.max_n_batches = typemax(Int)

        retval.learning_rate = 1.0
        retval.learning_rate_decay = 0.97
        retval.momentum_param = 0.5
        retval.regularization = 0.001
        retval.factor_weight_min = -8.0
        retval.factor_weight_max =  0.0

        retval.save_every = 5
        retval.save_dir = SAVE_DIR
        retval.same_name = "full_scene"

        retval.print_to_stdout = false
        retval.print_to_stdout_every = 5
        retval.print_to_file = ""
        retval.print_to_file_every = 5

        retval
    end
end
function pprint(params::StochasticGradientAscentParams, tup::Tuple)
    if params.print_to_stdout
        print(tup...)
    end
    if !isempty(params.print_to_file)
        open(params.print_to_file, "a") do fout
            print(fout, tup...)
        end
    end
    nothing
end
function pprintln(params::StochasticGradientAscentParams, tup::Tuple)
    if params.print_to_stdout
        println(tup...)
    end
    if !isempty(params.print_to_file)
        open(params.print_to_file, "a") do fout
            println(fout, tup...)
        end
    end
    nothing
end
function get_param_dict(params::StochasticGradientAscentParams)
    retval = Dict{Symbol, Any}()

    retval[:batch_size] = params.batch_size
    retval[:n_epochs] = params.n_epochs
    retval[:n_samples_monte_carlo_integration] = params.n_samples_monte_carlo_integration
    retval[:n_samples_monte_carlo_pseudolikelihood] = params.n_samples_monte_carlo_pseudolikelihood
    retval[:learning_rate] = params.learning_rate
    retval[:learning_rate_decay] = params.learning_rate_decay
    retval[:momentum_param] = params.momentum_param
    retval[:regularization] = params.regularization
    retval[:factor_weight_min] = params.factor_weight_min
    retval[:factor_weight_max] = params.factor_weight_max
    retval[:save_every] = params.save_every
    retval[:save_dir] = params.save_dir
    retval[:same_name] = params.same_name
    retval[:print_to_stdout] = params.print_to_stdout
    retval[:print_to_stdout_every] = params.print_to_stdout_every
    retval[:print_to_file_every] = params.print_to_file_every

    retval
end
function set_params!(params::StochasticGradientAscentParams, dict::Dict{Symbol, Any})
    fields = fieldnames(params)
    for (sym, v) in dict
        if in(sym, fields)
            setfield!(params, sym, v)
        end
    end
    params
end

function sample_batch!(batch::SceneDataset, dset::SceneDataset, rng::AbstractRNG)
    data_range = 1:length(dset)
    for i in 1 : length(batch)
        index = rand(rng, data_range)
        batch.scenes[i] = dset.scenes[index]
        batch.sources[i] = dset.sources[index]
        batch.structures[i] = dset.structures[index]
    end
    batch
end
function stochastic_gradient_ascent!(dset::SceneDataset, params::StochasticGradientAscentParams)

    start_time = now()

    nsamples = length(dset)
    nbatches_per_epoch = div(nsamples, params.batch_size)
    if nbatches_per_epoch < 1
        warn("insufficient number of batches $nbatches_per_epoch for dataset with $nsamples samples!")
        return
    end

    if !isempty(params.print_to_file)
        # create / clear the file
        fout = open(params.print_to_file, "w")
        close(fout)
    end

    pprintln(params, ("nsamples:           ", nsamples))
    pprintln(params, ("batch_size:         ", params.batch_size))
    pprintln(params, ("nbatches per epoch: ", nbatches_per_epoch))
    pprintln(params, ("nepochs:            ", params.n_epochs))
    pprintln(params, ("niterations:        ", params.n_epochs*nbatches_per_epoch))
    pprintln(params, (@sprintf("%15s %10s %20s %20s %20s %20s", "iteration", "epoch", "score", "Δscore", "time", "w₁"),))
    pprintln(params, ("-"^(15+10+4*20+5),))
    sleep(0.1)

    batch_count = 0
    epoch_count = 0
    prev_score = -Inf
    score = calc_pseudolikelihood(dset, n_samples_monte_carlo_integration=params.n_samples_monte_carlo_pseudolikelihood)
    Δt = (now()-start_time).value / (1000* 60) # to minutes
    pprintln(params, (@sprintf("%15d %10d %20.5f %20.5f %20.2fm %20.2f", batch_count, epoch_count, score, score-prev_score, Δt, dset.factors[1].weights[1]),))
    sleep(0.1)

    # allocate batch
    batch = SceneDataset(Array(Scene, params.batch_size),
                         Array(SceneSource, params.batch_size),
                         Array(SceneStructure, params.batch_size),
                         dset.factors)

    # momentum
    α = params.learning_rate
    γ = params.momentum_param

    grad_velocitities = Dict{Symbol, Vector{Float64}}()
    for ϕ in batch.factors
        grad_velocitities[ϕ.name] = zeros(Float64, length(ϕ.weights))
    end

    # run gradient ascent
    while epoch_count < params.n_epochs && batch_count < params.max_n_batches

        # sample a new batch
        sample_batch!(batch, dset, params.batch_rng)

        # apply the momentum speed update
        for factor_index in 1 : length(batch.factors)
            ϕ = batch.factors[factor_index]
            sym = ϕ.name
            grad_vel_arr = grad_velocitities[sym]
            for feature_index in 1 : length(grad_velocitities[sym])
                grad = calc_weight_gradient(factor_index, feature_index, batch, params.n_samples_monte_carlo_integration, params.regularization, params.mc_pseudolikelihood_rng)
                @assert(!isnan(grad))
                @assert(!isinf(grad))
                grad_vel = grad_vel_arr[feature_index]
                @assert(!isnan(grad_vel))
                @assert(!isinf(grad_vel))
                grad_velocitities[sym][feature_index] = γ*grad_vel + α*grad
            end
        end

        # apply gradient update
        for ϕ in batch.factors
            grads = grad_velocitities[ϕ.name]

            for i in 1 : length(ϕ.weights)
                grad = clamp(grads[i], -0.1, 0.1)
                ϕ.weights[i] = clamp(ϕ.weights[i] + grad, params.factor_weight_min, params.factor_weight_max)
                @assert(!isnan(ϕ.weights[i]))
            end
        end
        α *= params.learning_rate_decay

        # update count
        batch_count += 1
        epoch_count = div(batch_count-1, nbatches_per_epoch)
        print_to_stdout = mod(batch_count, params.print_to_stdout_every) == 0 && params.print_to_stdout
        print_to_fout = mod(batch_count, params.print_to_file_every) == 0 && !isempty(params.print_to_file)
        if print_to_stdout || print_to_fout
            prev_score = score
            score = calc_pseudolikelihood(dset, n_samples_monte_carlo_integration=params.n_samples_monte_carlo_pseudolikelihood)
        end

        Δt = (now()-start_time).value / (1000* 60) # to minutes

        if print_to_stdout
            @printf("%15d %10d %20.5f %20.5f %20.2fm %20.2f\n", batch_count, epoch_count, score, score-prev_score, Δt, dset.factors[1].weights[1])
            sleep(0.1)
        end
        if print_to_fout
            open(params.print_to_file, "a") do fout
                @printf(fout, "%15d %10d %20.5f %20.5f %20.2fm %20.5f\n", batch_count, epoch_count, score, score-prev_score, Δt, dset.factors[1].weights[1])
            end
        end

        save_model = params.save_every ≥ 0 && (mod(batch_count, params.save_every) == 0 || epoch_count ≥ params.n_epochs)
        if save_model
            save_count = div(batch_count+1, params.save_every)+1
            filename = @sprintf("%s_%04d.jld", params.same_name, save_count)
            filepath = joinpath(params.save_dir, filename)
            pprint(params, ("saving to file: ", filename, "   ")); tic()
            JLD.save(filepath, "weights", get_weight_dict(dset.factors), "time", string(now()), "params", get_param_dict(params))
            pprintln(params, (toq(),))
            sleep(0.1)
        end
    end

    dset
end

immutable StochasticGradientAscentEntry
    iteration::Int
    epoch::Int
    score::Float64
    Δscore::Float64
    time::Float64
    w::Float64
end
immutable StochasticGradientAscentResult
    nsamples::Int
    batch_size::Int
    nbatches_per_epoch::Int
    nepochs::Int
    niterations::Int
    entries::Vector{StochasticGradientAscentEntry}
    save_files::Vector{AbstractString}
end
function read_logfile(io::IO)

    #=
    nsamples:           1000
    batch_size:         50
    nbatches per epoch: 20
    nepochs:            2
    niterations:        40
          iteration      epoch                score               Δscore                 time                   w₂
    --------------------------------------------------------------------------------------------------------------
                  0          0           -174.82021                  Inf                 0.90m             -5.00000
                  5          1           -173.24016              1.58005                 2.72m             -4.96943
    saving to file: full_scene_0002.jld   0.00472586
                 10          1           -172.65669              0.58347                 5.17m             -4.93655
    saving to file: full_scene_0003.jld   0.004537875
                 15          1           -172.69833             -0.04163                 6.97m             -4.90828
    saving to file: full_scene_0004.jld   0.004710128
                 20          1           -172.63179              0.06654                 8.80m             -4.88414
    saving to file: full_scene_0005.jld   0.030870091
    saving to file: full_scene_0005.jld   0.07673332
    =#


    nsamples = -1
    batch_size = -1
    nbatches_per_epoch = -1
    nepochs = -1
    niterations = -1
    entries = StochasticGradientAscentEntry[]
    save_files = AbstractString[]

    for line in readlines(io)
        if startswith(line, "nsamples:")
            nsamples = parse(Int, match(r"\d+", line).match)
        elseif startswith(line, "batch_size: ")
            batch_size = parse(Int, match(r"\d+", line).match)
        elseif startswith(line, "nbatches per epoch: ")
            nbatches_per_epoch = parse(Int, match(r"\d+", line).match)
        elseif startswith(line, "nepochs: ")
            nepochs = parse(Int, match(r"\d+", line).match)
        elseif startswith(line, "niterations: ")
            niterations = parse(Int, match(r"\d+", line).match)
        elseif startswith(line, "saving to file: ")
            push!(save_files, match(r"\S+.jld", line).match)
        elseif ismatch(r"\d+\s+\d+\s+[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?\s", line)
            matches = collect(eachmatch(FLOATING_POINT_REGEX, line))
            iteration = parse(Int, matches[1].match)
            epoch = parse(Int, matches[2].match)
            score = parse(Float64, matches[3].match)
            Δscore = parse(Float64, matches[4].match)
            time = parse(Float64, matches[5].match)
            w = parse(Float64, matches[6].match)
            push!(entries, StochasticGradientAscentEntry(iteration, epoch, score, Δscore, time, w))
        end
    end

    StochasticGradientAscentResult(nsamples, batch_size, nbatches_per_epoch, nepochs, niterations, entries, save_files)
end