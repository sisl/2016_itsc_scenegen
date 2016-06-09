function get_first_colliding_pair(rec::SceneRecord, pastframe::Int, mem::CPAMemory = CPAMemory())
    j = NGSIM._state_ind(pastframe)
    for a in 1 : rec.n_vehicles[j]

        vehA = rec.vehicles[rec.ids[a,j]]
        vehA.state = rec.states[a,j]

        if !isnan(vehA.state.v) # ignore missing vehicles
            set_to_positioned_oriented_bounding_box!(mem.vehA, vehA)

            for b in a + 1 : rec.n_vehicles[j]

                vehB = rec.vehicles[rec.ids[b,j]]
                vehB.state = rec.states[b,j]

                if !isnan(vehB.state.v) # ignore missing vehicles
                    set_to_positioned_oriented_bounding_box!(mem.vehB, vehB)

                    if is_colliding(mem.vehA, mem.vehB, mem.mink)
                        return (a, b)
                    end
                end
            end
        end
    end

    return (0,0) # none found
end
function get_min_dist_btwn_vehicles(rec::SceneRecord, pastframe::Int, mem::CPAMemory = CPAMemory())

    a, b = get_first_colliding_pair(rec, pastframe, mem)
    if a != 0
        return 0.0
    end

    Inf
end
function get_first_collision_type(rec::SceneRecord, mem::CPAMemory = CPAMemory())
    for i in 1 : record_length(rec)
        pastframe = 1 - i
        (a,b) = get_first_colliding_pair(rec, pastframe, mem)
        if a != 0 # collision occured
            return classify_collision(rec[a, pastframe], rec[b, pastframe])
        end
    end
    return -1 # no collision
end
function get_time_of_first_collision(rec::SceneRecord, mem::CPAMemory=CPAMemory())
    for i in 1 : record_length(rec)
        pastframe = 1 - i
        (a,b) = get_first_colliding_pair(rec, pastframe, mem)
        if a != 0 # collision occured
            return (i-1)*NGSIM_TIMESTEP
        end
    end
    return NaN # no collision
end

function get_target_function_surrogate(rec::SceneRecord, mem::CPAMemory = CPAMemory())
    min_dist = Inf
    for i in 1 : record_length(rec)
        min_dist = min(get_min_dist_btwn_vehicles(rec, 1-i, mem), min_dist)
        if min_dist == 0.0
            return 0.0
        end
    end
    min_dist
end
target_function(min_distance_btwn_vehicles::Float64) = convert(Float64, min_distance_btwn_vehicles == 0.0)
function target_function(rec::SceneRecord)
    # true if there is a collision
    min_distance_btwn_vehicles = get_target_function_surrogate(rec)
    target_function(min_distance_btwn_vehicles)
end

const COLLISION_TYPE_NULL = 0 # when it doesn't fit the other categories
const COLLISION_TYPE_LONGITUDINAL = 1 # frontal / rear-end collision
const COLLISION_TYPE_LATERAL = 2 # left / right side collision


function classify_collision(A::Vehicle, B::Vehicle)

    a_rel_b = inertial2body(A.state.posG, B.state.posG)
    b_rel_a = inertial2body(B.state.posG, A.state.posG)

    if -30 < rad2deg(b_rel_a.θ) < 30 &&
       -30 < rad2deg(a_rel_b.θ) < 30
       # both vehicles are relatively facing parallel

        θb = rad2deg(atan2(b_rel_a.y, b_rel_a.x))
        θa = rad2deg(atan2(a_rel_b.y, a_rel_b.x))

        if (abs(θb) < 30 || abs(θb) > 125) &&
           (abs(θa) < 30 || abs(θa) > 125)

            return COLLISION_TYPE_LONGITUDINAL
        else
            return COLLISION_TYPE_LATERAL
        end
    end

    return COLLISION_TYPE_NULL
end

############################

function _get_closest_lead_follow_pair(scene::Scene)

    smallest_dist = Inf
    retval = (0,0)

    for (a, vehA) in enumerate(scene)
        b = get_neighbor_index_fore(scene, a)
        if b != 0
            dist = get_headway_dist_between(vehA, scene[b])
            if !isnan(dist) && dist < smallest_dist
                smallest_dist = dist
                retval = (a, b)
            end
        end
    end

    retval
end

function _get_drivers_dict(driver::DriverModel, scene::Scene)
    drivers = Dict{Int, DriverModel}()
    for veh in scene
        driver_copy = copy(driver)
        driver_copy.v_des = max(veh.state.v, 5.0)
        reset_hidden_state!(driver_copy)
        drivers[veh.id] = driver_copy
    end
    drivers
end

############################

immutable SamplingParams2
    factors::Vector{Factor}
    transition::ContinuousMultivariateDistribution
    nsteps_burnin::Int # number of burn-in steps for metropolis hastings
    dset::SceneDataset
    scene::Scene
end
function Base.rand!(sampler::SamplingParams2)
    range = 1:length(sampler.dset)
    scene_index = rand(range)
    scene, structure = get_scene_and_structure(sampler.dset, scene_index)
    copy!(sampler.scene, scene)
    metropolis_hastings!(sampler.scene, structure, sampler.factors, sampler.transition, sampler.nsteps_burnin)
    scene_source = dset.sources[scene_index]
    (sampler.scene, structure, scene_source)
end

############################

type SimParams
    frames_per_tick::Int
    sim_horizon::Int
    rec::SceneRecord
    actions::Vector{DriveAction}

    function SimParams(;
        frames_per_tick::Int=3,
        sim_horizon::Int=frames_per_tick*10,
        )

        rec = SceneRecord(sim_horizon+1)
        actions = Array(DriveAction, 100)
        new(frames_per_tick, sim_horizon, rec, actions)
    end
end

type DirectPropSampling <: ExpectationSolver

    sampler::SamplingParams2
    driver::DriverModel
    simparams::SimParams

    n::Int      # sample count
    μ::Float64  # online mean
    M::Float64  # sum of squares for differences of the current mean (helper variable)

    DirectPropSampling(sampler::SamplingParams2, driver::DriverModel, simparams::SimParams) = new(sampler, driver, simparams, 0, 0.0, 0.0)
end
function Base.empty!(sol::DirectPropSampling)
    sol.n = 0
    sol.μ = 0.0
    sol.M = 0.0
    sol
end
Base.copy(sol::DirectPropSampling) = DirectPropSampling(sol.sampler, sol.driver, sol.simparams)
get_nsamples(sol::DirectPropSampling) = sol.n
Base.mean(sol::DirectPropSampling) = sol.μ
Base.var(sol::DirectPropSampling) = sol.M / (sol.n - 1)
Base.std(sol::DirectPropSampling) = sqrt(var(sol))

function simulate_into_record!(rec::SceneRecord, scene::Scene, drivers::Dict{Int, DriverModel};
    frames_per_tick::Int = 3,
    sim_horizon::Int = record_length(rec)-1,
    actions::Vector{DriveAction} = Array(DriveAction, length(drivers)),
    )

    NGSIM.update!(rec, scene)
    for h in 1 : frames_per_tick : sim_horizon
        get_actions!(actions, scene, drivers)
        for t in 1 : frames_per_tick
            tick!(scene, drivers, actions)
            NGSIM.update!(rec, scene)
        end
    end
end

function tick!(sol::DirectPropSampling)

    scene, structure, source = rand!(sol.sampler)

    # specify driver models
    drivers = _get_drivers_dict(sol.driver, scene)

    # run simulation
    simparams = sol.simparams
    empty!(simparams.rec)
    NGSIM.update!(simparams.rec, scene)
    for h in 1 : simparams.frames_per_tick : simparams.sim_horizon
        get_actions!(simparams.actions, scene, drivers)
        for t in 1 : simparams.frames_per_tick
            tick!(scene, simparams.actions)
            NGSIM.update!(simparams.rec, scene)
        end
    end

    # run the target function on the record
    f = target_function(simparams.rec)

    n, μ, M = sol.n, sol.μ, sol.M

    # online update of mean and the variance-helper variable M
    n += 1
    μ_next = μ + (f - μ)/n
    M += (f - μ)*(f - μ_next)

    sol.n, sol.μ, sol.M = n, μ_next, M
    sol
end

############################

type ImportancePropSampling <: ExpectationSolver

    sampler::SamplingParams2
    driver::DriverModel
    simparams::SimParams

    n::Int
    w_arr::Vector{Float64}
    f_arr::Vector{Float64}

    function ImportancePropSampling(sampler::SamplingParams2, driver::DriverModel, simparams::SimParams, max_nsamples::Int)

        retval = new()
        retval.sampler = sampler
        retval.driver = driver
        retval.simparams = simparams
        retval.n = 0
        retval.w_arr = Array(Float64, max_nsamples)
        retval.f_arr = Array(Float64, max_nsamples)
        retval
    end
end
function Base.empty!(sol::ImportancePropSampling)
    sol.n = 0
    sol
end
Base.copy(sol::ImportancePropSampling) = ImportancePropSampling(sol.sampler, sol.driver, sol.simparams, length(sol.w_arr))

get_nsamples(sol::ImportancePropSampling) = sol.n
function Base.mean(sol::ImportancePropSampling)
    num = 0.0
    den = 0.0
    for i in 1 : get_nsamples(sol)
        f, w = sol.f_arr[i], sol.w_arr[i]
        num += f*w
        # den += w # for self-normalizing
        den += 1.0 # for standard
    end
    # @printf("den: %10.6f  %10d\n", den, get_nsamples(sol))
    num/den
end
function Base.var(sol::ImportancePropSampling, μ::Float64 = mean(sol))
    tot = 0.0
    n = get_nsamples(sol)
    for i in 1 : n
        f, w = sol.f_arr[i], sol.w_arr[i]
        Δ = (f - μ)
        tot += w*w*Δ*Δ
    end
    tot / (n-1)
end
Base.std(sol::ImportancePropSampling) = sqrt(var(sol))

function effective_sample_size(sol::ImportancePropSampling)
    sum_w = 0.0
    sum_w2 = 0.0
    for i in 1 : get_nsamples(sol)
        w = sol.w_arr[i]
        sum_w += w
        sum_w2 += w*w
    end
    sum_w*sum_w/sum_w2
end
function effective_sample_size_variance(sol::ImportancePropSampling)
    sum_w2 = 0.0
    sum_w4 = 0.0
    for i in 1 : get_nsamples(sol)
        w = sol.w_arr[i]
        w² = w*w
        sum_w2 += w²
        sum_w4 += w²*w²
    end
    sum_w2*sum_w2/sum_w4
end
function effective_sample_size_f(sol::ImportancePropSampling)
    den = 0.0
    for i in 1 : get_nsamples(sol)
        f, w = sol.f_arr[i], sol.w_arr[i]
        den += abs(f)*w
    end

    sum_w_tilde2 = 0.0
    for i in 1 : get_nsamples(sol)
        f, w = sol.f_arr[i], sol.w_arr[i]
        w_tilde = f*w / den
        sum_w_tilde2 += w_tilde*w_tilde
    end

    1 / sum_w_tilde2
end

function tick!(sol::ImportancePropSampling)

    scene, structure, source = rand!(sol.sampler)

    bias_a_fore = 0.0   # ft/s2 # bias to decel
    bias_a_rear = 0.01  # ft/s2 # bias to accel

    # specify driver models
    drivers = _get_drivers_dict(sol.driver, scene)

    # run simulation
    log_likelihood_ratio = 0.0
    simparams = sol.simparams
    empty!(simparams.rec)
    NGSIM.update!(simparams.rec, scene)
    for h in 1 : simparams.frames_per_tick : simparams.sim_horizon

        # get the two closest lead-trail vehicle pair
        a, b = _get_closest_lead_follow_pair(scene)

        # get actions
        for (veh_index, veh) in enumerate(scene)

            model = drivers[veh.id]
            observe!(model, scene, veh.id)

            if veh_index == a # follower
                N_a = Normal(model.P.μ[1] + bias_a_rear, sqrt(model.P.Σ.mat[1,1]))
                N_ϕdes = Normal(model.P.μ[2], sqrt(model.P.Σ.mat[2,2]))
                action = DriveAction(rand(N_a), rand(N_ϕdes))

                log_likelihood_ratio -= (logpdf(N_a, action.a) + logpdf(N_ϕdes, action.ϕdes))
                log_likelihood_ratio += get_action_logl(model, action)
                simparams.actions[veh_index] = action
            elseif veh_index == b # leader
                N_a = Normal(model.P.μ[1] + bias_a_fore, sqrt(model.P.Σ.mat[1,1]))
                N_ϕdes = Normal(model.P.μ[2], sqrt(model.P.Σ.mat[2,2]))
                action = DriveAction(rand(N_a), rand(N_ϕdes))

                log_likelihood_ratio -= (logpdf(N_a, action.a) + logpdf(N_ϕdes, action.ϕdes))
                log_likelihood_ratio += get_action_logl(model, action)
                simparams.actions[veh_index] = action
            else # normal - no bias
                simparams.actions[veh_index] = rand(model)
            end
        end

        for t in 1 : simparams.frames_per_tick
            tick!(scene, simparams.actions)
            NGSIM.update!(simparams.rec, scene)
        end
    end
    # run the target function on the record
    f = target_function(simparams.rec)
    w = exp(log_likelihood_ratio)

    sol.n += 1
    sol.w_arr[sol.n] = w
    sol.f_arr[sol.n] = f

    sol
end