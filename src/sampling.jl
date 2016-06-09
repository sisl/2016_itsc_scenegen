const TARGET_THRESHOLD = 16.4 # [ft/s²]
function get_required_acceleration(rear::Vehicle, fore::Vehicle)

    # the const. acceleration required to avoid a collision assuming
    # everyone drives with constant frenet-x velocity
    # NOTE: this assumes that fore and rear are in the same lane

    @assert(fore.state.posF.laneid == rear.state.posF.laneid)

    t_CPA, d_CPA = closest_time_of_approach_and_distance(fore.state.posG, fore.state.v,
                                                         rear.state.posG, rear.state.v)

    acceptance_gap = (fore.width + rear.width)/2

    if t_CPA < 0.0 || d_CPA > acceptance_gap
        return 0.0
    end

    Δv = fore.state.v - rear.state.v

    if Δv > 0.0
        0.0
    else
        Δx = fore.state.posF.s - rear.state.posF.s - fore.length
        # @assert(Δx > 0.0)

        a_req = -Δv*Δv / (2Δx)

        min(a_req, 0.0)
    end
end
function get_target_function_surrogate(scene::Scene, structure::SceneStructure; factor_index_following::Int=2)

    # value we try to maximize in Cross Entropy
    # S, the "performance metric"

    min_target = Inf

    for (factor_index, vehicle_indeces) in structure.factor_assignments
        if factor_index == factor_index_following
            rear = scene.vehicles[vehicle_indeces[1]]
            fore = scene.vehicles[vehicle_indeces[2]]
            a_req = get_required_acceleration(rear, fore)
            min_target = min(min_target, a_req)
        end
    end

    -min_target
end
target_function(target_function_surrogate::Float64) = convert(Float64, target_function_surrogate > TARGET_THRESHOLD)
function target_function(scene::Scene, structure::SceneStructure; factor_index_following::Int=2)
    target = get_target_function_surrogate(scene, structure, factor_index_following=factor_index_following)
    target_function(target)
end

############################

immutable SamplingParams
    factors::Vector{Factor}
    transition::ContinuousMultivariateDistribution
    nsteps_burnin::Int # number of burn-in steps for metropolis hastings
    dset::SceneDataset
    scene::Scene
    weights::Vector{Float64}

    function SamplingParams(
        factors::Vector{Factor},
        transition::ContinuousMultivariateDistribution,
        nsteps_burnin::Int, # number of burn-in steps for metropolis hastings
        dset::SceneDataset,
        scene::Scene=Scene(),
        )

        n = length(dset.scenes)
        weights = fill(1.0/n, n)

        new(factors, transition, nsteps_burnin, dset, scene, weights)
    end
    function SamplingParams(
        factors::Vector{Factor},
        transition::ContinuousMultivariateDistribution,
        nsteps_burnin::Int, # number of burn-in steps for metropolis hastings
        dset::SceneDataset,
        weights::Vector{Float64},
        scene::Scene=Scene(),
        )

        @assert(length(dset.scenes) == length(weights))
        tot = sum(weights)
        @assert(tot > 0.0)
        weights ./= tot

        new(factors, transition, nsteps_burnin, dset, scene, weights)
    end
    function SamplingParams(
        sampler::SamplingParams,
        weights::Vector{Float64},
        )

        @assert(length(sampler.dset.scenes) == length(weights))
        tot = sum(weights)
        @assert(tot > 0.0)
        weights ./= tot

        new(sampler.factors, sampler.transition, sampler.nsteps_burnin, sampler.dset, Scene(), weights)
    end
end
function draw_from_categorical(weights::Vector{Float64})
    r = rand()
    i = 1
    t = weights[i]
    while r > t && i < length(weights)
        i += 1
        t += weights[i]
    end
    @assert(i ≤ length(weights))
    i
end
function Base.rand!(sampler::SamplingParams)

    i = draw_from_categorical(sampler.weights)
    scene, structure = get_scene_and_structure(sampler.dset, i)
    copy!(sampler.scene, scene)
    metropolis_hastings!(sampler.scene, structure, sampler.factors, sampler.transition, sampler.nsteps_burnin)
    structure
end

############################

abstract ExpectationSolver
function tick!(exp::ExpectationSolver)
    warn("tick! NOT IMPLEMENTED FOR $(typeof(exp))")
    exp
end
function solve!(sol::ExpectationSolver, nsamples::Int)

    for i in 1 : nsamples
        tick!(sol)
    end
    sol
end

############################

type DirectSampling <: ExpectationSolver

    sampler::SamplingParams

    n::Int      # sample count
    μ::Float64  # online mean
    M::Float64  # sum of squares for differences of the current mean (helper variable)

    DirectSampling(sampler::SamplingParams) = new(sampler, 0, 0.0, 0.0)
end
Base.copy(sol::DirectSampling) = DirectSampling(sol.sampler)
get_nsamples(sol::DirectSampling) = sol.n
Base.mean(sol::DirectSampling) = sol.μ
Base.var(sol::DirectSampling) = sol.M / (sol.n - 1)
Base.std(sol::DirectSampling) = sqrt(var(sol))

function tick!(sol::DirectSampling)

    structure = rand!(sol.sampler)
    f = target_function(sol.sampler.scene, structure)

    n, μ, M = sol.n, sol.μ, sol.M

    # online update of mean and the variance-helper variable M
    n += 1
    μ_next = μ + (f - μ)/n
    M += (f - μ)*(f - μ_next)

    sol.n, sol.μ, sol.M = n, μ_next, M
    sol
end

###########

type ImportanceSampling <: ExpectationSolver

    P::SamplingParams
    Q::SamplingParams

    n::Int
    w_arr::Vector{Float64}
    f_arr::Vector{Float64}

    function ImportanceSampling(P::SamplingParams, Q::SamplingParams, max_nsamples::Int)

        retval = new()
        retval.P = P
        retval.Q = Q
        retval.n = 0
        retval.w_arr = Array(Float64, max_nsamples)
        retval.f_arr = Array(Float64, max_nsamples)
        retval
    end
end
Base.copy(sol::ImportanceSampling) = ImportanceSampling(sol.P, sol.Q, length(sol.w_arr))
function Base.empty!(sol::ImportanceSampling)
    sol.n = 0
end

get_nsamples(sol::ImportanceSampling) = sol.n
function Base.mean(sol::ImportanceSampling)
    num = 0.0
    den = 0.0
    for i in 1 : get_nsamples(sol)
        f, w = sol.f_arr[i], sol.w_arr[i]
        num += f*w
        den += w # for self-normalizing
        # den += 1.0 # for standard
    end
    @printf("den: %10.6f  %10d\n", den, get_nsamples(sol))
    num/den
end
function Base.var(sol::ImportanceSampling, μ::Float64 = mean(sol))
    tot = 0.0
    n = get_nsamples(sol)
    for i in 1 : n
        f, w = sol.f_arr[i], sol.w_arr[i]
        Δ = (f - μ)
        tot += w*w*Δ*Δ
    end
    tot / (n-1)
end
Base.std(sol::ImportanceSampling) = sqrt(var(sol))

function effective_sample_size(sol::ImportanceSampling)
    sum_w = 0.0
    sum_w2 = 0.0
    for i in 1 : get_nsamples(sol)
        w = sol.w_arr[i]
        sum_w += w
        sum_w2 += w*w
    end
    sum_w*sum_w/sum_w2
end
function effective_sample_size_variance(sol::ImportanceSampling)
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
function effective_sample_size_f(sol::ImportanceSampling)
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

function tick!(sol::ImportanceSampling)

    i = draw_from_categorical(sol.Q.weights)

    scene, structure = get_scene_and_structure(sol.Q.dset, i)
    copy!(sol.Q.scene, scene)
    metropolis_hastings!(sol.Q.scene, structure, sol.Q.factors, sol.Q.transition, sol.Q.nsteps_burnin)
    f = target_function(sol.Q.scene, structure)

    probP = sol.P.weights[i] * evaluate(sol.P.factors, sol.Q.scene, structure)
    probQ = sol.Q.weights[i] * evaluate(sol.Q.factors, sol.Q.scene, structure)
    w = probP/probQ

    sol.n += 1
    sol.w_arr[sol.n] = w
    sol.f_arr[sol.n] = f

    sol
end