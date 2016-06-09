include("piecewise_uniform.jl")

type UnivariateSceneGenerator <: OldSceneGenerator
    d_front :: ContinuousUnivariateDistribution
    v       :: ContinuousUnivariateDistribution
    d_cl    :: ContinuousUnivariateDistribution
    yaw     :: ContinuousUnivariateDistribution
end

function prob_gt(y::Float64, edges::Vector{Float64}, pdfs::Vector{Float64})
    # compute p(x > y) = ∫ pdf(x) dx from x = y to ∞

    a = y
    b = edges[end]

    if isapprox(a,b)
        return 0.0
    end
    @assert(b > a)

    retval = 0.0
    for i = 1 : length(edges)-1
        lo, hi = edges[i], edges[i+1]
        if a ≤ lo ≤ hi ≤ b # completely contained
            retval += pdfs[i] * (hi-lo)
        elseif lo ≤ a ≤ b ≤ hi # completely surrounds
            retval += pdfs[i] * (b-a)
            return retval
        elseif lo ≤ a ≤ hi # encapsulates a
            retval += pdfs[i] * (hi-a)
        elseif lo ≤ b ≤ hi # encapsulated b
            retval += pdfs[i] * (b-lo)
            return retval
        end
    end

    retval
end
function prob_gen(a::Float64, b::Float64, edges::Vector{Float64}, pdfs::Vector{Float64})
    # compute p(x) = ∫ pdf(u)/u du from u = a to b

    #=
    generative model:
        u ∼ f
        x ∼ U(0,u)

    p(x) = ∫f(u)p(x|u)du
         = ∫f(u)/u du from x to ∞

    p(x < y) = ∫∫f(u)/u du dx from x to ∞, from 0 to y
    =#

    if isapprox(a,b)
        return 0.0
    end
    @assert(b > a)

    retval = 0.0
    for i = 1 : length(edges)-1
        lo, hi = edges[i], edges[i+1]
        if a ≤ lo ≤ hi ≤ b # completely contained
            retval += pdfs[i] * log(hi/lo)
        elseif lo ≤ a ≤ b ≤ hi # completely surrounds
            retval += pdfs[i] * log(b/a)
            return retval
        elseif lo ≤ a ≤ hi # encapsulates a
            retval += pdfs[i] * log(hi/a)
        elseif lo ≤ b ≤ hi # encapsulated b
            retval += pdfs[i] * log(b/lo)
            return retval
        end
    end

    retval
end
function prob_uniform_sample_less_than(
    y       :: Float64,
    edges   :: Vector{Float64},
    pdfs    :: Vector{Float64};
    n_steps :: Int = 99
    )

    a = edges[1]
    c = edges[end]

    if y ≥ c
        return 1.0
    end
    b = y

    retval = prob_gen(b,c, edges, pdfs) # + prob_gen(a,c, edges, pdfs) NOTE(tim): this last part is 0

    if mod(n_steps,2) == 0
        n_steps +- 1
    end

    h = (b-a)/n_steps
    odd = true
    x = a
    for i = 1 : n_steps-1
        x += h
        retval += prob_gen(b,c,edges, pdfs) * (2+2*odd)
        odd = !odd
    end
    retval *= h/3

    return retval
end
prob_gt(d::PiecewiseUniform, y::Float64) = prob_gt(y, d.disc.binedges,d.pdfs)
prob_gen(d::PiecewiseUniform, a::Float64, b::Float64) = prob_gen(a, b, d.disc.binedges,d.pdfs)
prob_uniform_sample_less_than(d::PiecewiseUniform, y::Float64; n_steps :: Int = 99) =
    prob_uniform_sample_less_than(y, d.disc.binedges, d.pdfs, n_steps=n_steps)

function generate_scene(
    sg  :: UnivariateSceneGenerator,
    sti :: SceneTranslationInfo
    )

    road = sti.road
    nlanes = road.nlanes
    centers = lanecenters(road)
    roadway = NGSIM.ROADWAY_DICT[sti.roadway_name]

    vehicles = OldVehicle[]

    for lane in 1 : nlanes

        lti = sti.lanes[lane]
        curve = roadway.centerlines[lti.laneid]
        s_lo = curve[lti.curveind_lo].s
        s_hi = curve[lti.curveind_hi].s

        lane_center = centers[lane]
        offset_y = rand()*rand(sg.d_front) + s_lo
        v        = rand(sg.v)
        y        = offset_y
        push!(vehicles, OldVehicle(lane_center+rand(sg.d_cl),y,0.0,v))

        # now generate vehicles out forward
        done = false
        while !done

            d_front = rand(sg.d_front)

            if y + d_front+VEHICLE_MEAN_LENGTH < s_hi
                y += d_front+VEHICLE_MEAN_LENGTH
                push!(vehicles, OldVehicle(lane_center+rand(sg.d_cl),y,rand(sg.yaw),rand(sg.v)))
            else
                done = true
            end
        end
    end

    scene = RoadScene(road, vehicles)
    calc_vehicle_extra_info!(scene)
    scene
end
function loglikelihood(
    sg    :: UnivariateSceneGenerator,
    scene :: RoadScene
    )

    vehicles = scene.vehicles
    if isempty(scene.info)
        calc_vehicle_extra_info!(scene)
    end
    sceneinfo = scene.info

    score = 0.0

    for (i,veh) in enumerate(vehicles)
        info = sceneinfo[i]

        score += logpdf(sg.d_cl, clamp(info.d_cl, -5.75, 5.75))
        score += logpdf(sg.yaw,  clamp(veh.ϕ, -0.1, 0.1))
        score += logpdf(sg.v,    veh.v)

        if isinf(logpdf(sg.d_cl, clamp(info.d_cl, -5.75, 5.75)))
            println("d_cl: ", clamp(info.d_cl, -5.75, 5.75))
        end
        if isinf(logpdf(sg.yaw,  clamp(veh.ϕ, -0.1, 0.1)))
            println("yaw: ", clamp(veh.ϕ, -0.1, 0.1))
        end
        if isinf(logpdf(sg.v,    veh.v))
            println("v: ", veh.v)
        end

        if info.carind_fore != 0
            d_front = max(calc_d_front(scene, i), 0.0)
            score += logpdf(sg.d_front, d_front)
        else
            # incorporate score for offset_end
            # P(d_front > offset_end_observed)
            offset_end = scene.road.horizon_fore - veh.y - VEHICLE_MEAN_LENGTH
            if offset_end > 0.0
                score += log(prob_gt(sg.d_front, offset_end))
                if isinf(log(prob_gt(sg.d_front, offset_end)))
                    println("HERE 2!: ", offset_end)
                end
            end
        end

        if info.carind_rear == 0
            # incorporate score for offset_y
            # P(offset_y)
            offset_y = veh.y + scene.road.horizon_rear
            @assert(offset_y > 0.0)
            if isinf(log(prob_gen(sg.d_front, offset_y, Inf)))
                println("HERE!: ", offset_y)
            end
            score += log(prob_gen(sg.d_front, offset_y, Inf))
        end
    end

    @assert(!isinf(score))

    score
end
function train_univariate_piecewise_uniform(
    x::Vector{Float64},
    nbins::Int,
    bounds::Tuple{Float64, Float64},
    bin_spikes::Vector{Float64},
    spike_width::Float64 = 0.01;
    laplace_smoothing_counts::Int = 0
    )

    @assert(nbins > 0)

    m = length(x)

    bound_lo = bounds[1]
    bound_hi = bounds[2]

    edges = linspace(bound_lo, bound_hi, nbins+1)
    # now we need to intersect forced bins onto this
    for spike in bin_spikes
        lo = clamp(spike - spike_width/2, bound_lo, bound_hi)
        hi = clamp(spike + spike_width/2, bound_lo, bound_hi)

        @assert(!isapprox(lo,hi))

        for i = 1 : length(edges)-1
            edge_lo, edge_hi = edges[i], edges[i+1]
            if hi < edge_lo
                continue
            elseif edge_lo < lo && hi < edge_hi
                # spike is contained with an existing bin
                #  ->   just split within
                insert!(edges, i+1, lo)
                insert!(edges, i+2, hi)
                break
            elseif isapprox(edge_lo, lo) && hi < edge_hi
                insert!(edges, i+1, hi)
                break
            elseif isapprox(edge_lo, lo) && hi == edge_hi
                # do nothing
                break
            elseif edge_lo < lo && isapprox(hi, edge_hi)
                insert!(edges, i+1, lo)
                break
            elseif lo < edge_lo && hi > edge_lo
                # spike straddles an edge
                # -> remove existing edge and then split bit
                edges[i] = lo
                insert!(edges, i+1, hi)
                break
            end
        end
    end

    fit_mle(PiecewiseUniform, x, edges, laplace_smoothing_counts=laplace_smoothing_counts)
end
function train_univariatescenegenerator_piecewise_uniform(scenes::Vector{RoadScene},
    nbins_v        :: Int,
    nbins_d_front  :: Int,
    nbins_d_cl     :: Int,
    nbins_yaw      :: Int
    )

    if (BOUNDS_V_HI- BOUNDS_V_LO[1])/nbins_v < 0.1 ||
       (BOUNDS_ΔD_HI - BOUNDS_ΔD_LO)/nbins_d_front < 1.0 ||
       (BOUNDS_T_HI - BOUNDS_T_LO)/nbins_d_cl < 0.1 ||
       (BOUNDS_ϕ_HI - BOUNDS_ϕ_LO)/nbins_yaw < deg2rad(1.0)

        nbins_v = nbins_d_front = nbins_d_cl = nbins_yaw = 2
    end

    (arr_v, arr_d_front, arr_d_cl, arr_yaw) = get_base_arrays(scenes)

    # distr_v       = train_univariate_piecewise_uniform(arr_v, nbins_v, BOUNDS_V, [0.0:5.0:75.0], laplace_smoothing_counts=1)
    distr_v       = fit_mle(PiecewiseUniform, arr_v,       nbins_v,       bounds = (BOUNDS_V_LO,       BOUNDS_V_HI),       laplace_smoothing_counts=1)
    distr_d_front = fit_mle(PiecewiseUniform, arr_d_front, nbins_d_front, bounds = (BOUNDS_ΔD_LO, BOUNDS_ΔD_HI), laplace_smoothing_counts=1)
    distr_d_cl    = fit_mle(PiecewiseUniform, arr_d_cl,    nbins_d_cl,    bounds = (BOUNDS_T_LO,       BOUNDS_T_HI),       laplace_smoothing_counts=1)
    distr_yaw     = fit_mle(PiecewiseUniform, arr_yaw,     nbins_yaw,     bounds = (BOUNDS_ϕ_LO,       BOUNDS_ϕ_HI),       laplace_smoothing_counts=1)

    UnivariateSceneGenerator(distr_d_front, distr_v, distr_d_cl, distr_yaw)
end

function train{T<:UnivariateSceneGenerator}(res::ModelOptimizationResults{T}, scenes::Vector{RoadScene})
    p = res.params
    nbins_v         = p[:nbins_v]::Int
    nbins_d_front = p[:nbins_d_front]::Int
    nbins_d_cl      = p[:nbins_d_cl]::Int
    nbins_yaw     = p[:nbins_yaw]::Int
    train_univariatescenegenerator_piecewise_uniform(scenes, nbins_v, nbins_d_front, nbins_d_cl, nbins_yaw)
end

function cross_validate_scenegenerator(
    ::Type{UnivariateSceneGenerator},
    scenes::Vector{RoadScene};
    CV_nfolds :: Int = 10, # k in k-fold cross validations
    CV_rounds :: Int = 10, # number of CV evals to do
    nbins_v         :: Int = 20,
    nbins_d_front :: Int = 20,
    nbins_d_cl      :: Int = 20,
    nbins_yaw     :: Int = 20,
    )

    nscenes = length(scenes)

    (arr_v, arr_d_front, arr_d_cl, arr_yaw) = get_base_arrays(scenes)

    # arr_v += (0.02*rand(length(arr_v))-0.01)

    # bounds_v         = extrema(arr_v)
    # bounds_d_front = extrema(arr_d_front)
    # bounds_d_cl      = extrema(arr_d_cl)
    # bounds_yaw     = extrema(arr_yaw)

    f_train(inds::Vector{Int}) = train_univariatescenegenerator_piecewise_uniform(scenes[inds],
        nbins_v, nbins_d_front, nbins_d_cl, nbins_yaw
        # bounds_v, bounds_d_front, bounds_d_cl, bounds_yaw
        )
    f_test(sg::SceneGenerator, inds::Vector{Int}) = begin
        score = 0.0
        for ind in inds
            score += loglikelihood(sg, scenes[ind])
        end
        score
    end

    scores = Array(Float64, CV_nfolds, CV_rounds)
    for round = 1 : CV_rounds
        for (fold, train_inds) in enumerate(Kfold(nscenes, CV_nfolds))
            test_inds = setdiff(1:nscenes, train_inds)
            model = f_train(train_inds)
            scores[fold,round] = f_test(model, test_inds)
        end
    end
    # for round = 1 : CV_rounds
    #     scores[:,round] = MLBase.cross_validate(f_train, f_test, nscenes, Kfold(nscenes, CV_nfolds))
    # end

    scores
end

function cross_validate(
    ::Type{UnivariateSceneGenerator},
    scenes::Vector{RoadScene},
    nbins::Dict{Symbol,Int};
    CV_nfolds :: Int = 10, # k in k-fold cross validations
    CV_rounds :: Int = 10  # number of CV evals to do
    )

    scores = cross_validate_scenegenerator(UnivariateSceneGenerator, scenes, CV_rounds = CV_rounds, CV_nfolds =CV_nfolds,
        nbins_v = nbins[:speed],
        nbins_d_front = nbins[:d_front],
        nbins_d_cl = nbins[:d_cl],
        nbins_yaw = nbins[:yaw]
        )
    ave_crossvalidated_likelihood(scores)
end
function cyclic_coordinate_ascent_parallel(
    T::Type{UnivariateSceneGenerator},
    scenes::Vector{RoadScene};
    CV_nfolds :: Int = 10,
    CV_rounds :: Int = 10
    )

    param_options = (
            [3,5,10,50,100,250,500,750,1000,2500,5000], # :speed
            [3,5,10,50,100,250,500,750,1000,2500,5000], # :d_cl
            [3,5,10,50,100,250,500,750,1000,2500,5000], # :yaw
            [3,5,10,50,100,250,500,750,1000,2500,5000], # :d_front
        )

    n_params = length(param_options)
    paraminds = ones(Int, n_params)
    params_tried = Set{Vector{Int}}()

    converged = false
    best_μ = -Inf
    best_σ = NaN
    iter = 0
    while !converged
        iter += 1
        println("iteration: ", iter)

        converged = true
        for coordinate = 1 : n_params
            println("\tcoordinate: ", coordinate)

            to_try = Vector{Int}[]
            ip_arr = Int[]
            for ip in 1 : length(param_options[coordinate])

                newparams = copy(paraminds)
                newparams[coordinate] = ip

                if !in(newparams, params_tried)
                    push!(params_tried, newparams)
                    push!(ip_arr, ip)
                    push!(to_try, newparams)
                end
            end

            if !isempty(to_try)

                tic()
                scores = pmap(to_try) do newparams
                    cross_validate(UnivariateSceneGenerator, scenes,
                       Dict(
                        :speed      => param_options[1][newparams[1]],
                        :d_cl       => param_options[2][newparams[2]],
                        :yaw        => param_options[3][newparams[3]],
                        :d_front  => param_options[4][newparams[4]],
                        ),CV_rounds = CV_rounds, CV_nfolds =CV_nfolds)
                end
                toc()

                run_ind = 0
                run_μ = -Inf
                run_σ = NaN
                for i = 1 : length(scores)
                    μ = scores[i][1]
                    σ = scores[i][2]
                    if μ > run_μ
                        run_ind, run_μ, run_σ = i, μ, σ
                    end
                end

                if run_μ > best_μ
                    best_μ = run_μ
                    best_σ = run_σ
                    paraminds[coordinate] = ip_arr[run_ind]
                    converged = false
                    println("\tfound better: ", coordinate, " -> ", param_options[coordinate][ip_arr[run_ind]])
                    println("\t\tscore: ", best_μ)
                end
            end
        end
    end

    println("done!")
    println("optimal params: ", )
    println("\tnbins_speed     = ", param_options[1][paraminds[1]])
    println("\tnbins_d_cl      = ", param_options[2][paraminds[2]])
    println("\tnbins_yaw       = ", param_options[3][paraminds[3]])
    println("\tnbins_d_front = ", param_options[4][paraminds[4]])
    println("score mean: ", best_μ)
    println("score std:  ", best_σ)
    println("num iter:   ", iter)

    opt_params = Dict{Symbol,Any}()
    opt_params[:nbins_v        ] = param_options[1][paraminds[1]]
    opt_params[:nbins_d_front] = param_options[2][paraminds[2]]
    opt_params[:nbins_d_cl     ] = param_options[3][paraminds[3]]
    opt_params[:nbins_yaw      ] = param_options[4][paraminds[4]]

    ModelOptimizationResults{T}(opt_params, best_μ, best_σ, best_μ, iter, CV_nfolds, CV_rounds)
end