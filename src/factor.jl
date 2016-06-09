const VEHICLE_FEATURE_SIZE = 4 # {s, t, v, ϕ}
const IND_S = 1
const IND_T = 2
const IND_V = 3
const IND_ϕ = 4

const BOUNDS_V_LO =  -2.5 # [ft/s]
const BOUNDS_V_HI =  95.1 # [ft/s]
const BOUNDS_T_LO = -17.2 # [ft]
const BOUNDS_T_HI =  11.9 # [ft]
const BOUNDS_ϕ_LO = -0.46 # [rad]
const BOUNDS_ϕ_HI =  0.82 # [rad]
const BOUNDS_ΔD_LO = 0.0   # [ft]
const BOUNDS_ΔD_HI = 320.0 # [ft]
const BOUNDS_ΔV_LO = -51.6 # [ft/s]
const BOUNDS_ΔV_HI =  50.2 # [ft/s]

###############################
# Feature

abstract Feature
get_nvehicles(f::Feature) = div(length(scope(f)), VEHICLE_FEATURE_SIZE)
macro get_bitvector(var_expr::Expr, nvehicles::Int)

    vars = eval(var_expr)::Vector{Symbol}

    # construct the bitvector for the given symbols and number of vehicles
    # the bit vector is of length VEHICLE_FEATURE_SIZE*nvehicles
    # vars contains symbols of the form:
    #       {s, t, v, ϕ}{number}
    #       ex: [:v1, ϕ2] would be a factor that uses velocity of vehicle 1 and heading of vehicle 2
    #           and get_bitvector([:v1, ϕ2], 2) -> [F, F, T, F, F, F, F, T]

    nbits = nvehicles*VEHICLE_FEATURE_SIZE
    retval = falses(nbits)

    for var in vars
        chars = collect(string(var))
        bit_offset = 0
        if chars[1] == 's'
            bit_offset = IND_S
        elseif chars[1] == 't'
            bit_offset = IND_T
        elseif chars[1] == 'v'
            bit_offset = IND_V
        elseif chars[1] == 'ϕ'
            bit_offset = IND_ϕ
        else
            error("unrecognized first token in var $var")
        end

        veh_index = parse(Int, chars[2])
        @assert(1 ≤ veh_index ≤ nvehicles)

        retval[VEHICLE_FEATURE_SIZE*(veh_index-1) + bit_offset] = true
    end

    retval
end

immutable RoadFeature <: Feature
    t_bin::Int
    v_bin::Int
    ϕ_bin::Int
end
# const _RoadFeature_BINS_T = [-Inf; linspace(BOUNDS_T_LO, BOUNDS_T_HI, 8)[2:end-1]; Inf]
const _RoadFeature_BINS_V = [-Inf; linspace(BOUNDS_V_LO, BOUNDS_V_HI, 8)[2:end-1]; Inf]
# const _RoadFeature_BINS_ϕ = [-Inf; linspace(BOUNDS_ϕ_LO, BOUNDS_ϕ_HI, 8)[2:end-1]; Inf]
const _RoadFeature_BINS_T = [-Inf, -6.0, -2.0, -1.0, 0.0, 1.0, 2.0, Inf]
const _RoadFeature_BINS_ϕ = [-Inf, deg2rad(-4.0), deg2rad(-2.0), deg2rad(-1.0), deg2rad(1.0), deg2rad(2.0), deg2rad(4.0), Inf]
const _RoadFeatureScope = @get_bitvector([:t1, :v1, :ϕ1], 1)
scope(::RoadFeature) = _RoadFeatureScope
function evaluate(f::RoadFeature, x::Vector{Float64})

    # DEBUG
    # @printf("\t%3d  %3d  %3d\n", f.t_bin, f.v_bin, f.ϕ_bin); sleep(0.1)
    # @printf("\t%10.3f ≤ %10.3f < %10.3f : %s\n", _RoadFeature_BINS_T[f.t_bin], x[IND_T], _RoadFeature_BINS_T[f.t_bin+1], string(_RoadFeature_BINS_T[f.t_bin] ≤ x[IND_T] < _RoadFeature_BINS_T[f.t_bin+1])); sleep(0.1)
    # @printf("\t%10.3f ≤ %10.3f < %10.3f : %s\n", _RoadFeature_BINS_V[f.v_bin], x[IND_V], _RoadFeature_BINS_V[f.v_bin+1], string(_RoadFeature_BINS_V[f.v_bin] ≤ x[IND_V] < _RoadFeature_BINS_V[f.v_bin+1])); sleep(0.1)
    # @printf("\t%10.3f ≤ %10.3f < %10.3f : %s\n", _RoadFeature_BINS_T[f.ϕ_bin], x[IND_ϕ], _RoadFeature_BINS_ϕ[f.ϕ_bin+1], string(_RoadFeature_BINS_ϕ[f.ϕ_bin] ≤ x[IND_ϕ] < _RoadFeature_BINS_ϕ[f.ϕ_bin+1])); sleep(0.1)
    # println(""); sleep(0.1)

    convert(Float64, _RoadFeature_BINS_T[f.t_bin] ≤ x[IND_T] < _RoadFeature_BINS_T[f.t_bin+1] &&
                     _RoadFeature_BINS_V[f.v_bin] ≤ x[IND_V] < _RoadFeature_BINS_V[f.v_bin+1] &&
                     _RoadFeature_BINS_ϕ[f.ϕ_bin] ≤ x[IND_ϕ] < _RoadFeature_BINS_ϕ[f.ϕ_bin+1])
end

immutable RoadFeatureCont <: Feature
    t_pow::Float64
    v_pow::Float64
    ϕ_pow::Float64
end
scope(::RoadFeatureCont) = _RoadFeatureScope
standardize(v::Float64, μ::Float64, σ::Float64) = (v-μ)/σ
function evaluate(f::RoadFeatureCont, x::Vector{Float64})

    #=
    v:      -2.436      95.096      29.883      13.480
    t:     -17.100      11.889      -0.332       1.629
    ϕ:      -0.459       0.814      -0.001       0.021
    =#

    μ_t = -0.332
    μ_v = 29.883
    μ_ϕ = -0.001
    σ_t =  1.629
    σ_v = 13.480
    σ_ϕ =  0.021*5

    standardize(x[IND_T],μ_t,σ_t)^f.t_pow *
    standardize(x[IND_V],μ_v,σ_v)^f.v_pow *
    exp(-abs(standardize(x[IND_ϕ],μ_ϕ,σ_ϕ)^f.ϕ_pow))
end

immutable FollowFeature <: Feature
    Δs_bin::Int
    Δv_bin::Int
end
const _FollowFeature_BINS_ΔS = [-Inf, 0.0, 25.0, 50.0, 100.0, 200.0, Inf] # headway distance [ft]
const _FollowFeature_BINS_ΔV = [-Inf, -25.0, -10.0, -5.0, 0.0, 5.0, 10.0, 25.0, Inf] # velocity difference [ft/s]
const _FollowFeatureScope = @get_bitvector([:s1, :s2, :v1, :v2], 2)
scope(::FollowFeature) = _FollowFeatureScope
function evaluate(f::FollowFeature, x::Vector{Float64})

    Δs = x[IND_S+VEHICLE_FEATURE_SIZE] - x[IND_S]
    Δv = x[IND_V+VEHICLE_FEATURE_SIZE] - x[IND_V]

    convert(Float64, _FollowFeature_BINS_ΔS[f.Δs_bin] ≤ Δs < _FollowFeature_BINS_ΔS[f.Δs_bin+1] &&
                     _FollowFeature_BINS_ΔV[f.Δv_bin] ≤ Δv < _FollowFeature_BINS_ΔV[f.Δv_bin+1])
end

immutable FollowFeatureCont <: Feature
    Δs_pow::Float64
    Δv_pow::Float64
    T_pow::Float64 # timegap
end
scope(::FollowFeatureCont) = _FollowFeatureScope
function evaluate(f::FollowFeatureCont, x::Vector{Float64}, vehicles::Vector{Vehicle})

    veh_fore = vehicles[2]
    Δs = x[IND_S+VEHICLE_FEATURE_SIZE] - x[IND_S] - veh_fore.length
    Δv = x[IND_V+VEHICLE_FEATURE_SIZE] - x[IND_V]
    T = sqrt(clamp(Δs / max(x[IND_V],1.0), 0.0, 100.0) + 1e-6)

    #=
    Δd:      0.057     276.738      49.443      30.552
    Δv:    -51.560      50.177       0.184       4.203
    =#

    μ_Δs, σ_Δs = 49.443, 30.552
    μ_Δv, σ_Δv =  0.184,  4.203*10
    μ_T,  σ_T  =  1.373,  0.466

    standardize(Δs,μ_Δs,σ_Δs)^f.Δs_pow *
    standardize(Δv,μ_Δv,σ_Δv)^f.Δv_pow *
    standardize( T,μ_T, σ_T)^f.T_pow
end

immutable NeighborFeature <: Feature
    index::Int # index of neighbor
    mem::CPAMemory
end
const _NeighborFeatureScope = @get_bitvector([:t1, :t2, :ϕ1, :ϕ2, :s1, :s2, :v1, :v2], 2)
const NeighborFeature_MAX_INDEX = 5
scope(::NeighborFeature) = _NeighborFeatureScope
function get_relative_posG(state::VehicleState, s_fut, ϕ_fut, t_fut, roadway)

    C = state.posG
    F = state.posF

    Δs = s_fut - F.s
    Δϕ = ϕ_fut - F.ϕ

    curve = roadway.centerlines[F.laneid]
    extind_fut = move_extind_along(F.extind, curve, Δs)
    footpoint_fut = curve_at(curve, extind_fut)
    footpoint_fut.pos + polar(-t_fut, footpoint_fut.pos.θ - π/2) + VecSE2(0.0,0.0,Δϕ)
end
function evaluate(f::NeighborFeature, x::Vector{Float64}, vehicles::Vector{Vehicle}, o::Vector{VehicleState}, roadway::Roadway)

    posG1_fut = get_relative_posG(o[1], x[IND_S], x[IND_ϕ], x[IND_T], roadway)
    posG2_fut = get_relative_posG(o[2], x[VEHICLE_FEATURE_SIZE + IND_S], x[VEHICLE_FEATURE_SIZE + IND_ϕ], x[VEHICLE_FEATURE_SIZE + IND_T], roadway)

    vehA = vehicles[1]
    vehB = vehicles[2]
    state_origA = vehA.state
    state_origB = vehB.state
    vehA.state = VehicleState(posG1_fut, o[1].v)
    vehB.state = VehicleState(posG2_fut, o[2].v)

    t_CPA, d_CPA = get_time_and_dist_of_closest_approach(vehA, vehB, f.mem)

    vehA.state = state_origA
    vehB.state = state_origB

    if f.index == 1
        convert(Float64, t_CPA == 0.0 && d_CPA == 0.0)
    elseif f.index == 2
        convert(Float64, 0.0 < t_CPA ≤ 1.0 && d_CPA ≤ 1.64)
    elseif f.index == 3
        convert(Float64, 1.0 < t_CPA ≤ 4.0 && d_CPA ≤ 1.64)
    elseif f.index == 4
        convert(Float64, 4.0 < t_CPA ≤ 10.0 && d_CPA ≤ 1.64)
    elseif f.index == 5
        convert(Float64, 10.0 < t_CPA || d_CPA > 1.64)
    else
        warn("Unknown NeighborFactor Index: $(f.index)")
        0.0
    end
end

###############################
# Factor

type Factor
    name::Symbol
    features::Vector{Feature}         # a set of n features
    weights::Vector{Float64}          # a set of n weights
    feature_values::Vector{Float64}   # current feature value
    x::Vector{Float64}                # input vector, length is a multiple of VEHICLE_FEATURE_SIZE
    vehicles::Vector{Vehicle}         # set of original vehicles
    orig_states::Vector{VehicleState} # set of original vehicle states

    function Factor(name::Symbol, features::Vector{Feature})
        weights = rand(Float64, length(features))
        feature_values = Array(Float64, length(features))

        nvehicles = get_nvehicles(features[1])
        for f in features
            @assert(get_nvehicles(f) == nvehicles)
        end
        x = Array(Float64, nvehicles*VEHICLE_FEATURE_SIZE)
        vehicles = Array(Vehicle, nvehicles)
        orig_states = Array(VehicleState, nvehicles)

        retval = new()
        retval.name = name
        retval.features = features
        retval.weights = weights
        retval.feature_values = feature_values
        retval.x = x
        retval.vehicles = vehicles
        retval.orig_states = orig_states
        retval
    end
end
function insert_vehicle!(ϕ::Factor, veh::Vehicle, index::Int)
    i = (index-1)*VEHICLE_FEATURE_SIZE
    ϕ.x[i+IND_S] = veh.state.posF.s # s, position along lane
    ϕ.x[i+IND_T] = veh.state.posF.t # position perp to lane
    ϕ.x[i+IND_V] = veh.state.v
    ϕ.x[i+IND_ϕ] = veh.state.posF.ϕ # ϕ, lane-relative heading
    ϕ.orig_states[index] = veh.state
    ϕ.vehicles[index] = veh # NOTE: I assume that we are not changing its properties
    ϕ
end
function extract_features!(ϕ::Factor, roadway::Roadway)
    for i in 1 : length(ϕ.features)
        if isa(ϕ.features[i], RoadFeature)
            ϕ.feature_values[i] = evaluate(ϕ.features[i]::RoadFeature, ϕ.x)
        elseif isa(ϕ.features[i], RoadFeatureCont)
            ϕ.feature_values[i] = evaluate(ϕ.features[i]::RoadFeatureCont, ϕ.x)
        elseif isa(ϕ.features[i], FollowFeature)
            ϕ.feature_values[i] = evaluate(ϕ.features[i]::FollowFeature, ϕ.x)
        elseif isa(ϕ.features[i], FollowFeatureCont)
            ϕ.feature_values[i] = evaluate(ϕ.features[i]::FollowFeatureCont, ϕ.x, ϕ.vehicles)
        elseif isa(ϕ.features[i], NeighborFeature)
            ϕ.feature_values[i] = evaluate(ϕ.features[i]::NeighborFeature, ϕ.x, ϕ.vehicles, ϕ.orig_states, roadway)
        end
    end
    ϕ
end
set_variable!(ϕ::Factor, variable_index::Int, v::Float64) = ϕ.x[variable_index] = v
get_variable(ϕ::Factor, variable_index::Int) = ϕ.x[variable_index]
get_feature(ϕ::Factor, factor_index::Int) = ϕ.features[feature_index]
function factor_dot(ϕ::Factor)
    retval = 0.0
    for i in 1 : length(ϕ.weights)
        retval += ϕ.weights[i]*ϕ.feature_values[i]
    end
    retval
end
eval_factor(ϕ::Factor) = exp(factor_dot(ϕ))

###############################

function create_core_factors()
    retval = Factor[]

    road_features = Feature[]
    # for i in 1 : length(_RoadFeature_BINS_T)-1
    #     for j in 1 : length(_RoadFeature_BINS_V)-1
    #         for k in 1 : length(_RoadFeature_BINS_ϕ)-1
    #             push!(road_features, RoadFeature(i, j, k))
    #         end
    #     end
    # end
    max_pow = 3
    for i in 0:max_pow
        for j in 0:max_pow-i
            for k in 0:max_pow-i-j
                if !(i == j == k == 0)
                    push!(road_features, RoadFeatureCont(i*1.0,j*1.0,k*1.0))
                end
            end
        end
    end
    push!(retval, Factor(:road, road_features))

    follow_features = Feature[]
    # for i in 1 : length(_FollowFeature_BINS_ΔS)-1
    #     for j in 1 : length(_FollowFeature_BINS_ΔV)-1
    #         push!(follow_features, FollowFeature(i, j))
    #     end
    # end
    max_pow = 3
    for i in 0 : max_pow
        for j in 0 : max_pow-i
            for k in 0:max_pow-i-j
                if !(i == j == k == 0)
                    push!(follow_features, FollowFeatureCont(i*1.0,j*1.0,k*1.0))
                end
            end
        end
    end
    push!(retval, Factor(:follow, follow_features))

    neighbor_features = Feature[]
    for i in 1 : NeighborFeature_MAX_INDEX
        push!(neighbor_features, NeighborFeature(i, CPAMemory()))
    end
    push!(retval, Factor(:neighbor, neighbor_features))

    retval
end
function get_weight_dict(factors::Vector{Factor})
    retval = Dict{Symbol, Vector{Float64}}()
    for ϕ in factors
        retval[ϕ.name] = deepcopy(ϕ.weights)
    end
    retval
end
function set_weights!(factors::Vector{Factor}, weights::Dict{Symbol, Vector{Float64}})
    for (key, v) in weights

        factor_index = findfirst(ϕ->ϕ.name == key, factors)

        if factor_index != 0
            copy!(factors[factor_index].weights, v)
        end
    end
    factors
end