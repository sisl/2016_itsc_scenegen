using PGFPlots
using DataFrames
using ProfileView
using Interact
using Discretizers
using Distributions
using OnlineStats
using NGSIM

using JLD
using Vec
using Cairo
using Reel
Reel.set_output_type("gif")
Reel.extension(m::MIME"image/svg+xml") = "svg"

const SAVE_DIR = "/home/tim/Documents/papers/2016_itsc_scenegen_wheeler/code/full_scene/output"
const FLOATING_POINT_REGEX = r"(([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)|Inf)"
const FRAME_TIME_STEP = 0.1 # frame step rate

if !isdefined(:Renderer)
    include("viz/Renderer.jl")
end
if !isdefined(:ColorScheme)
    include("viz/ColorScheme.jl")
    pushPGFPlotsPreamble("\\pgfplotsset{compat=1.10}")
    pushPGFPlotsPreamble("\\usepgfplotslibrary{fillbetween}")
    define_color("monokai1", [0xCF/0xFF,0xBF/0xFF,0xAD/0xFF])
    define_color("monokai2", [0x27/0xFF,0x28/0xFF,0x22/0xFF])
    define_color("monokai3", [0x52/0xFF,0xE3/0xFF,0xF6/0xFF])
    define_color("monokai4", [0xA7/0xFF,0xEC/0xFF,0x21/0xFF])
    define_color("monokai5", [0xFF/0xFF,0x00/0xFF,0x7F/0xFF])
    define_color("monokai6", [0xF9/0xFF,0x97/0xFF,0x1F/0xFF])
    define_color("monokai7", [0x79/0xFF,0xAB/0xFF,0xFF/0xFF])
end

using Renderer
using ColorScheme

###########################

include("probdrive/standardization.jl")
include("probdrive/driver_models.jl")

###########################

include("minkowski.jl")
include("factor.jl")
include("scene_extract_params.jl")

###########################

include("probdrive/simulate.jl")
include("probdrive/features_scenegen.jl")
include("probdrive/linear_bayesian.jl")
include("probdrive/IDM.jl")

###########################

function pull_subscene(trajdata::Trajdata, frame::Int, scene_params::SceneExtractParams, extractor::Scene)

    get!(extractor, trajdata, frame)

    vehicles = Vehicle[]
    for i in 1 : length(extractor)
        veh = extractor.vehicles[i]
        if contains(scene_params.box, veh.state.posG)
            push!(vehicles, deepcopy(veh))
        end
    end

    Scene(trajdata.roadway.name, vehicles)
end
function pull_subscenes(trajdata::Trajdata, scene_params::SceneExtractParams, extractor::Scene=Scene())

    # pull all of the scenes

    subscenes = Scene[]

    for frame in 1 : nframes(trajdata)
        if is_scene_in_bounds(trajdata, frame, scene_params)
            scene = pull_subscene(trajdata, frame, scene_params, extractor)
            push!(subscenes, scene)
        end
    end

    subscenes
end
function pull_subscenes(trajdata::Trajdata, scene_extract_params_arr::Vector{SceneExtractParams}, extractor::Scene=Scene())

    subscenes = Scene[]

    for scene_params in scene_extract_params_arr
        append!(subscenes, pull_subscenes(trajdata, scene_params, extractor))
    end

    subscenes
end

##############

type SceneStructure
    nvehicles::Int
    factor_assignments::Vector{Tuple{Int, Vector{Int}}} # list of factor index -> vehicles
    active_vehicles::Set{Int} # set of vehicles that can be manipulated (vehicle index)
    index_fore::Vector{Int} # maps vehicle index to the index of its leader, 0 if none
    index_rear::Vector{Int} # maps vehicle index to the index of its follower, 0 if none
end
_ordered_tup_pair(a::Int, b::Int) = a < b ? (a,b) : (b,a)
function gen_scene_structure(scene::Scene, factors::Vector{Factor})

    nvehicles = length(scene)
    factor_assignments = Tuple{Int, Vector{Int}}[]
    active_vehicles = Set{Int}()
    index_fore = Array(Int, nvehicles)
    index_rear = Array(Int, nvehicles)

    # add all road factors
    # - these only get added to active cars, cars which have a lead and a rear

    factor_index_road = findfirst(f->f.name == :road, factors)
    factor_index_follow = findfirst(f->f.name == :follow, factors)
    factor_index_neighbor = findfirst(f->f.name == :neighbor, factors)

    pair_factors_added = Set{Tuple{Int, Int}}() # (a,b) s.t. a < b

    for vehicle_index in 1 : nvehicles

        veh_fore_index = get_neighbor_index_fore(scene, vehicle_index)
        veh_rear_index = get_neighbor_index_rear(scene, vehicle_index)

        index_fore[vehicle_index] = veh_fore_index
        index_rear[vehicle_index] = veh_rear_index

        if veh_fore_index != 0 && veh_rear_index != 0
            push!(factor_assignments, (factor_index_road, Int[vehicle_index]))
            push!(active_vehicles, vehicle_index)

            tup = _ordered_tup_pair(vehicle_index, veh_fore_index)
            if !in(tup, pair_factors_added)
                push!(factor_assignments, (factor_index_follow, Int[vehicle_index, veh_fore_index]))
                push!(pair_factors_added, tup)
            end

            tup = _ordered_tup_pair(vehicle_index, veh_rear_index)
            if !in(tup, pair_factors_added)
                push!(factor_assignments, (factor_index_follow, Int[veh_rear_index, vehicle_index]))
                push!(pair_factors_added, tup)
            end
        end
    end

    # add all neighbor features
    try_to_add_neighbor_factor = (vehicle_index, ind) -> begin
        if ind != 0
            tup = _ordered_tup_pair(vehicle_index, ind)
            if !in(tup, pair_factors_added)
                push!(factor_assignments, (factor_index_neighbor, Int[vehicle_index, ind]))
                push!(pair_factors_added, tup)
            end
        end
        nothing
    end

    for vehicle_index in active_vehicles
        try_to_add_neighbor_factor(vehicle_index, get_neighbor_index_right(scene, vehicle_index, 0.0, 100.0))
        try_to_add_neighbor_factor(vehicle_index, get_neighbor_index_right(scene, vehicle_index, -100.0, 0.0))
        try_to_add_neighbor_factor(vehicle_index, get_neighbor_index_left(scene, vehicle_index, 0.0, 100.0))
        try_to_add_neighbor_factor(vehicle_index, get_neighbor_index_left(scene, vehicle_index, -100.0, 0.0))
    end

    SceneStructure(nvehicles, factor_assignments, active_vehicles, index_fore, index_rear)
end

###############

immutable SceneSource
    trajdata_id::Int
    frame::Int
end

type SceneDataset
    scenes::Vector{Scene}
    sources::Vector{SceneSource}
    structures::Vector{SceneStructure}
    factors::Vector{Factor}
end
Base.length(dset::SceneDataset) = length(dset.scenes)
function nvehicles(dset::SceneDataset)
    count = 0
    for scene in dset.scenes
        count += length(scene)
    end
    count
end

function load_trajdatas(re_extract::Bool=false)

    trajdatas = Dict{Int, Trajdata}()

    for (path_index, trajdata_input_path) in enumerate(NGSIM.TRAJDATA_INPUT_PATHS)

        if re_extract
            trajdatas[path_index] = Trajdata(trajdata_input_path, path_index)
        else
            extracted_trajdata_csv = input_path_to_extracted_trajdata_csv(trajdata_input_path)
            trajdatas[path_index] = Trajdata(extracted_trajdata_csv, path_index)
        end
    end

    trajdatas
end

function get_scene_and_structure(dset::SceneDataset, index::Int)
    scene = dset.scenes[index]
    structure = dset.structures[index]
    (scene, structure)
end
function pull_scene_dataset(
    trajdata::Trajdata,
    scene_params::SceneExtractParams,
    factors::Vector{Factor};
    extractor::Scene=Scene(),
    max_sample_size::Int=typemax(Int),
    scene_skip::Int=1, # number of frames to skip between scenes
    frames_per_tick::Int = 3,
    rec::SceneRecord = SceneRecord(frames_per_tick*10+1),
    mem::CPAMemory=CPAMemory(),
    )

    scenes = Scene[]
    sources = SceneSource[]
    structures = SceneStructure[]

    for frame in 1 : scene_skip : nframes(trajdata)
        scene = pull_subscene(trajdata, frame, scene_params, extractor)

        if is_scene_well_behaved(scene_params, extractor, rec, frames_per_tick, mem)

            source = SceneSource(trajdata.id, frame)
            structure = gen_scene_structure(scene, factors)

            if !isempty(structure.factor_assignments) # we can learn from it
                push!(scenes, scene)
                push!(sources, source)
                push!(structures, structure)
            end
        end

        if length(scenes) ≥ max_sample_size
            break
        end
    end

    SceneDataset(scenes, sources, structures, factors)
end

function reset_weights!(dset::SceneDataset)
    for ϕ in dset.factors
        for j in 1 : length(ϕ.weights)
            ϕ.weights[j] = 0.01*rand()
        end
    end
    dset
end
function pull_vehicle_data!(ϕ::Factor, scene::Scene, vehicle_indeces::AbstractVector{Int})
    for i in 1 : length(vehicle_indeces)
        vehicle_index = vehicle_indeces[i]
        veh = scene.vehicles[vehicle_index]
        insert_vehicle!(ϕ, veh, i)
    end
    ϕ
end

include("dataset_export.jl")

function pull_exportable_subscene_dataset(;
    max_sample_size_per_trajdata::Int = typemax(Int),
    scene_skip::Int = 10,
    )

    core_factors = create_core_factors()
    dsetexp = DatasetExport()
    extractor = Scene()
    mem = CPAMemory()

    for trajdata_id in 1 : length(NGSIM.TRAJDATA_INPUT_PATHS)

        trajdata = Trajdata(trajdata_id)
        csvfile = NGSIM.TRAJDATA_INPUT_PATHS[trajdata_id]

        for (key, extract) in REGIONS
            if (contains(csvfile, "101") && contains(key, "101")) ||
               (contains(csvfile, "80") && contains(key, "80"))

                println("pulling scenes for ", key)
                tic()
                dset = pull_scene_dataset(trajdata, extract, core_factors,
                            max_sample_size=max_sample_size_per_trajdata, scene_skip=scene_skip,
                            extractor=extractor, mem=mem)

                append!(dsetexp, dset, trajdata_id)
                toc()
            end
        end
    end

    dsetexp
end
function get_most_resent_file(dir::AbstractString, regex::Regex=r"full_scene_(\d+).jld")
    best = "none"
    best_time = 0.0
    for file in readdir(dir)
        file = joinpath(dir, file)
        if isfile(file) && ismatch(regex, file)
            stats = stat(file)
            if stats.mtime > best_time
                best_time = stats.mtime
                best = file
            end
        end
    end
    best
end

type OnlineExtrema
    lo::Float64
    hi::Float64
    μ::Float64
    M::Float64
    n::Int
    OnlineExtrema() = new(Inf, -Inf, 0.0, 0.0, 0)
end
function update!(o::OnlineExtrema, v::Float64)
    o.lo = min(o.lo, v)
    o.hi = max(o.hi, v)

    n, μ, M = o.n, o.μ, o.M

    n += 1
    μ_next = μ + (v - μ)/n
    M += (v - μ)*(v - μ_next)

    o.n, o.μ, o.M = n, μ_next, M

    o
end
function get_dset_bounds(dset::SceneDataset)

    v = OnlineExtrema()
    t = OnlineExtrema()
    ϕ = OnlineExtrema()
    Δd = OnlineExtrema()
    Δv = OnlineExtrema()
    T = OnlineExtrema()

    for scene in dset.scenes
        for (i,veh) in enumerate(scene)
            s = veh.state
            update!(v, s.v)
            update!(t, s.posF.t)
            update!(ϕ, s.posF.ϕ)

            ind_fore = get_neighbor_index_fore(scene, i)
            if ind_fore != 0
                Δs = get_headway_dist_between(scene[i], scene[ind_fore])
                update!(Δd, Δs)
                update!(Δv, scene[ind_fore].state.v - s.v)
                update!(T, sqrt(clamp(Δs / max(s.v,1.0), 0.0, 100.0) + 1e-6))
            end
        end
    end

    @printf("v:  %10.3f  %10.3f  %10.3f  %10.3f\n",  v.lo,  v.hi,  v.μ, sqrt( v.M/( v.n - 1)))
    @printf("t:  %10.3f  %10.3f  %10.3f  %10.3f\n",  t.lo,  t.hi,  t.μ, sqrt( t.M/( t.n - 1)))
    @printf("ϕ:  %10.3f  %10.3f  %10.3f  %10.3f\n",  ϕ.lo,  ϕ.hi,  ϕ.μ, sqrt( ϕ.M/( ϕ.n - 1)))
    @printf("Δd: %10.3f  %10.3f  %10.3f  %10.3f\n", Δd.lo, Δd.hi, Δd.μ, sqrt(Δd.M/(Δd.n - 1)))
    @printf("Δv: %10.3f  %10.3f  %10.3f  %10.3f\n", Δv.lo, Δv.hi, Δv.μ, sqrt(Δv.M/(Δv.n - 1)))
    @printf("T:  %10.3f  %10.3f  %10.3f  %10.3f\n",  T.lo,  T.hi,  T.μ, sqrt( T.M/( T.n - 1)))

    (v, t, ϕ, Δd, Δv, T)
end


# ################

abstract SceneGenerator

include("learning.jl")
include("weighted_GMM.jl")
include("metropolis_hastings.jl")
include("sampling.jl")
include("probdrive/sample_and_prop.jl")

# ################

include("old_models/old_trajdata.jl")
include("old_models/old_new_interface.jl")
include("old_models/univariate_scene_generator.jl")
include("old_models/joint_bn_chain_scene_generator.jl")

# ################

include("viz/viz.jl")
include("viz/viz_factor.jl")
include("viz/viz_convergence.jl")




