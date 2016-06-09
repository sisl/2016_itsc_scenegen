using MLBase

abstract OldSceneGenerator <: SceneGenerator

const VEHICLE_MEAN_LENGTH = 14.25 # [ft]
const VEHICLE_MEAN_WIDTH  =  6.15 # [ft]

####################

immutable StraightRoadway
    # NOTE(tim): all values should be positive
    nlanes       :: Int
    lanewidth    :: Float64
    horizon_fore :: Float64
    horizon_rear :: Float64

    function StraightRoadway(
        nlanes       :: Int     =   5,
        lanewidth    :: Float64 =  AVE_LANE_WIDTH,
        horizon_fore :: Float64 = 316.0,
        horizon_rear :: Float64 = 0.0
        )
        @assert(nlanes > 0)
        @assert(lanewidth > 0)
        @assert(horizon_fore > 0)
        @assert(horizon_rear ≥ 0)
        new(nlanes, lanewidth, horizon_fore, horizon_rear)
    end
end

lanecenters(road::StraightRoadway) = (collect(1:road.nlanes)-0.5)*road.lanewidth
roadwidth(road::StraightRoadway) = road.nlanes*road.lanewidth
roadlength(road::StraightRoadway) = road.horizon_rear + road.horizon_fore

####################

type OldVehicle
    x :: Float64
    y :: Float64
    ϕ :: Float64
    v :: Float64

    class  :: Int # ∈ (1-motorcycle, 2-auto, 3-truck)
    length :: Float64
    width  :: Float64

    OldVehicle() = new()
    function OldVehicle(x::Float64, y::Float64, ϕ::Float64=0.0, v::Float64=0.0;
        class  :: Int     = CLASS_AUTOMOBILE,
        length :: Float64 = VEHICLE_MEAN_LENGTH,
        width  :: Float64 = VEHICLE_MEAN_WIDTH
        )
        new(x, y, ϕ, v, class, length, width)
    end
end
Base.show(io::IO, veh::OldVehicle) = @printf(io, "OldVehicle(x=%.3f, y=%.3f, ϕ=%.4f, v=%.3f, class=%d, length=%.3f, width=%.3f)",
                                             veh.x, veh.y, veh.ϕ, veh.v, veh.class, veh.length, veh.width)
Base.print(io::IO, veh::OldVehicle) = show(io, veh)



type VehicleExtraInfo
    laneindex   :: Int # what lane it is assigned to
    laneorder   :: Int # what index within the lane it is, increasing in longitudinal direction
    d_cl        :: Float64
    carind_fore :: Int
    carind_rear :: Int
    fake_x      :: Float64

    VehicleExtraInfo() = new()
    VehicleExtraInfo(laneindex::Int, laneorder::Int, d_cl::Float64, carind_fore::Int, carind_rear::Int) =
        new(laneindex, laneorder, d_cl, carind_fore, carind_rear, 0.0)
end

Base.show(io::IO, i::VehicleExtraInfo) = @printf(io, "VehicleExtraInfo(%d, %d, %.3f, %d, %d, %.3f)",
                                             i.laneindex, i.laneorder, i.d_cl, i.carind_fore, i.carind_rear, i.fake_x)
Base.print(io::IO, i::VehicleExtraInfo) = show(io, i)


type RoadScene
    road     :: StraightRoadway
    vehicles :: Vector{OldVehicle}
    info     :: Vector{VehicleExtraInfo}

    lanespeeds   :: Vector{Float64}
    lanedensities :: Vector{Float64}

    RoadScene(road::StraightRoadway, vehicles::Vector{OldVehicle}, info::Vector{VehicleExtraInfo} = VehicleExtraInfo[]) =
        new(road, vehicles, info, Float64[], Float64[])
end

if !isdefined(:ModelOptimizationResults)
    type ModelOptimizationResults{T<:SceneGenerator}
        params     :: Dict{Symbol, Any}
        logl_mean  :: Float64
        logl_stdev :: Float64
        score      :: Float64
        iter       :: Int
        CV_nfolds  :: Int
        CV_rounds  :: Int
    end
end

#####################

function calc_lanespeeds(scene::RoadScene)
    nlanes = scene.road.nlanes
    retval = zeros(Float64, nlanes)
    counts = zeros(Int, nlanes)

    for (i,veh) in enumerate(scene.vehicles)
        lane = scene.info[i].laneindex
        retval[lane] += veh.v
        counts[lane] += 1
    end

    retval ./ counts
end
function calc_lanedensities(scene::RoadScene)

    nlanes = scene.road.nlanes
    counts = zeros(Int, nlanes)

    for (i,info) in enumerate(scene.info)
        counts[info.laneindex] += 1
    end

    counts / roadlength(scene.road) # [vehicles / ft]
end
function calc_d_front(scene::RoadScene, carind::Int)
    @assert(scene.info[carind].carind_fore != 0)
    veh_front = scene.vehicles[scene.info[carind].carind_fore]
    veh_front.y - scene.vehicles[carind].y - veh_front.length
end
function calc_d_rear(scene::RoadScene, carind::Int)
    @assert(scene.info[carind].carind_rear != 0)
    veh_rear = scene.vehicles[scene.info[carind].carind_rear]
    veh = scene.vehicles[carind]
    veh.y - veh_rear.y - veh.length
end

function calc_laneindex(centers::Vector{Float64}, x::Float64)
    bestind = 0
    bestval = Inf
    for (i,center) in enumerate(centers)
        Δ = x - center
        if abs(Δ) < bestval
            bestval, bestind = abs(Δ), i
        end
    end
    bestind
end
function calc_vehicle_extra_info!(scene::RoadScene)

    nlanes = scene.road.nlanes
    centers = lanecenters(scene.road)

    vehicles = scene.vehicles
    nvehicles = length(vehicles)

    lanesets = Array(Set{Int}, nlanes)
    for i = 1 : nlanes
        lanesets[i] = Set{Int}()
    end

    sceneinfo = Array(VehicleExtraInfo, nvehicles)

    for (i,veh) in enumerate(vehicles)
        laneindex = calc_laneindex(centers, veh.x)
        d_cl      = veh.x - centers[laneindex]
        sceneinfo[i] = VehicleExtraInfo(laneindex, 0, d_cl, 0, 0)
        push!(lanesets[laneindex],i)
    end

    for i = 1 : nlanes
        carinds  = collect(lanesets[i])
        ncarinds = length(carinds)
        arr_y    = Array(Float64, ncarinds)
        for (j, carind) in enumerate(carinds)
            arr_y[j] = vehicles[carind].y
        end
        p = sortperm(arr_y)
        for laneorder = 1 : ncarinds
            carind = carinds[p[laneorder]]
            veh = vehicles[carind]
            info = sceneinfo[carind]
            info.laneorder = laneorder
            if laneorder > 1
                info.carind_rear = carinds[p[laneorder-1]]

                veh_rear = vehicles[sceneinfo[carind].carind_rear]
                info_rear = sceneinfo[sceneinfo[carind].carind_rear]
                # @assert(info_rear.carind_fore == carind)
                @assert(info_rear.laneorder == laneorder-1)
                @assert(veh_rear.y < veh.y)
                @assert(info_rear.laneindex == info.laneindex)
            end
            if laneorder < ncarinds
                info.carind_fore = carinds[p[laneorder+1]]

                veh_fore = vehicles[sceneinfo[carind].carind_fore]
                info_fore = sceneinfo[sceneinfo[carind].carind_fore]
                @assert(veh_fore.y > veh.y)
                @assert(info_fore.laneindex == info.laneindex)
            end
        end
    end

    scene.info = sceneinfo
    scene.lanespeeds = calc_lanespeeds(scene)
    scene.lanedensities = calc_lanedensities(scene)

    sceneinfo
end
function calc_sorted_cars_in_lane(scene::RoadScene, lane::Int)
    n_cars = 0
    for info in scene.info
        if info.laneindex == lane
            n_cars += 1
        end
    end

    if n_cars == 0
        return Int[]
    end

    retval = Array(Int, n_cars)
    for (i,info) in enumerate(scene.info)
        if info.laneindex == lane
            retval[info.laneorder] = i
        end
    end
    retval
end

function get_scenes(sg::SceneGenerator, road::StraightRoadway, nscenes::Int)
    retval = Array(RoadScene, nscenes)
    for i = 1 : nscenes
        retval[i] = generate_scene(sg, road)
    end
    retval
end
function loglikelihood(
    sg     :: SceneGenerator,
    scenes :: Vector{RoadScene}
    )

    retval = 0.0
    for scene in scenes
        retval += loglikelihood(sg, scene)
    end
    retval
end
function ave_crossvalidated_likelihood(scores::Matrix{Float64})
    arr = mean(scores, 1)
    μ = mean(arr)
    σ = stdm(arr, μ)
    (μ,σ)
end

function get_base_arrays(scenes::Vector{RoadScene})
    tot_nvehicles = 0
    for scene in scenes
        tot_nvehicles += length(scene.vehicles)
    end

    arr_v       = fill(NaN, tot_nvehicles)
    arr_d_front = fill(NaN, tot_nvehicles)
    arr_d_cl    = fill(NaN, tot_nvehicles)
    arr_yaw     = fill(NaN, tot_nvehicles)

    count = 0
    count_d_front = 0
    for scene in scenes
        for (i,veh) in enumerate(scene.vehicles)
            count += 1
            info = scene.info[i]
            arr_v[count]    = veh.v
            arr_d_cl[count] = info.d_cl
            arr_yaw[count]  = veh.ϕ
            if info.carind_fore != 0
                arr_d_front[count_d_front+=1] = calc_d_front(scene, i)
                if !(arr_d_front[count_d_front] > 0.0)
                    println("rear: ", veh)
                    println("      ", info)
                    println("fore: ", scene.vehicles[info.carind_fore])
                    println("      ", scene.info[info.carind_fore])
                end
                @assert(arr_d_front[count_d_front] > -scene.vehicles[info.carind_fore].length)
            end
        end
    end
    @assert(count == tot_nvehicles)

    println("v: ", extrema(arr_v))
    println("d_front: ", extrema(arr_d_front))
    println("t: ", extrema(arr_d_cl))
    println("ϕ: ", extrema(arr_yaw))

    # arr_v         += (rand(tot_nvehicles)-0.5)*0.02
    # arr_d_front += (rand(tot_nvehicles)-0.5)*0.02
    # arr_d_cl      += (rand(tot_nvehicles)-0.5)*0.02
    # arr_yaw       += (rand(tot_nvehicles)-0.5)*0.02

    (arr_v, arr_d_front[1:count_d_front], arr_d_cl, arr_yaw)
end