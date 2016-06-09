const AVE_LANE_WIDTH = 11.44 # [ft]

immutable SceneLaneTranslationInfo
    laneid::Int
    curveind_lo::Int
    curveind_hi::Int
end
type SceneTranslationInfo
    extract::SceneExtractParams
    road::StraightRoadway
    roadway_name::Symbol
    lanes::Vector{SceneLaneTranslationInfo}
end
function get_lane_translation_info(sti::SceneTranslationInfo, laneid::Int)
    for s in sti.lanes
        if s.laneid == laneid
            return s
        end
    end

    error("lane translation info not found for laneid $laneid, $(sti.lanes)") # DEBUG
    sti.lanes[1]
end

function get_scene_translation_info(extract::SceneExtractParams, roadway::Roadway)
    #=
    1 - find all centerlines that intersect with the scene_extract box
    2 - for each, get its SceneLaneTranslationInfo
    3 - construct a StraightRoadway
    =#

    lanes = SceneLaneTranslationInfo[]

    for (laneid, centerline) in enumerate(roadway.centerlines)

        curveind_lo = -1
        curveind_hi = -1

        for i in 1 : length(centerline)
            if contains(extract.box, centerline[i].pos)

                if curveind_lo == -1
                    # we found the first extind that is in the box
                    curveind_lo = i
                else
                    # update curvein_hi
                    curveind_hi = i
                end
            end
        end

        if curveind_lo != -1

            if curveind_hi == -1
                curveind_hi = length(centerline)
            end

            # back up 10 ft
            s₀ = centerline[curveind_lo].s
            while curveind_lo > 1 && s₀ - centerline[curveind_lo-1].s < 10.0
                curveind_lo -= 1
            end
            # advance 10 ft
            s₀ = centerline[curveind_hi].s
            while curveind_hi < length(centerline) && centerline[curveind_hi+1].s - s₀ < 10.0
                curveind_hi += 1
            end

            push!(lanes, SceneLaneTranslationInfo(laneid, curveind_lo, curveind_hi))
        end
    end

    road = StraightRoadway(length(lanes))
    SceneTranslationInfo(extract, road, roadway.name, lanes)
end
function get_scene_translation_infos()
    retval = Dict{ASCIIString, SceneTranslationInfo}()
    for (key, extract) in REGIONS
        if contains(key, "101")
            retval[key] = get_scene_translation_info(extract, NGSIM.ROADWAY_101)
        else
            retval[key] = get_scene_translation_info(extract, NGSIM.ROADWAY_80)
        end
    end
    retval
end

OLD_ROADWAYS = get_scene_translation_infos()

function get_scene_translation_info(scene::Scene)

    veh = scene[1]
    for sti in values(OLD_ROADWAYS)
        if contains(sti.extract.box, veh.state.posG)
            return sti
        end
    end

    error("scene translation info not found")
    OLD_ROADWAYS[1]
end

function Base.convert(::Type{RoadScene}, scene::Scene)
    roadway = get_roadway(scene)
    sti = get_scene_translation_info(scene)
    road = sti.road

    vehicles = Array(OldVehicle, length(scene.vehicles))
    info = Array(VehicleExtraInfo, length(vehicles))

    for (i,veh) in enumerate(scene)

        carind_fore = get_neighbor_index_fore(scene, i)
        carind_rear = get_neighbor_index_rear(scene, i)

        if carind_fore != 0
            @assert(scene[carind_fore].state.posF.laneid == veh.state.posF.laneid)
        end
        if carind_rear != 0
            @assert(scene[carind_rear].state.posF.laneid == veh.state.posF.laneid)
        end

        laneorder = 1
        index_target = carind_rear
        index_of_last_car = carind_rear == 0 ? i : carind_rear
        while index_target != 0
            laneorder += 1
            index_of_last_car = index_target
            index_target = get_neighbor_index_rear(scene, index_target)
        end

        # NOTE: x is lateral and y is longitudinal
        # will set all lanes such that s₀ = 0

        lti = get_lane_translation_info(sti, veh.state.posF.laneid)
        s_lo = roadway.centerlines[lti.laneid][lti.curveind_lo].s
        x = veh.state.posF.t
        y = veh.state.posF.s - s_lo
        if y < 0.0
            println("y: ", y)
            println("posF: ", veh.state.posF)
            println("s_lo: ", s_lo)
        end
        @assert(y ≥ 0.0)

        @assert(!isnan(veh.state.posF.s))
        @assert(!isnan(veh.state.posF.t))
        @assert(!isnan(veh.state.posF.ϕ))
        vehicles[i] = OldVehicle(x, y, veh.state.posF.ϕ, veh.state.v,
                                 class=veh.class, length=veh.length, width=veh.width)

        info[i] = VehicleExtraInfo(veh.state.posF.laneid, laneorder, veh.state.posF.t, carind_fore, carind_rear)
    end

    RoadScene(road, vehicles, info)
end
function Base.convert(::Type{Scene}, roadscene::RoadScene, sti::SceneTranslationInfo)

    roadway = NGSIM.ROADWAY_DICT[sti.roadway_name]
    vehicles = Array(Vehicle, length(roadscene.vehicles))

    for i in 1:length(roadscene.vehicles)

        veh_old = roadscene.vehicles[i]
        info = roadscene.info[i]

        veh = Vehicle()
        veh.id = -999
        veh.class = veh_old.class
        veh.length = veh_old.length
        veh.width = veh_old.width

        lti = get_lane_translation_info(sti, info.laneindex)
        curve = roadway.centerlines[lti.laneid]

        s = veh_old.y

        extind = get_extind(curve, s)

        posF = Frenet(info.laneindex, extind, s, info.d_cl, veh_old.ϕ)
        footpoint = curve_at(curve, extind).pos
        posG = VecSE2(convert(VecE2, footpoint) + polar(posF.t, footpoint.θ+π/2), footpoint.θ + posF.ϕ)

        veh.state = VehicleState(posG, posF, veh_old.v)

        vehicles[i] = veh
    end

    retval = Scene(sti.roadway_name, vehicles)
end

function get_roadscene_dset(dset::SceneDataset)
    roadscenes = Array(RoadScene, length(dset.scenes))
    for (i,scene) in enumerate(dset.scenes)
        roadscenes[i] = convert(RoadScene, scene)
    end
    roadscenes
end

function generate_scene(sg::OldSceneGenerator, dset::SceneDataset)

    # 1 - sample a random scene / structure
    scene_index = rand(1:length(dset))
    scene, structure = get_scene_and_structure(dset, scene_index)
    source = dset.sources[scene_index]

    # 2 - generate the corresponding StraightRoadway
    sti = get_scene_translation_info(scene)

    # 3 - sample a RoadScene from the scene generator
    roadscene = generate_scene(sg, sti)

    # 4 - convert to a Scene
    ret_scene = convert(Scene, roadscene, sti)
    (ret_scene, source)
end
function generate_scene_and_structure(sg::OldSceneGenerator, dset::SceneDataset)
    scene, source = generate_scene(sg, dset)
    structure = gen_scene_structure(scene, dset.factors)
    (scene, source, structure)
end