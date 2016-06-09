type SceneExtractParams
    center::VecSE2 # center position of scene
    length::Float64
    box::ConvexPolygon

    function SceneExtractParams(
        center::VecSE2,
        length::Float64=300.0, # distance along scene orientation to extend the scene to [ft]
        width::Float64=200.0, # distance to either side to extend the scene to [ft]
        )

        x = polar(length/2, center.θ)
        y = polar(width/2, center.θ+π/2)

        o = convert(VecE2, center)
        box = ConvexPolygon(4)
        push!(box, o + x - y)
        push!(box, o + x + y)
        push!(box, o - x + y)
        push!(box, o - x - y)

        ensure_pts_sorted_by_min_polar_angle!(box)

        new(center, length, box)
    end
end

const REGIONS = Dict{ASCIIString, SceneExtractParams}(
  "80A"  => SceneExtractParams(VecSE2(6042777.824, 2133302.509, 1.684)),
  "101B" => SceneExtractParams(VecSE2(6451853.000, 1872681.000, 2.435)),
  "101A" => SceneExtractParams(VecSE2(6451293.000, 1873251.000, 2.382)),
  "80B"  => SceneExtractParams(VecSE2(6042697.824, 2134372.509, 1.745)),
  "101C" => SceneExtractParams(VecSE2(6452453.000, 1872211.000, 2.417)),
  )

function is_in_bounds(scene::Scene, scene_params::SceneExtractParams)

    # at least one vehicle has passed the scene and at least one vehicle has yet to enter it

    max_dist_front = 0.0
    max_dist_rear = 0.0

    for i in 1 : length(scene)
        veh = scene.vehicles[i]
        p_rel = inertial2body(veh.state.posG, scene_params.center)
        max_dist_front = max(max_dist_front, p_rel.x)
        max_dist_rear = min(max_dist_rear, p_rel.x - veh.length)
    end

    max_dist_front ≥ scene_params.length/2 &&
    -max_dist_rear ≥ scene_params.length/2
end
function is_collision_free(scene::Scene, mem::CPAMemory=CPAMemory())

    for i in 1 : length(scene)-1
        set_to_positioned_oriented_bounding_box!(mem.vehA, scene[i])

        for j in i+1 : length(scene)
            set_to_positioned_oriented_bounding_box!(mem.vehB, scene[j])

            if is_colliding(mem.vehA, mem.vehB, mem.mink)
                return false
            end
        end
    end

    true
end
function is_collision_free_to_horizon(scene::Scene, rec::SceneRecord, frames_per_tick::Int, mem::CPAMemory=CPAMemory();
    model::DriverModel = IntelligentDriverModel(σ₁=1e-4, σ₂=1e-6),
    )

    empty!(rec)
    drivers = _get_drivers_dict(model, scene)
    simulate_into_record!(rec, scene, drivers, frames_per_tick=frames_per_tick)

    isnan(get_time_of_first_collision(rec, mem))
end
function is_there_longitudinal_room(scene::Scene)

    #=
    True if there is is always headway separation between all vehicles
    =#

    roadway = get_roadway(scene)

    for vehicle_index in 1 : length(scene)-1

        veh = scene.vehicles[vehicle_index]

        ind2 = get_neighbor_index_fore(scene, vehicle_index)
        if ind2 != 0
            veh_fore = scene.vehicles[ind2]
            dist = get_headway_dist_between(veh, veh_fore)
            if !isnan(dist) && dist < 0.0
                return false
            end
        end

        ind2 = get_neighbor_index_rear(scene, vehicle_index)
        if ind2 != 0
            veh_rear = scene.vehicles[ind2]
            dist = get_headway_dist_between(veh_rear, veh)
            if !isnan(dist) && dist < 0.0
                return false
            end
        end
    end

    true
end
function is_scene_well_behaved(
    scene_params::SceneExtractParams,
    scene::Scene,
    rec::SceneRecord,
    frames_per_tick::Int,
    mem::CPAMemory=CPAMemory(),
    )

    is_in_bounds(scene, scene_params) &&
        is_there_longitudinal_room(scene) &&
        # is_collision_free(scene, mem) &&
        is_collision_free_to_horizon(scene, rec, frames_per_tick, mem)
end