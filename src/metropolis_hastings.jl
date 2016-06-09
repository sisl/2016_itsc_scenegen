type FactorGraphSceneGenerator <: SceneGenerator
    factors::Vector{Factor}
    propsal_distribution::ContinuousMultivariateDistribution
    n_steps_metrohaste::Int
end

function shift_state(
    scene::Scene,
    structure::SceneStructure,
    veh_index::Int,
    Δ_propose::Vector{Float64},
    )

    s = scene.vehicles[veh_index].state
    curve = get_roadway(scene).centerlines[s.posF.laneid]
    new_extind = move_extind_along(s.posF.extind, curve, Δ_propose[IND_S])
    posF = Frenet(s.posF.laneid, new_extind, s.posF.s+Δ_propose[IND_S], s.posF.t + Δ_propose[IND_T], s.posF.ϕ + Δ_propose[IND_ϕ])

    posG = s.posG + polar(Δ_propose[IND_T], s.posG.θ + π/2)
    posG = posG + polar(Δ_propose[IND_S], s.posG.θ)
    posG = posG + VecSE2(0.0, 0.0, Δ_propose[IND_ϕ])

    v = s.v + Δ_propose[IND_V]

    VehicleState(posG, posF, v)
end
function adheres_to_structure(
    state_propose::VehicleState,
    veh_index::Int,
    scene::Scene,
    scene_structure::SceneStructure,
    )

    #=
    Returns whether the proposed state adheres to the given scene structure
    =#

    if state_propose.v      < BOUNDS_V_LO || state_propose.v      > BOUNDS_V_HI ||
       state_propose.posF.t < BOUNDS_T_LO || state_propose.posF.t > BOUNDS_T_HI ||
       state_propose.posF.ϕ < BOUNDS_ϕ_LO || state_propose.posF.ϕ > BOUNDS_ϕ_HI

        return false
    end

    # check whether proposed vehicle is between lead and lag
    veh = scene.vehicles[veh_index]
    state_orig = veh.state
    veh.state = state_propose
    ind_fore = scene_structure.index_fore[veh_index]
    ind_rear = scene_structure.index_rear[veh_index]
    fails_fore = ind_fore != 0 && get_headway_dist_between(veh, scene.vehicles[ind_fore]) ≤ 0.0
    fails_rear = ind_rear != 0 && get_headway_dist_between(scene.vehicles[ind_rear], veh) ≤ 0.0
    veh.state = state_orig

    # TODO: check for lateral collision?

    !(fails_fore || fails_rear)
end
function draw_active_vehicle_index(scene::Scene, structure::SceneStructure)
    veh_index = rand(1:length(structure.active_vehicles))
    for i in 1 : length(scene.vehicles)
        if in(i, structure.active_vehicles)
            veh_index -= 1
            if veh_index == 0
                veh_index = i
                break
            end
        end
    end
    veh_index
end

function evaluate(factors::Vector{Factor}, scene::Scene, structure::SceneStructure)

    roadway = get_roadway(scene)

    numerator = 0.0
    for (factor_index, vehicle_indeces) in structure.factor_assignments
        ϕ = factors[factor_index]
        pull_vehicle_data!(ϕ, scene, vehicle_indeces)
        extract_features!(ϕ, roadway)
        fd = factor_dot(ϕ)
        numerator += fd
    end

    exp(numerator)
end
function evaluate_subset(factors::Vector{Factor}, scene::Scene, structure::SceneStructure, veh_index::Int)
    #=
    Only evaluate the factors of which veh_index is a member
    =#

    roadway = get_roadway(scene)

    numerator = 0.0
    for (factor_index, vehicle_indeces) in structure.factor_assignments
        if in(veh_index, vehicle_indeces)
            ϕ = factors[factor_index]
            pull_vehicle_data!(ϕ, scene, vehicle_indeces)
            extract_features!(ϕ, roadway)
            fd = factor_dot(ϕ)
            numerator += fd
        end
    end

    exp(numerator)
end

function calc_acceptance_probability(
    scene::Scene,
    structure::SceneStructure,
    factors::Vector{Factor},
    veh_index::Int,
    state_propose::VehicleState,
    )

    if !adheres_to_structure(state_propose, veh_index, scene, structure)
        return 0.0 # do not accept out-of-bounds scenes
    end

    p_current = evaluate_subset(factors, scene, structure, veh_index)


    veh = scene.vehicles[veh_index]
    state_current = veh.state
    veh.state = state_propose
    p_propose = evaluate_subset(factors, scene, structure, veh_index)
    veh.state = state_current

    if !(p_propose ≥ 0.0)
        println("p_propose: ", p_propose)
    end
    @assert(p_propose ≥ 0.0)
    if p_propose == 0.0
        return 0.0
    end

    if p_current == 0.0
        # accept if the current one is not acceptable
        return 1.0
    end

    min(1.0, p_propose / p_current)
end
function metropolis_hastings_step!(
    scene::Scene,
    structure::SceneStructure,
    factors::Vector{Factor},
    propsal_distribution::ContinuousMultivariateDistribution,
    Δ_propose::Vector{Float64},
    )

    #=
    pick a random vehicle and shift it
    =#

    rand!(propsal_distribution, Δ_propose)

    veh_index = draw_active_vehicle_index(scene, structure)
    veh = scene.vehicles[veh_index]

    state_propose = shift_state(scene, structure, veh_index, Δ_propose)

    if rand() ≤ calc_acceptance_probability(scene, structure, factors, veh_index, state_propose)
        veh.state = state_propose
    end

    scene
end
function metropolis_hastings!(
    scene::Scene,
    structure::SceneStructure,
    factors::Vector{Factor},
    propsal_distribution::ContinuousMultivariateDistribution,
    n_steps::Int
    )

    #=
    Runs Metropolis Hastings for n steps
    =#

    Δ_propose = fill(NaN, VEHICLE_FEATURE_SIZE)

    for i in 1 : n_steps
        metropolis_hastings_step!(scene, structure, factors, propsal_distribution, Δ_propose)
    end
    scene
end

function generate_scene_and_structure(sg::FactorGraphSceneGenerator, dset::SceneDataset)
    starting_scene_index = rand(1:length(dset))
    scene_orig, structure_orig = get_scene_and_structure(dset, starting_scene_index)

    scene = deepcopy(scene_orig)
    source = dset.sources[starting_scene_index]
    structure = deepcopy(structure_orig)

    metropolis_hastings!(scene, structure, sg.factors, sg.propsal_distribution, sg.n_steps_metrohaste)

    (scene, source, structure)
end

function generate_dset(sg::SceneGenerator, dset::SceneDataset, n_scenes::Int)

    scenes = Array(Scene, n_scenes)
    sources = Array(SceneSource, n_scenes)
    structures = Array(SceneStructure, n_scenes)

    for i in 1 : n_scenes
        scene, source, structure = generate_scene_and_structure(sg, dset)
        scenes[i] = scene
        sources[i] = source
        structures[i] = structure
    end

    SceneDataset(scenes, sources, structures, deepcopy(dset.factors))
end