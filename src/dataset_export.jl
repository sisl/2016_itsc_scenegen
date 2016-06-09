immutable ExportVehicleDescription
    class::Int # ∈ (1-motorcycle, 2-auto, 3-truck)
    length::Float64
    width::Float64
end

immutable ExportState
    id::Int # vehicle ID
    state::VehicleState
end

immutable ExportScene
    lo::Int # index in states of first vehicle in scene
    hi::Int # index in states of last vehicle in scene
    trajdata_id::Int
    frame::Int
end

immutable VehicleSet
    trajdata_id::Int
    roadway_name::Symbol
    descriptions::Dict{Int, ExportVehicleDescription} # id -> desc

    VehicleSet(trajdata_id::Int, roadway_name::Symbol) =
        new(trajdata_id, roadway_name, Dict{Int, ExportVehicleDescription}())
end

type DatasetExport
    vehiclesets::Dict{Int, VehicleSet} # contains base vehicle information, keyed by trajdata_id
    states::Vector{ExportState} # list of vehicle states (for each scene)
    scenes::Vector{ExportScene}

    DatasetExport() = new(Dict{Symbol, Dict{Int, ExportVehicleDescription}}(), ExportState[], ExportScene[])
end

function Base.append!(dsetexp::DatasetExport, dset::SceneDataset, trajdata_id::Int)

    state_ind = length(dsetexp.states)
    scene_ind = length(dsetexp.scenes)
    resize!(dsetexp.states, state_ind + nvehicles(dset))
    resize!(dsetexp.scenes, scene_ind + length(dset))

    for (scene, source) in zip(dset.scenes, dset.sources)

        lo = state_ind+1

        if !haskey(dsetexp.vehiclesets, trajdata_id)
            dsetexp.vehiclesets[trajdata_id] = VehicleSet(trajdata_id, scene.roadway_name)
        end

        vehicleset = dsetexp.vehiclesets[trajdata_id]
        @assert vehicleset.roadway_name == scene.roadway_name
        descriptions = vehicleset.descriptions

        for veh in scene
            if !haskey(descriptions, veh.id)
                descriptions[veh.id] = ExportVehicleDescription(veh.class, veh.length, veh.width)
            end

            let # DEBUG
                stored = descriptions[veh.id]
                @assert(stored.class == veh.class)
                @assert(stored.length == veh.length)
                @assert(stored.width == veh.width)
            end

            state_ind += 1
            dsetexp.states[state_ind] = ExportState(veh.id, veh.state)
        end

        hi = state_ind
        scene_ind += 1
        @assert(trajdata_id == source.trajdata_id)
        dsetexp.scenes[scene_ind] = ExportScene(lo, hi, trajdata_id, source.frame)
    end

    dsetexp
end

function reconstruct_vehicle(desc::ExportVehicleDescription, expstate::ExportState)
    veh = Vehicle()
    veh.id = expstate.id
    veh.class = desc.class
    veh.length = desc.length
    veh.width = desc.width
    veh.state = expstate.state
    veh
end
function reconstruct_scene(sceneexp::ExportScene, dsetexp::DatasetExport)

    vehicleset = dsetexp.vehiclesets[sceneexp.trajdata_id]
    descriptions = vehicleset.descriptions

    state_index_range = sceneexp.lo : sceneexp.hi
    vehicles = Array(Vehicle, length(state_index_range))
    for (a, state_index) in enumerate(state_index_range)
        expstate = dsetexp.states[state_index]
        vehicles[a] = reconstruct_vehicle(descriptions[expstate.id], expstate)
    end

    Scene(vehicleset.roadway_name, vehicles)
end
function reconstruct_dataset(dsetexp::DatasetExport, factors::Vector{Factor}, extractor::Scene=Scene())

    scenes = Array(Scene, length(dsetexp.scenes))
    sources = Array(SceneSource, length(scenes))
    structures = Array(SceneStructure, length(scenes))

    for (scene_ind, sceneexp) in enumerate(dsetexp.scenes)
        scene = reconstruct_scene(sceneexp, dsetexp)
        @assert(is_there_longitudinal_room(scene))
        source = SceneSource(sceneexp.trajdata_id, sceneexp.frame)
        structure = gen_scene_structure(scene, factors)::SceneStructure

        @assert(!isempty(structure.factor_assignments))
        scenes[scene_ind] = scene
        sources[scene_ind] = source
        structures[scene_ind] = structure
    end

    SceneDataset(scenes, sources, structures, factors)
end