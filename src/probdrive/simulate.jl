const KP_DESIRED_ANGLE = 1.0

function propagate(veh::Vehicle, action::DriveAction, roadway::Roadway, Δt::Float64=NGSIM_TIMESTEP, n_integration_steps::Int=4)

    a = action.a # accel
    ϕdes = action.ϕdes # desired heading angle

    @assert(!isinf(a))
    @assert(!isnan(a))
    @assert(!isinf(ϕdes))
    @assert(!isnan(ϕdes))

    center = get_center(veh)

    x = center.x
    y = center.y
    θ = veh.state.posG.θ
    v = veh.state.v

    δt = Δt/n_integration_steps

    for i in 1 : n_integration_steps

        posG = VecSE2(x, y, θ) + polar(veh.length/2, θ)
        posF = project_posG_to_frenet(posG, roadway)
        ω = (ϕdes - posF.ϕ)*KP_DESIRED_ANGLE

        x += v*cos(θ)*δt
        y += v*sin(θ)*δt
        θ += ω*δt
        v += a*δt
    end

    posG = VecSE2(x, y, θ) + polar(veh.length/2, θ) # project from center to new car front
    posF = project_posG_to_frenet(posG, roadway)
    vehstate = VehicleState(posG, posF, v)
end

function get_actions!(
    actions::Vector{DriveAction},
    scene::Scene,
    models::Dict{Int, DriverModel}, # id → model
    )


    i = 0
    for veh in scene
        if haskey(models, veh.id)
            model = models[veh.id]
            observe!(model, scene, veh.id)
            actions[i+=1] = rand(model)
        end
    end

    actions
end

function tick!(
    scene::Scene,
    models::Dict{Int, DriverModel}, # id → model
    trajdata::Trajdata, # for propagating non-modeled cars
    frame::Int,
    actions::Vector{DriveAction},
    )

    roadway = get_roadway(scene)

    j = 0
    for veh in scene
        if haskey(models, veh.id)
            veh.state = propagate(veh, actions[j+=1], roadway)
        else
            if iscarinframe(trajdata, veh.id, frame)
                veh.state = get_vehiclestate(trajdata, veh.id, frame)
            else
                veh.state = VehicleState(VecSE2(NaN,NaN,NaN), NaN) # car disappears
            end
        end
    end

    scene
end
function tick!(
    scene::Scene,
    actions::Vector{DriveAction},
    )

    roadway = get_roadway(scene)

    for (veh_index, veh) in enumerate(scene)
        veh.state = propagate(veh, actions[veh_index], roadway)
    end

    scene
end