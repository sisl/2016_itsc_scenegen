immutable TargetsAccelDesAng
    frames_per_tick::Int # number of NGSIM frames per tick
end
ntargets(::TargetsAccelDesAng) = 2
required_horizon(ttype::TargetsAccelDesAng) = ttype.frames_per_tick

function extract_targets!(ttype::TargetsAccelDesAng,
    targets::Vector{Float64},
    trajdata::Trajdata,
    frame::Int,
    egoid::Int,
    )

    s₀ = get_vehiclestate(trajdata, egoid, frame)
    v₀ = s₀.v
    ϕ₀ = s₀.posF.ϕ

    s₁ = get_vehiclestate(trajdata, egoid, frame+ttype.frames_per_tick)
    v₁ = s₁.v
    ϕ₁ = s₁.posF.ϕ

    KP_DESIRED_ANGLE = 1.0
    ΔT = NGSIM_TIMESTEP*ttype.frames_per_tick
    expconst = exp(-KP_DESIRED_ANGLE*ΔT)

    a = (v₁ - v₀) / ΔT
    ϕdes = (ϕ₁ - ϕ₀*expconst) / (1.0 - expconst)

    targets[1] = a
    targets[2] = ϕdes
end