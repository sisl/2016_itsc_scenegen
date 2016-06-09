immutable DriveAction
    a::Float64
    ϕdes::Float64
end

#################

abstract FeatureSet

#################

function _draw_action_from_MvNormal!(a::Vector{Float64}, P::MvNormal)
    # this assumes it is 2x2

    μ::Vector{Float64} = P.μ
    Σ::Matrix{Float64} = P.Σ.mat

    a[1] = randn()*sqrt(Σ[1,1]) + μ[1]

    # compute conditional values
    μ₂ = μ[2] + Σ[1,2]*(a[1] - μ[1])/Σ[1,1]
    Σ₂ = Σ[2,2] - Σ[1,2]*Σ[1,2]/Σ[1,1]
    a[2] = randn()*sqrt(Σ₂) + μ₂

    a
end
function _get_e_and_denom(P::MvNormal, a::Vector{Float64})

    μ::Vector{Float64} = P.μ
    Σ::Matrix{Float64} = P.Σ.mat

    Δx = a[1] - μ[1]
    Δy = a[2] - μ[2]
    det = Σ[1,1]*Σ[2,2] - Σ[1,2]*Σ[1,2]

    e = exp(-0.5*((Σ[2,2]*Δx - Σ[1,2]*Δy)*Δx + (Σ[1,1]*Δy - Σ[1,2]*Δx)*Δy) / det)
    denom = sqrt(4*π*π*det)

    (e,denom)
end
function _pdf_of_MvNormal(P::MvNormal{PDMats.PDMat{Float64,Matrix{Float64}},Vector{Float64}}, a::Vector{Float64})
    # this assumes it is 2x2

    e, denom = _get_e_and_denom(P, a)
    e / denom
end
function _logpdf_of_MvNormal(P::MvNormal{PDMats.PDMat{Float64,Matrix{Float64}},Vector{Float64}}, a::Vector{Float64})
    # this assumes it is 2x2

    e, denom = _get_e_and_denom(P, a)
    log(e) - log(denom)
end

#################

abstract DriverModel

get_name(::DriverModel) = "???"

reset_hidden_state!(model::DriverModel) = model # do nothing by default
function prime_with_history!(model::DriverModel, trajdata::Trajdata, frame_start::Int, frame_end::Int, egoid::Int, scene::Scene=Scene())

    reset_hidden_state!(model)

    for frame in frame_start : frame_end
        get!(scene, trajdata, frame)
        observe!(model, scene, egoid)
    end

    model
end

#################

type StaticGaussian <: DriverModel
    P::MvNormal
end

get_name(::StaticGaussian) = "StaticGaussian"
observe!(model::StaticGaussian, scene::Scene, egoid::Int) = model # do nothing
function Base.rand(model::StaticGaussian)
    a = rand(model.P)
    DriveAction(a[1], a[2])
end
Distributions.pdf(model::StaticGaussian, a::DriveAction) = pdf(model.P, [a.a, a.ϕdes])
Distributions.logpdf(model::StaticGaussian, a::DriveAction) = logpdf(model.P, [a.a, a.ϕdes])

#################

type CleanDriver <: DriverModel
    normal_a::Normal
    normal_ϕdes::Normal
    CleanDriver() = new(Normal(0.0,1.0), Normal(0.0,1.0))
end

get_name(::CleanDriver) = "CleanDriver"
Base.copy(::CleanDriver) = CleanDriver()
function observe!(model::CleanDriver, scene::Scene, egoid::Int)
    ego_index = get_index_of_first_vehicle_with_id(scene, egoid)
    veh = scene[ego_index]
    t = veh.state.posF.t
    ϕdes = clamp(-0.04*t, -0.2, 0.2)
    model.normal_ϕdes = Normal(ϕdes, 0.01)

    a = 0.0
    ind_fore = get_neighbor_index_fore(scene, ego_index)
    if ind_fore > 0
        veh_fore = scene[ind_fore]
        distance_headway = get_headway_dist_between(veh, veh_fore)
        @assert(!isnan(distance_headway))
        relative_speed = get_vel_s(veh_fore.state) - get_vel_s(veh.state)
        @assert(!isnan(relative_speed))
        a = (distance_headway - 30.0)*0.2 + relative_speed*0.4
    end
    a = clamp(a, -15.0, 10.0)

    model.normal_a = Normal(a, 0.1)

    model
end
Base.rand(model::CleanDriver) = DriveAction(rand(model.normal_a), rand(model.normal_ϕdes))
Distributions.pdf(model::CleanDriver, a::DriveAction) = pdf(model.normal_a, a.a) * pdf(model.normal_ϕdes, a.ϕdes)
Distributions.logpdf(model::CleanDriver, a::DriveAction) = logpdf(model.normal_a, a.a) + logpdf(model.normal_ϕdes, a.ϕdes)