
type IntelligentDriverModel <: DriverModel
    P::MvNormal # static covariance on top of mean accel based on IDM, mean turnrate a proportional controller to track lane center

    k_spd::Float64 # proportional constant for speed tracking [s⁻¹]
    k_lat::Float64 # proportional constant for lane tracking [1/ft⋅s]

    δ::Float64 # acceleration exponent [-]
    T::Float64 # desired time headway [s]
    v_des::Float64 # desired speed [ft/s]
    s_min::Float64 # minimum acceptable gap [ft]
    a_max::Float64 # maximum acceleration ability [ft/s²]
    d_cmf::Float64 # comfortable deceleration [ft/s²]

    function IntelligentDriverModel(;
        σ₁::Float64    =     0.1,
        σ₂::Float64    =     0.001,
        k_spd::Float64 =   1.0,
        k_lat::Float64 =     0.05,
        δ::Float64     =   4.0,
        T::Float64     =   1.5,
        v_des::Float64 =  10.0, # typically overwritten
        s_min::Float64 =   6.0,
        a_max::Float64 =  10.0,
        d_cmf::Float64 =  10.0,
        )

        retval = new()
        retval.P = MvNormal([0.0,0.0], [σ₁^2 0.0; 0.0 σ₂^2])
        retval.k_spd = k_spd
        retval.k_lat = k_lat
        retval.δ     = δ
        retval.T     = T
        retval.v_des = v_des
        retval.s_min = s_min
        retval.a_max = a_max
        retval.d_cmf = d_cmf
        retval
    end
end

get_name(::IntelligentDriverModel) = "IDM"
function Base.copy(model::IntelligentDriverModel)
    IntelligentDriverModel(σ₁ = sqrt(model.P.Σ.mat[1,1]),
                           σ₂ = sqrt(model.P.Σ.mat[2,2]),
                           k_spd = model.k_spd,
                           k_lat = model.k_lat,
                           δ = model.δ,
                           T = model.T,
                           v_des = model.v_des,
                           s_min = model.s_min,
                           a_max = model.a_max,
                           d_cmf = model.d_cmf)
end
function observe!(model::IntelligentDriverModel, scene::Scene, carid::Int)

    # update the mean of P based on IDM prediction

    ego_index = get_index_of_first_vehicle_with_id(scene, carid)
    veh_ego = scene[ego_index]
    v = veh_ego.state.v

    ind_fore = get_neighbor_index_fore(scene, ego_index)
    if ind_fore > 0
        veh_fore = scene[ind_fore]
        s_gap = get_headway_dist_between(veh_ego, veh_fore)

        Δv = veh_fore.state.v - v
        s_des = model.s_min + v*model.T - v*Δv / (2*sqrt(model.a_max*model.d_cmf))
        a_idm = model.a_max * (1.0 - (v/model.v_des)^model.δ - (s_des/s_gap)^2)
        if isnan(a_idm)
            println("v:     ", v)
            println("Δv:    ", Δv)
            println("v_des: ", v_des)
            println("s_des: ", s_des)
            println("s_gap: ", s_gap)
            println("v_for: ", veh_fore.state.v)
            println("a_idm: ", a_idm)
        end
        @assert(!isnan(a_idm))

        model.P.μ[1] = a_idm
    else
        # no lead vehicle, just drive to match desired speed
        Δv = model.v_des - v
        model.P.μ[1] = Δv*model.k_spd # predicted accel to match target speed
    end

    model.P.μ[1] = clamp(model.P.μ[1], -model.d_cmf, model.a_max)
    model.P.μ[2] = -veh_ego.state.posF.t * model.k_lat # track lane centerline

    model
end
function Base.rand(model::IntelligentDriverModel)
    a = [0.0,0.0]
    _draw_action_from_MvNormal!(a, model.P)
    DriveAction(a[1], a[2])
    # DriveAction(model.P.μ[1], model.P.μ[2])
end
get_action_logl(model::IntelligentDriverModel, a::DriveAction) = _logpdf_of_MvNormal(model.P, [a.a, a.ϕdes])
get_action_likelihood(model::IntelligentDriverModel, a::DriveAction) = _pdf_of_MvNormal(model.P, [a.a, a.ϕdes])