function _cycle_linear_plot!(xarr::Vector{Float64}, yarr::Vector{Float64})
    append!(yarr, [0.0, 0.0, yarr[1]])
    append!(xarr, [xarr[end], xarr[1], xarr[1]])
    nothing
end
function plot_marginal_roadfactor(f_road::Factor;
    nbins_v::Int = 101,
    nbins_t::Int = 101,
    nbins_ϕ::Int = 101,
    roadway::Roadway = NGSIM.ROADWAY_80,
    )

    pts_v = collect(linspace(BOUNDS_V_LO, BOUNDS_V_HI, nbins_v))
    pts_t = collect(linspace(BOUNDS_T_LO, BOUNDS_T_HI, nbins_t))
    pts_ϕ = collect(linspace(BOUNDS_ϕ_LO, BOUNDS_ϕ_HI, nbins_ϕ))
    v_values = zeros(Float64, nbins_v)
    t_values = zeros(Float64, nbins_t)
    ϕ_values = zeros(Float64, nbins_ϕ)
    tot = 0.0

    for (i,v) in enumerate(pts_v)
        set_variable!(f_road, IND_V, v)
        for (j,t) in enumerate(pts_t)
            set_variable!(f_road, IND_T, t)
            for (k,ϕ) in enumerate(pts_ϕ)
                set_variable!(f_road, IND_ϕ, ϕ)
                extract_features!(f_road, roadway)
                e = eval_factor(f_road)
                v_values[i] += e
                t_values[j] += e
                ϕ_values[k] += e
                tot += e
            end
        end
    end

    # produce a probability density
    v_values ./= tot
    t_values ./= tot
    ϕ_values ./= tot

    plot_width = "20cm"
    plot_height = "4cm"
    style = "solid, thick, mark=none, black"

    _cycle_linear_plot!(pts_v, v_values)
    _cycle_linear_plot!(pts_t, t_values)
    _cycle_linear_plot!(pts_ϕ, ϕ_values)

    g = GroupPlot(1, 3, groupStyle = "vertical sep = 1cm")
    push!(g, Axis(PGFPlots.Plots.Linear(pts_v, v_values, style=style*", fill=monokai3"), xlabel="speed [ft/s]", ylabel="pdf", width=plot_width, height=plot_height, enlargelimits="false"))
    push!(g, Axis(PGFPlots.Plots.Linear(pts_t, t_values, style=style*", fill=monokai4"), xlabel="lane centerline offset [ft]", ylabel="pdf", width=plot_width, height=plot_height, enlargelimits="false"))
    push!(g, Axis(PGFPlots.Plots.Linear(rad2deg(pts_ϕ), ϕ_values, style=style*", fill=monokai5"), xlabel="lane-relative heading [deg]", ylabel="pdf", width=plot_width, height=plot_height, enlargelimits="false"))
    g
end
function plot_marginal_roadfactor_empirical(dset::SceneDataset;
    nbins_v::Int = 199,
    nbins_t::Int = 99,
    nbins_ϕ::Int = 99,
    )

    min_v, max_v = Inf, -Inf
    min_t, max_t = Inf, -Inf
    min_ϕ, max_ϕ = Inf, -Inf

    disc_v = LinearDiscretizer(collect(linspace(BOUNDS_V_LO, BOUNDS_V_HI, nbins_v+1)))
    disc_t = LinearDiscretizer(collect(linspace(BOUNDS_T_LO, BOUNDS_T_HI, nbins_t+1)))
    disc_ϕ = LinearDiscretizer(collect(linspace(BOUNDS_ϕ_LO, BOUNDS_ϕ_HI, nbins_ϕ+1)))
    v_values = zeros(Float64, nbins_v)
    t_values = zeros(Float64, nbins_t)
    ϕ_values = zeros(Float64, nbins_ϕ)
    tot = 0.0

    m = length(dset)
    for i in 1 : m
        scene, structure = get_scene_and_structure(dset, i)

        for vehicle_index in structure.active_vehicles

            veh = scene[vehicle_index]

            v = veh.state.v
            t = veh.state.posF.t
            ϕ = veh.state.posF.ϕ

            min_v, max_v = min(min_v, v), max(max_v, v)
            min_t, max_t = min(min_t, t), max(max_t, t)
            min_ϕ, max_ϕ = min(min_ϕ, ϕ), max(max_ϕ, ϕ)

            v_values[encode(disc_v, v)] += 1.0
            t_values[encode(disc_t, t)] += 1.0
            ϕ_values[encode(disc_ϕ, ϕ)] += 1.0
        end
    end

    # produce a probability density
    v_values ./= m
    t_values ./= m
    ϕ_values ./= m

    plot_width = "20cm"
    plot_height = "4cm"
    style = "solid, thick, mark=none, black"

    pts_v = bincenters(disc_v)
    pts_t = bincenters(disc_t)
    pts_ϕ = bincenters(disc_ϕ)

    _cycle_linear_plot!(pts_v, v_values)
    _cycle_linear_plot!(pts_t, t_values)
    _cycle_linear_plot!(pts_ϕ, ϕ_values)

    g = GroupPlot(1, 3, groupStyle = "vertical sep = 1cm")
    push!(g, Axis([PGFPlots.Plots.Linear(pts_v, v_values, style=style*", fill=monokai3"),
                   PGFPlots.Plots.Node(@sprintf("[%.3f, %.3f]", min_v, max_v),0.8,0.8,axis="axis description cs"),
        ], xlabel="speed [ft/s]", ylabel="pdf", width=plot_width, height=plot_height, enlargelimits="false"))
    push!(g, Axis([PGFPlots.Plots.Linear(pts_t, t_values, style=style*", fill=monokai4"),
                   PGFPlots.Plots.Node(@sprintf("[%.3f, %.3f]", min_t, max_t),0.8,0.8, axis="axis description cs"),
        ], xlabel="lane centerline offset [ft]", ylabel="pdf", width=plot_width, height=plot_height, enlargelimits="false"))
    push!(g, Axis([PGFPlots.Plots.Linear(rad2deg(pts_ϕ), ϕ_values, style=style*", fill=monokai5"),
                   PGFPlots.Plots.Node(@sprintf("[%.3f, %.3f]", min_ϕ, max_ϕ),0.8,0.8,axis="axis description cs"),
        ], xlabel="lane-relative heading [deg]", ylabel="pdf", width=plot_width, height=plot_height, enlargelimits="false"))
    g
end
function plot_marginal_roadfactor_empirical(scenes::Vector{Scene}, factors::Vector{Factor};
    nbins_v::Int = 199,
    nbins_t::Int = 99,
    nbins_ϕ::Int = 99,
    )

    min_v, max_v = Inf, -Inf
    min_t, max_t = Inf, -Inf
    min_ϕ, max_ϕ = Inf, -Inf

    disc_v = LinearDiscretizer(collect(linspace(BOUNDS_V_LO, BOUNDS_V_HI, nbins_v+1)))
    disc_t = LinearDiscretizer(collect(linspace(BOUNDS_T_LO, BOUNDS_T_HI, nbins_t+1)))
    disc_ϕ = LinearDiscretizer(collect(linspace(BOUNDS_ϕ_LO, BOUNDS_ϕ_HI, nbins_ϕ+1)))
    v_values = zeros(Float64, nbins_v)
    t_values = zeros(Float64, nbins_t)
    ϕ_values = zeros(Float64, nbins_ϕ)
    tot = 0.0

    m = length(dset)
    for scene in scenes
        structure = gen_scene_structure(scene, factors)

        for vehicle_index in structure.active_vehicles

            veh = scene[vehicle_index]

            v = veh.state.v
            t = veh.state.posF.t
            ϕ = veh.state.posF.ϕ

            min_v, max_v = min(min_v, v), max(max_v, v)
            min_t, max_t = min(min_t, t), max(max_t, t)
            min_ϕ, max_ϕ = min(min_ϕ, ϕ), max(max_ϕ, ϕ)

            v_values[encode(disc_v, v)] += 1.0
            t_values[encode(disc_t, t)] += 1.0
            ϕ_values[encode(disc_ϕ, ϕ)] += 1.0
        end
    end

    # produce a probability density
    v_values ./= m
    t_values ./= m
    ϕ_values ./= m

    plot_width = "20cm"
    plot_height = "4cm"
    style = "solid, thick, mark=none, black"

    pts_v = bincenters(disc_v)
    pts_t = bincenters(disc_t)
    pts_ϕ = bincenters(disc_ϕ)

    _cycle_linear_plot!(pts_v, v_values)
    _cycle_linear_plot!(pts_t, t_values)
    _cycle_linear_plot!(pts_ϕ, ϕ_values)

    g = GroupPlot(1, 3, groupStyle = "vertical sep = 1cm")
    push!(g, Axis([PGFPlots.Plots.Linear(pts_v, v_values, style=style*", fill=monokai3"),
                   PGFPlots.Plots.Node(@sprintf("[%.3f, %.3f]", min_v, max_v),0.8,0.8,axis="axis description cs"),
        ], xlabel="speed [ft/s]", ylabel="pdf", width=plot_width, height=plot_height, enlargelimits="false"))
    push!(g, Axis([PGFPlots.Plots.Linear(pts_t, t_values, style=style*", fill=monokai4"),
                   PGFPlots.Plots.Node(@sprintf("[%.3f, %.3f]", min_t, max_t),0.8,0.8, axis="axis description cs"),
        ], xlabel="lane centerline offset [ft]", ylabel="pdf", width=plot_width, height=plot_height, enlargelimits="false"))
    push!(g, Axis([PGFPlots.Plots.Linear(rad2deg(pts_ϕ), ϕ_values, style=style*", fill=monokai5"),
                   PGFPlots.Plots.Node(@sprintf("[%.3f, %.3f]", min_ϕ, max_ϕ),0.8,0.8,axis="axis description cs"),
        ], xlabel="lane-relative heading [deg]", ylabel="pdf", width=plot_width, height=plot_height, enlargelimits="false"))
    g
end