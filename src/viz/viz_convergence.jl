const EXPECTATION_SOLVER_COLOR_DICT = Dict(
            DirectSampling => "monokai3",
            DirectPropSampling => "monokai3",
            ImportanceSampling => "monokai5",
            ImportancePropSampling => "monokai5",
        )

function get_fill_region_arrays{R}(x_arr::Vector, y_lo::Vector{R}, y_hi::Vector{R})
    x_arr2 = [x_arr; reverse(x_arr); x_arr[1]]
    y_arr2 = [y_lo; reverse(y_hi); y_lo[1]]
    (x_arr2, y_arr2)
end

function get_convergence_data{T<:ExpectationSolver}(
    sol::T,
    arr_nsamples::AbstractVector{Int} = round(Int, logspace(2.0, 4.0, 3)), # 10
    )

    arr_μ = Array(Float64, length(arr_nsamples))
    arr_σ = Array(Float64, length(arr_nsamples)) # variance of the mean, σ²/n

    sol2 = copy(sol)

    for (i,nsamples) in enumerate(arr_nsamples)

        solve!(sol2, nsamples - get_nsamples(sol2))

        arr_μ[i] = mean(sol2)
        arr_σ[i] = sqrt(var(sol2)/nsamples)
    end

    (arr_μ, arr_σ)
end
function get_convergence_data{T<:ExpectationSolver}(
    sol::T,
    n_to_average::Int,
    arr_nsamples::AbstractVector{Int} = round(Int, logspace(2.0, 4.0, 3)), # 10
    )

    arr_updaters = Array(OnlineExtrema, length(arr_nsamples))
    for i in 1 : length(arr_updaters)
        arr_updaters[i] = OnlineExtrema()
    end

    for i in 1 : n_to_average
        sol2 = copy(sol)
        for (j,nsamples) in enumerate(arr_nsamples)
            solve!(sol2, nsamples - get_nsamples(sol2))
            update!(arr_updaters[j], mean(sol2))
        end
    end

    arr_μ = Array(Float64, length(arr_updaters))
    arr_σ = Array(Float64, length(arr_μ))
    for (i,v) in enumerate(arr_updaters)
        arr_μ[i] = v.μ
        arr_σ[i] = sqrt( v.M/( v.n - 1))
    end

    (arr_μ, arr_σ)
end
function get_convergence_plots(
    sol::ExpectationSolver,
    f_true::Float64;
    arr_nsamples::Vector{Int} = round(Int, logspace(2.0, 4.0, 10)),
    color_dict::Dict{DataType, ASCIIString} = EXPECTATION_SOLVER_COLOR_DICT,
    )

    n_to_average = 5
    (arr_μ, arr_σ) = get_convergence_data(sol, n_to_average, arr_nsamples)

    arr_hi = arr_μ + arr_σ
    arr_lo = arr_μ - arr_σ

    rel_errors = abs((arr_μ ./ -f_true) .+ 1.0)

    rel_errors_hi = Array(Float64, length(rel_errors))
    for i in 1 : length(rel_errors)
        μ, σ = arr_μ[i], arr_σ[i]
        ε₁ = abs(((μ + σ) / -f_true) + 1.0)
        ε₂ = abs(((μ - σ) / -f_true) + 1.0)
        rel_errors_hi[i] = max(ε₁, ε₂, rel_errors[i])
    end


    color = color_dict[typeof(sol)]

    p_estimates = [PGFPlots.Plots.Linear(get_fill_region_arrays(arr_nsamples, arr_lo, arr_hi)...,
                    style="draw="*color*", mark=none, fill="*color*", fill opacity=0.3, forget plot"),
                   PGFPlots.Plots.Linear(arr_nsamples, arr_μ,
                    style="solid, thick, "*color*", mark=none", legendentry=string(typeof(sol)))]

    p_rel_errors = [PGFPlots.Plots.Linear(get_fill_region_arrays(arr_nsamples, rel_errors, rel_errors_hi)...,
                    style="draw=none, mark=none, fill="*color*", fill opacity=0.3, forget plot"),
                   PGFPlots.Plots.Linear(arr_nsamples, rel_errors,
                    style="solid, ultra thick, "*color*", mark=none, forget plot")]

    (p_estimates, p_rel_errors)
end

function plot_convergence{E<:ExpectationSolver}(
    solvers::Vector{E},
    f_true::Float64;
    arr_nsamples::Vector{Int} = round(Int, logspace(1.0, 4.0, 101)),
    color_dict::Dict{DataType, ASCIIString} = EXPECTATION_SOLVER_COLOR_DICT,
    )

    p_estimates = PGFPlots.Plots.Plot[]
    p_rel_errors = PGFPlots.Plots.Plot[]

    for sol in solvers
        a, b = get_convergence_plots(sol, f_true,
                                     arr_nsamples=arr_nsamples, color_dict=color_dict)
        append!(p_estimates, a)
        append!(p_rel_errors, b)
    end

    g = GroupPlot(2, 1, groupStyle="horizontal sep=2.0cm")
    push!(g, Axis(p_estimates, xlabel="{number of samples}", ylabel="{probability of adverse event}",
                  width="10cm", xmode="log", style="log ticks with fixed point"))

    push!(g, Axis(p_rel_errors, xlabel="{number of samples}", ylabel="{relative error}",
                  width="10cm", xmode="log", ymode="log", style="log ticks with fixed point"))
    g
end