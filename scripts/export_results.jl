# jl scripts/export_results.jl
include("../full_scene.jl")

push!(LOAD_PATH, "/home/tim/Documents/wheelerworkspace/UtilityCode/")
using LaTeXeXport

const TEXFILE = "/home/tim/Documents/papers/2016_itsc_scenegen_wheeler/2016_itsc_scenegen.tex"
const TEXDIR = splitdir(TEXFILE)[1]
typealias DataSets Dict{Symbol, SceneDataset}

function load_datasets()
    retval = Dict{Symbol, SceneDataset}()

    dsetexp = JLD.load("../output/dsetexport.jld", "dsetexp")
    dset_real = reconstruct_dataset(dsetexp, create_core_factors())

    retval[:real] = dset_real
    retval[:FG] = JLD.load(joinpath(SAVE_DIR, "dset_FG.jld"), "dset")
    retval[:BN] = JLD.load(joinpath(SAVE_DIR, "dset_BN.jld"), "dset")
    retval
end

type DistrMetric
    label::AbstractString # ex: "speed [\\si{m/s}]"
    disc::AbstractDiscretizer
    extract::Function # (val::Float64, worked::Bool) ← metric.extract(scene, structure, vehicle_index)

    function DistrMetric(label::AbstractString, lo::Float64, hi::Float64, nbins::Int, extract::Function)
        disc = LinearDiscretizer(collect(linspace(lo, hi, nbins)))
        new(label, disc, extract)
    end
end

extract_speed(scene::Scene, structure::SceneStructure, vehicle_index::Int) = (scene[vehicle_index].state.v, true)
extract_dcl(scene::Scene, structure::SceneStructure, vehicle_index::Int) = (scene[vehicle_index].state.posF.t, true)
extract_phi(scene::Scene, structure::SceneStructure, vehicle_index::Int) = (scene[vehicle_index].state.posF.ϕ, true)
function extract_timegap(scene::Scene, structure::SceneStructure, vehicle_index::Int)
    ind_fore = get_neighbor_index_fore(scene, vehicle_index)
    if ind_fore != 0
        timegap = get_headway_time_between(scene[vehicle_index], scene[ind_fore])
        @assert(!isnan(timegap))
        (timegap, true)
    else
        (NaN, false)
    end
end
function extract_headway(scene::Scene, structure::SceneStructure, vehicle_index::Int)
    ind_fore = get_neighbor_index_fore(scene, vehicle_index)
    if ind_fore != 0
        dist = get_headway_dist_between(scene[vehicle_index], scene[ind_fore])
        @assert(!isnan(dist))
        (dist, true)
    else
        (NaN, false)
    end
end
function extract_a_req(scene::Scene, structure::SceneStructure, vehicle_index::Int)
    ind_fore = get_neighbor_index_fore(scene, vehicle_index)
    if ind_fore != 0
        veh_rear = scene[vehicle_index]
        veh_fore = scene[ind_fore]

        sr = veh_rear.state.posF.s
        vr = veh_rear.state.v
        sf = veh_fore.state.posF.s
        vf = veh_fore.state.v

        dv = vf - vr
        if dv ≥ 0.0 # they are pulling away; we are good
            (0.0,true)
        end

        dx = sf - sr
        a_req = -dv*dv / (2dx)
        (a_req, true)
    else
        (NaN, false)
    end
end
function extract_ittc(scene::Scene, structure::SceneStructure, vehicle_index::Int)
    ind_fore = get_neighbor_index_fore(scene, vehicle_index)
    if ind_fore != 0
        veh_rear = scene[vehicle_index]
        veh_fore = scene[ind_fore]

        sr = veh_rear.state.posF.s
        vr = veh_rear.state.v
        sf = veh_fore.state.posF.s
        vf = veh_fore.state.v

        dv = vf - vr
        if dv ≥ 0.0 # they are pulling away; we are good
            (0.0,true)
        end

        dx = sf - sr
        ittc = -dv/dx
        (ittc, true)
    else
        (NaN, false)
    end
end

function get_counts(metric::DistrMetric, dset::SceneDataset)

    disc = metric.disc
    retval = zeros(Int, nlabels(disc))

    for i in 1 : length(dset)
        scene, structure = get_scene_and_structure(dset, i)
        for vehicle_index in structure.active_vehicles
            val, worked = metric.extract(scene, structure, vehicle_index)::Tuple{Float64, Bool}
            if worked
                retval[encode(disc, val)] += 1
            end
        end
    end

    retval
end

function KLdiv(countsP::Vector{Int}, countsQ::Vector{Int}, disc::LinearDiscretizer)
    totP = sum(countsP)
    totQ = sum(countsQ)

    retval = 0.0
    for bin in 1 : length(countsP)
        w = binwidth(disc, bin)
        P = countsP[bin] / totP
        Q = countsQ[bin] / totQ
        if P > 0.0
            retval += P * log(P/Q)
        end
    end
    retval
end

function _export_next_groupplot_line(io::IO, metric::DistrMetric, include_legend::Bool)
    println(io, "\\nextgroupplot[")
    # println(io, "\tylabel={$(metric.label)},")
    println(io, "\twidth=7.75cm, height=3cm,")
    println(io, "\tenlarge x limits=0.0,")
    println(io, "\tymin=0.0,")
    # println(io, "\tx label style={at={(axis description cs:0.5,-0.15)},anchor=north},")
    println(io, "\ty tick label style={")
    println(io, "\t\t/pgf/number format/.cd,")
    println(io, "\t\tfixed,")
    println(io, "\t\tfixed zerofill,")
    println(io, "\t\tprecision=2,")
    println(io, "\t/tikz/.cd")
    println(io, "},")
    if include_legend
        println(io, "\tlegend style={")
        println(io, "\t\tdraw=none,")
        println(io, "\t\tat={(0.5,-0.35)},")
        println(io, "\t\tanchor=north,")
        println(io, "\t\tlegend columns=3,")
        println(io, "\t\tfont=\\scriptsize,")
        println(io, "\t},")
    end
    println(io, "]")
end
function _export_next_groupplot_kldiv(io::IO)
    println(io, "\\nextgroupplot[")
    println(io, "\twidth=3cm, height=3cm,")
    println(io, "\txbar, xmin=0,")
    println(io, "\taxis x line*=bottom,")
    println(io, "\tenlarge y limits=6.0,")
    println(io, "\tvisualization depends on=x \\as \\rawx,")
    println(io, "\tsymbolic y coords={0,1},")
    println(io, "\tytick=\\empty,")
    println(io, "\tnodes near coords, nodes near coords align={horizontal},")
    println(io, "\tevery node near coord/.append style={")
    println(io, "\t\tfont=\\tiny,")
    println(io, "\t\tshift={(axis direction cs:-\\rawx,0)}")
    println(io, "\t},")
    println(io, "\taxis line style={white},")
    println(io, "% \txtick pos = left")
    println(io, "\txtick=\\empty,")
    println(io, "]")
end
function _write_line_plot(io::IO, options::AbstractString, counts::Vector{Int}, disc::LinearDiscretizer)
    println(io, "\\addplot[$(options)] coordinates{")
    tot_count = sum(counts)
    for (bin,count) in enumerate(counts)
        lo, hi = extrema(disc, bin)
        val = count/tot_count
        @printf(io, "(%.3f,%.3f) (%.3f,%.3f) ", lo, val, hi, val)
    end
    println(io, "};")
end
function export_distribution_compare(io::IO, dsets::DataSets, metrics::Vector{DistrMetric})

    println(io, "\\begin{groupplot}[")
    println(io, "\tgroup style={group size= 2 by $(length(metrics)), vertical sep=1.5em, horizontal sep=0.25em},")
    println(io, "\tticklabel style = {font=\\scriptsize},")
    println(io, "\ttitle style = {font=\\scriptsize, yshift=-1.5ex},")
    println(io, "\tlabel style = {font=\\scriptsize},")
    println(io, "]\n")

    for (i,metric) in enumerate(metrics)
        _export_next_groupplot_line(io, metric, i==length(metrics))

        binedges = metric.disc.binedges
        counts_real = get_counts(metric, dsets[:real])
        counts_FG = get_counts(metric, dsets[:FG]) .+ 1
        counts_BN = get_counts(metric, dsets[:BN]) .+ 1

        _write_line_plot(io, "solid, ultra thick, mark=none, colorA", counts_real, metric.disc)
        _write_line_plot(io, "solid, thick, mark=none, colorB", counts_BN, metric.disc)
        _write_line_plot(io, "solid, mark=none, black", counts_FG, metric.disc)
        println(io, "\\node at (axis description cs:0.95,0.95) [anchor=north east] {\\scriptsize $(metric.label)};")

        if i == length(metrics)
            println(io, "\\legend{True, Bayesian Network, Factor Graph}")
        end

        _export_next_groupplot_kldiv(io)

        kldiv_FG = KLdiv(counts_real, counts_FG, metric.disc)
        kldiv_BN = KLdiv(counts_real, counts_BN, metric.disc)

        println(io, "\\addplot+[xbar, draw=black, fill=black!40, text=black] coordinates {(", @sprintf("%.3f", kldiv_FG), ",0)};")
        println(io, "\\addplot+[xbar, draw=colorB,  fill=colorB!40, text=colorB] coordinates {(", @sprintf("%.3f", kldiv_BN), ",1)};")
    end

    println(io, "\\end{groupplot}")
end

dsets = load_datasets()

metrics = DistrMetric[]
push!(metrics, DistrMetric("speed [\\si{m/s}]",            BOUNDS_V_LO, BOUNDS_V_HI, 31, extract_speed))
push!(metrics, DistrMetric("lane offset [\\si{m}]",       -7.5,         7.5,         31, extract_dcl)) # BOUNDS_T_LO, BOUNDS_T_HI, 31, extract_dcl))
push!(metrics, DistrMetric("heading [\\si{rad}]",         -0.3,         0.3,         31, extract_phi)) # BOUNDS_ϕ_LO, BOUNDS_ϕ_HI, 31, extract_phi))
push!(metrics, DistrMetric("timegap [\\si{s}]",            0.0,         10.0,        31, extract_timegap))
push!(metrics, DistrMetric("headway [\\si{ft}]",           0.0,         200.0,       31, extract_headway))
push!(metrics, DistrMetric("inverse time to collision [\\si{1/s}]",             0.0,         0.6,         31, extract_ittc))

fh = STDOUT

write_to_texthook(TEXFILE, "distribution-compare") do fh
    export_distribution_compare(fh, dsets, metrics)
end

println("DONE EXPORTING RESULTS TO TEX")