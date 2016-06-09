using Distributions
using Discretizers

import Distributions: pdf, logpdf, fit_mle, rand

type PiecewiseUniform <: ContinuousUnivariateDistribution
    cat  :: Categorical # defines probability of picking a particular bin
    disc :: LinearDiscretizer # allows for quick discretization stuff
    pdfs :: Vector{Float64}

    function PiecewiseUniform(cat::Categorical, disc::LinearDiscretizer)
        @assert(ncategories(cat) == nlabels(disc))

        pdfs = deepcopy(cat.p)
        pdfs ./= binwidths(disc)

        new(cat, disc, pdfs)
    end
end

function Base.show(io::IO, d::PiecewiseUniform)
    show(io, d.cat)
    println(io, "")
    show(io, d.disc)
end

Base.extrema(d::PiecewiseUniform) = extrema(d.disc)

#### Evaluation

function pdf(d::PiecewiseUniform, x::Float64)
    if x < d.disc.binedges[1] || x > d.disc.binedges[end]
        return 0.0
    end
    d.pdfs[encode(d.disc, x)]
end
function logpdf(d::PiecewiseUniform, x::Float64)
    if x < d.disc.binedges[1] || x > d.disc.binedges[end]
        return -Inf
    end
    log(d.pdfs[encode(d.disc, x)])
end

#### Sampling

rand(d::PiecewiseUniform) = decode(d.disc, rand(d.cat))

#### Fitting

function fit_mle{T<:Real}(::Type{PiecewiseUniform}, x_arr::AbstractArray{T}, nbins::Int;
    bounds::Tuple{Float64,Float64} = (NaN,NaN),
    laplace_smoothing_counts :: Int = 0
    )

    if isnan(bounds[1])
        binedges = binedges(DiscretizeUniformWidth(nbins), x_arr)
    else
        lo = bounds[1]
        hi = bounds[2]
        if hi < lo
            hi, lo = lo, hi
        end
        binedges = collect(linspace(lo, hi, nbins+1))
    end

    fit_mle(PiecewiseUniform, x_arr, binedges, laplace_smoothing_counts=laplace_smoothing_counts)
end
function fit_mle{T<:Real}(::Type{PiecewiseUniform}, x_arr::Vector{T}, binedges::Vector{Float64};
    laplace_smoothing_counts :: Int = 0
    )

    disc = LinearDiscretizer(binedges)
    nbins = length(binedges)-1

    counts = fill(laplace_smoothing_counts, nbins)
    for bin in encode(disc, x_arr)
        counts[bin] += 1
    end
    cat = Categorical(counts ./ sum(counts))
    PiecewiseUniform(cat, disc)
end