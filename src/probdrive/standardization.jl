immutable Standardization
    μ::Float64 # mean
    σ::Float64 # standard deviation
end
function calc_standardization(arr::Vector{Float64})
    μ = mean(arr)
    σ = stdm(arr, μ)
    Standardization(μ, σ)
end
standardize(v::Real, s::Standardization) = (v - s.μ)/s.σ
function standardize!{R<:Real}(arr::Vector{R}, s::Standardization)
    for i in 1 : length(arr)
        arr[i] = (arr[i] - s.μ) / s.σ
    end
    arr
end
destandardize(v::Real, s::Standardization) = v*s.σ + s.μ
destandardize_var(v::Real, s::Standardization) = v*s.σ*s.σ
function destandardize!{R<:Real}(arr::Vector{R}, s::Standardization)
    for i in 1 : length(arr)
        arr[i] = arr[i]*s.σ + s.μ
    end
    arr
end