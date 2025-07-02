function compute_vertical_shift_for_zero_on_volume(x::Vector{Float64},y::Vector{Float64})::Float64
    return -trapz(x.^2,y)/(x[end]-x[begin]).^2
    #return -2trapz(x,y.*x)/(x[end]-x[begin]).^2
end