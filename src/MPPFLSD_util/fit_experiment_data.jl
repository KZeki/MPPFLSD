

function fit_experiment_data(
    df::DataFrame;
    flip_x::Bool=false,
    flip_y::Bool=false,
    shift_x::Symbol=:minimum_y,
    shift_y::Symbol=:volume,
    max_x_width::Union{Float64,Nothing}=nothing)

    if flip_x
        x_coords = -df.X;
    else
        x_coords = df.X;
    end

    if flip_y
        y_coords = -df.Z;
    else
        y_coords = df.Z;
    end
    
    if shift_x == :minimum_y
        x_zero = minimum(x_coords[y_coords .== minimum(y_coords)])
        mask = (minimum(x_zero) .<= x_coords )
        # shift x
        x_coords = x_coords .- x_zero

    elseif shift_x == :minimum_x
        x_coords = x_coords .- minimum(x_coords)
    else
        @warn "kwarg `shift_x::Symbol=:$align_y_on` of function `fit_experiment_data` is not recognized"
    end

    mask = (mask .&& ( max_x_width !== nothing ? x_coords .<= max_x_width : true))

    if shift_y == :volume
        # calculate shift
        y_shift = compute_vertical_shift_for_zero_on_volume(x_coords[mask].- minimum(x_coords[mask]),y_coords[mask]);
        # shift 
        y_coords = y_coords .+ y_shift;
    else
        @warn "kwarg `shift_y::Symbol=:$align_y_on` of function `fit_experiment_data` is not recognized"
    end
   
    return x_coords, y_coords, mask
end