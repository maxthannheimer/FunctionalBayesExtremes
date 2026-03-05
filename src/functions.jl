""" Here different functions are defined, which are used in the rest of the code. """

""" variogram of the process 
 c⋅||x||^β """
function vario(x::Vector{Float64},param::Parameter)::Float64
    return param.c*sqrt(
        (x[1])^2
        + (x[2])^2
        )^param.β
end

""" covariance function for two locations x and y 
 c⋅||x-x0||^β + c⋅||y-x0||^β -c⋅||x-y||^β """
function cov_fun_vario(param::Parameter,x::Vector{Float64},y::Vector{Float64},grid::Grid)::Float64
    return vario(x-grid.coord_x0,param) + vario(y-grid.coord_x0,param)-vario(x-y,param)
end

