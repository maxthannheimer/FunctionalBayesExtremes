""" Here different structs are defined, which are used in the rest of the code. """


""" Struct for the parameters of the model. """
struct Parameter
    α::Float64
    β::Float64
    c::Float64
    # inner constructor with keyword arguments only 
    function Parameter(; α::Float64, β::Float64, c::Float64)
        new(α, β, c)
    end
end
#testing Parameter struct
p=Parameter(α=1.0, β=2.0, c=3.0)


""" Struct for the Observation Simulation. """
struct Grid
    gridsize::Int
    N_coarse::Int
    coord_fine::Matrix{Float64}
    coord_coarse::Matrix{Float64}
    rows_coord_coarse::Vector{Int}
    x0::Int
    coord_x0::Vector{Float64}
    function Grid(;gridsize::Int, N_coarse::Int)
        coord_fine = ones(gridsize*gridsize,2)
        for x in 0:(gridsize-1) 
            for y in 0:(gridsize-1)
                coord_fine[y+x*gridsize+1,2] = x/(gridsize)
                coord_fine[y+x*gridsize+1,1] = y/(gridsize)
            end
        end
        rows_coord_coarse_and_x0 = sample(1:gridsize^2,N_coarse+1,replace=false)
        rows_coord_coarse = rows_coord_coarse_and_x0[2:end]
        coord_coarse = coord_fine[rows_coord_coarse,:]
        x0 = rows_coord_coarse_and_x0[1] # index of the normalization point
        coord_x0 = coord_fine[x0,:]
        new(gridsize, N_coarse, coord_fine, coord_coarse, rows_coord_coarse, x0,coord_x0)
    end
    function Grid()
        gridsize=9
        N_coarse=3
        coord_fine=ones(gridsize*gridsize,2)
        for x in 0:(gridsize-1) 
            for y in 0:(gridsize-1)
                coord_fine[y+x*gridsize+1,2]=x/(gridsize)
                coord_fine[y+x*gridsize+1,1]=y/(gridsize)
            end
        end
        rows_coord_coarse = [11,14,17,38,44,65,68,71]
        coord_coarse = coord_fine[rows_coord_coarse,:]
        x0 = 41 # index of the normalization point
        new(gridsize, N_coarse, coord_fine, coord_coarse, rows_coord_coarse, x0,coord_fine[x0,:])
    end
end


#Testing Grid struct
grid=Grid(gridsize=10, N_coarse=5)
grid=Grid()
grid.coord_fine


struct Observation
    sim_data::Matrix{Float64}
    obs_data::Matrix{Float64}
    obs_x0::Vector{Float64}
    function Observation(;sim_data::Matrix{Float64}, obs_data::Matrix{Float64}, obs_x0::Vector{Float64})
        new(sim_data, obs_data, obs_x0)
    end
    function Observation(;param::Parameter, grid::Grid, num_runs::Int, num_sim::Int)
        sim_data, obs_data, obs_x0 = simulate_pareto_process(param=param, grid=grid, num_runs=num_runs, num_sim=num_sim)
        new(sim_data, obs_data, obs_x0)
    end
end

struct Sum
    sumand::Float64
    a::Float64
    function Sum(;a,b)
         new(a+b,a)
    end
end
Sum(a=1.0, b=2.0).a

