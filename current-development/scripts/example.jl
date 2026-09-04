using DrWatson
@quickactivate "FunctionalBayesExtremes"

projectname()
projectdir()
datadir()


""" Generate some data and save it to the data directory. """

""" Prepare Simulation
In Julia, @. is the broadcasting macro. 
It applies operations element-wise to arrays or collections within the expression that follows it. 
For example, in a function like fakesim, if you have something like @. x + y, 
it computes x .+ y (element-wise addition) without needing to write the dot explicitly for each operator.

"""
function fakesim(a,b,v,method="linear")
    if method == "linear"
        return r = @. a+ b*v
    elseif method == "cubic"
        return r = @. a+ b*v^3
    end
    y = sqrt(b)
    return r,y
end

a, b = 2, 3 
v =rand(5)
method = "linear"
r, y = fakesim(a,b,v,method)


params= @strdict a b v method

allparams = Dict(
    "a" => [1,2],
    "b" => [3,4],
    "v" => [rand(1:2, 2)],
    "method" => "linear",
)

dicts=dict_list(allparams)


""" Run and save the results of the simulation. """

function makesim(d::Dict)
    @unpack a, b, v, method = d
    r, y = fakesim(a,b,v,method)
    fulld= copy(d)
    fulld["r"] = r
    fulld["y"] = y
    return fulld
end

#manual run and save loop
for (i,d) in enumerate(dicts)
    f = makesim(d)
    wsave(datadir("simulations", "sim_$(i).jld2"), f)
end
safesave(datadir("ana", "linear.jld2"), @strdict analysis)


readdir(datadir("simulations"))

#savename is a helper function that generates a filename based on the parameters of the simulation by itself
savename(params)
savename(dicts[1])
savename(dicts[1],"jld2")

# smarter run and save loop using savename
for (i,d) in enumerate(dicts)
    f = makesim(d)
    wsave(datadir("simulations", savename(d,"jld2")), f)
end
readdir(datadir("simulations"))


# another run with use of tagsave



"""
@tagsave is a macro that not only saves the data but also tags it with metadata like the gitcommit and the line in the script where it is saved.
"""   

for (i,d) in enumerate(dicts)
    f = makesim(d)
    @tagsave(datadir("simulations", savename(d,"jld2")), f)
end 
readdir(datadir("simulations"))

firstsim = readdir(datadir("simulations"))[1]

wload(datadir("simulations", firstsim))


""" analyse results with the collect_results function. """

using DataFrames

df = collect_results(datadir("simulations"))



""" with safesave one cannot overwrite files, it saves with a counter if the file already exists. """

analysis = 42

safesave(datadir("ana", "linear.jld2"), @strdict analysis)

