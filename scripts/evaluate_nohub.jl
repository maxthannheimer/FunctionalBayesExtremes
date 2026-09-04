

date_string="2026_07_17_nohup"

folder = datadir("exp_raw",date_string)

files = sort(filter(f -> startswith(basename(f), "nohup_"),
                    readdir(folder, join=true)))

single_times = Float64[]
double_times = Float64[]

for file in files
    lines = readlines(file)

    mode = :none

    for line in lines
        if occursin("Time for single param", line)
            mode = :single
        elseif occursin("Time for double param", line)
            mode = :double
        else
            m = match(r"([\d\.]+)\s+seconds", line)
            if m !== nothing
                t = parse(Float64, m.captures[1])
                if mode == :single
                    push!(single_times, t)
                    mode = :none
                elseif mode == :double
                    push!(double_times, t)
                    mode = :none
                end
            end
        end
    end
end

println("Found $(length(single_times)) single-param runs.")
single_times_mean = mean(single_times)
println("Found $(length(double_times)) double-param runs.")
double_times_mean = mean(double_times)






# Times (seconds)
single_times = Float64[]
double_times = Float64[]

# Allocations (billions)
single_allocs = Float64[]
double_allocs = Float64[]

# Memory (TiB)
single_memory = Float64[]
double_memory = Float64[]

mode = :none

for file in files
    for line in eachline(file)

        if occursin("Time for single param", line)
            mode = :single
            continue
        elseif occursin("Time for double param", line)
            mode = :double
            continue
        end

        m = match(r"([\d\.]+)\s+seconds\s+\(([\d\.]+)\s+G allocations:\s+([\d\.]+)\s+TiB", line)

        if m !== nothing
            t = parse(Float64, m.captures[1])
            alloc = parse(Float64, m.captures[2])
            mem = parse(Float64, m.captures[3])

            if mode == :single
                push!(single_times, t)
                push!(single_allocs, alloc)
                push!(single_memory, mem)
            elseif mode == :double
                push!(double_times, t)
                push!(double_allocs, alloc)
                push!(double_memory, mem)
            end

            mode = :none
        end
    end
end

println("Found $(length(single_times)) single-parameter runs.")
println("Found $(length(double_times)) double-parameter runs.")

mean(single_times)
mean(double_times)

# Allocations (billions)
mean(single_allocs)
mean(double_allocs)

# Memory (TiB)
mean(single_memory)
mean(double_memory) 
