#!/bin/julia
info("Received $(length(ARGS)) job(s):\n\t" * join(ARGS, "\n\t"))

for arg in ARGS
    if !endswith(arg, ".yml")
        warn("Skipping $arg, which does not end on '.yml'")
    elseif !isfile(arg)
        warn("Skipping $arg, which is not a file")
    else
        
        # run the job with the given config
        Job.run(arg)
        
    end
end
