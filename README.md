# FunctionalBayesExtremes

This code base is using the [Julia Language](https://julialang.org/) to make a reproducible scientific project named
> FunctionalBayesExtremes

It is authored by Max Thannheimer.

To (locally) reproduce this project, do the following:

0. Download this code base.
1. Open Terminal and clone repository:
   ```
   git clone https://github.com/maxthannheimer/FuncionalBayesExtremes.git
   cd FuncionalBayesExtremes/biometrika-code
   julia --project=.
   ```
2. In Julia
   ```
   >julia using Pkg
   >julia Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.

In the biometrika-code folder you can find the code necessary for reproducing our simulation study. One can just run the 'run_MCMC_server.jl' to run the same simulations as on the study, without parallelization, or directly the 'start_julia_instances.sh' to run the code in parallel on a linux machine.


