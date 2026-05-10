# script to get the notebook organized
using Printf
Base.show(io::IO, f::Float64) = @printf(io, "%1.5e", f)
Base.show(io::IO, f::Float32) = @printf(io, "%1.5e", f)
Base.show(io::IO, f::Float16) = @printf(io, "%1.5e", f)
push!(LOAD_PATH, "./src")
ENV["JULIA_CONDAPKG_VERBOSITY"] = "-1"
ENV["JULIA_CONDAPKG_LOG"] = "error"
ENV["JULIA_CONDAPKG_BACKEND"] = "MicroMamba"
using NonLipPTC
