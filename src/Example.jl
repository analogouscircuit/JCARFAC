module Example

include("JCARFAC.jl")
using .JCARFAC
using Plots; pyplot()
using Formatting

################################################################################
# Plotting Functions 
################################################################################
function plotnap(nap::Array{Float64,2}, fcs::Array{Float64,1}, fs;
                offset=0.05)
    numchannels, numpoints = size(nap)
    ts = (0:numpoints) .* (1/fs)
    ytickvals = 1:4:numchannels
    channels = fcs[ytickvals]
    ytickvals = Float64.(ytickvals) .* offset
    yticklabels = [format("{:.1f}",val) for val ∈ channels]
    p = plot(;
             # ylim=(fcs[1],fcs[end]),
             xlabel="Time (s)",
             ylabel="CF (Hz)", 
             tickfontsize=10,
             guidefontsize=14,
             yticks = (ytickvals, yticklabels),
             legend=false,
             size = (600,400))
    zeroline = zeros(numpoints)
    for c ∈ numchannels:-1:1
        plot!(p, ts, [ (nap[c,:] .+ 0.05*(c-1)) (nap[c,:] .+ 0.05*(c-1))];
              fillrange = [zeroline (nap[c,:] .+ 0.05*(c-1))],
              fillcolor=:white,
              lw=0.4,
              lc=:black)
    end
    p
end


################################################################################
# Main Script
################################################################################

## Generate in input signal -- here a harmonic complex
fs = 44100.0
f0 = 120.0
dur = 0.250
numharmonics = 12
ts = 0:1/fs:dur
sig = zeros(length(ts))
for k ∈ 1:numharmonics
    sig .+= (1.0/k)*sin.(2π*f0*k .* ts)
end
windur = 0.020
numsteps = Int(round(fs*windur))
steps = 0:numsteps-1
sig[1:numsteps] .*= sin.((π/2)*steps ./ numsteps)
sig[end-numsteps+1:end] .*= cos.((π/2)*steps ./ numsteps)


## Generate NAP for entire signal
nap, fcs = calcnap(sig, fs)
p_nap = plotnap(nap, fcs, fs)


## Generate SAI movie
sai, fcs = calcsai(sig, fs, defaultparams(72)...;
                   trigwindowtime = 0.02, 
                   advancestep = 0.001,
                   numtrigwindows=3)
anim = @animate for k ∈ Int(round(sai.numframes/2)):sai.numframes
    heatmap(sai.images[k,:,:];
            xticks=:none,
            yticks=:none,
            legend=false,
            c=:thermal, clims=(0,0.1))
end
gif(anim, "saimovie.gif"; fps=25);




end
