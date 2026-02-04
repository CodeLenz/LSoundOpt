using Plots
using LaTeXStrings
using Plots.Measures # Required for specifying margins (mm, cm, px)

# 1. Setup for Publication-Quality Plots
default(
    fontfamily="Computer Modern", 
    framestyle=:box,              
    grid=:alpha,                 
    gridalpha=0.3,
    linewidth=2.0,               
    tickfontsize=11,
    guidefontsize=12,
    legendfontsize=10,
    size=(600, 400),              
    dpi=300,
    # --- FIX FOR CLIPPING ---
    margin=5mm # Adds a 5mm buffer around the entire figure
)

# 2. Extract Data
freqs = x2[1]
spl_initial = x2[2]
spl_opt = x2[3]

# 3. Create the Plot
p = plot()

# Highlight the optimization range [700, 900]
vspan!([700, 900], 
    color=:gray, 
    alpha=0.2, 
    label="Target Band [700-900 Hz]"
)

# Plot the curves
plot!(freqs, spl_initial, 
    label="Initial Design", 
    xlabel="Frequency [Hz]",       
    ylabel="SPL [dB]", 
    color=:dodgerblue,
    legend=:topright,
    xlims=(minimum(freqs), maximum(freqs))
)

plot!(freqs, spl_opt, 
    label="Optimized Topology", 
    color=:firebrick
)

# 4. Save
savefig("response_comparison_highlighted.pdf") 
display(p)