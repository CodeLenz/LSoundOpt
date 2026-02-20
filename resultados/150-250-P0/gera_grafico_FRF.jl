using Plots
using LaTeXStrings
using Plots.Measures 

# espaçamento das frquências 
esp = 50

# Faixa de frequência 
banda = [150,250]

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
    margin=5mm 
)

# 2. Extract Data
freqs = x2[1]
spl_initial = x2[2]
spl_opt = x2[3]



# 3. Create the Plot
p = plot()

# Highlight the optimization range 
vspan!(banda, 
    color=:gray, 
    alpha=0.2, 
    label="Target Band [$(banda[1])-$(banda[2]) Hz]"
)

# Plot the curves
plot!(freqs, spl_initial, 
    label="Initial Design", 
    xlabel="Frequency [Hz]",       
    ylabel="SPL [dB]", 
    color=:dodgerblue,
    legend=:topright,
    xlims=(minimum(freqs), maximum(freqs)),
    xticks=minimum(freqs):esp:maximum(freqs)
)

plot!(freqs, spl_opt, 
    label="Optimized Topology", 
    color=:firebrick
)

# 4. Save
savefig("FRF_$(banda[1])_$(banda[2]).pdf") 
display(p)