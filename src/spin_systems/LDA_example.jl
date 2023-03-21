include("LDA.jl")

n_qubits = 3

J_perp_list = fill(-1.0, n_qubits-1)
J_perp = nearest_neighbor(J_perp_list)

uniform_Z = LinRange(-1.0,1.0,100)
scaling = [-10.0, -1.0, 0.0, 1.0, 10.0, 100.0]
E = energy_functional(n_qubits, J_perp, scaling, uniform_Z)
kinetic = E[3]

E_hc = []
for i in eachindex(E)
    push!(E_hc, E[i] .- kinetic)
end


using Plots

energy_func_plt = plot(uniform_Z, E[1], label=scaling[1])
for i in 2:(size(E)[1])
    plot!(uniform_Z, E[i], label=scaling[i])
end
xlabel!("σᶻ")
ylabel!("Energy")
title!("Exact energy")
display(energy_func_plt) 

hc_func_plt = plot(uniform_Z, E_hc[1], label=scaling[1])
for i in 2:(size(E_hc)[1])
    plot!(uniform_Z, E_hc[i], label=scaling[i])
end
xlabel!("σᶻ")
ylabel!("E_hc")
title!("Heisenberg plus correlation energy")
display(hc_func_plt) 

"""
using CurveFit


quadratic = []
for i in eachindex(E)
    push!(quadratic, poly_fit(uniform_Z, E[i], 2))
end

E_approx = []
for i in eachindex(quadratic)
    coeff = quadratic[i]
    E_quad = coeff[1] .+ coeff[2] *uniform_Z .+ coeff[3] * (uniform_Z).^2
    push!(E_approx,E_quad)
end

coeff_plt = plot(scaling, getindex.(quadratic,1), label="Constant")
plot!(scaling, getindex.(quadratic,2), label="Linear")
plot!(scaling, getindex.(quadratic,3), label="Quadratic")

"""

E_LDA = []
for i in eachindex(scaling)
    J_par_list = scaling[i] * J_perp_list
    push!(E_LDA, LDA(uniform_Z, J_perp_list, J_par_list))
end


diff = E_LDA - E
for i in eachindex(diff)
    diff[i] = abs.(diff[i])
end

LDA_plt = plot(uniform_Z, E_LDA[1], label=scaling[1])
for i in 2:length(scaling)
    plot!(uniform_Z, E_LDA[i], label=scaling[i])
end
ylabel!("Energy")
xlabel!("σᶻ")
title!("LDA total energy")

diff_plt = plot(uniform_Z, diff[1], label=scaling[1])
for i in 2:length(scaling)
    plot!(uniform_Z, diff[i], label=scaling[i])
end
ylabel!("Energy error")
xlabel!("σᶻ")
title!("Error")

plt = plot(energy_func_plt, LDA_plt, diff_plt, layout=(3,1))

savefig(plt, "LDA_energy.png")