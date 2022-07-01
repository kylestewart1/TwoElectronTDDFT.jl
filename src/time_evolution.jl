using LinearAlgebra

mutable struct TimeEvolution{P <: AbstractPropagator}
    ψ::Ket
    H::Hamiltonian
    U::P
    Δt::Float64
    n_time_steps::Int
    iter::Int
    t₀::Float64
    t::Float64
    function TimeEvolution{P}(H::Hamiltonian, U::P, Δt::Float64, t₀::Float64, n_time_steps::Int)
        N = size(H)[0]
        Ψ = Ket{Vector{Float64}}(zeros(N))
        t=t₀
        iter = 0
        new(Ψ, H, U, Δt, n_time_steps, iter, t₀, t)
    end
end


function update!(TE::TimeEvolution, H_new::Hamiltonian)
    TE.U(TE.ψ)
    update!(TE.U, H_new)
    TE.H = H_new
    TE.iter += 1
    TE.t += TE.Δt
end


function ham_1D(grid_size, Δx, V)
    off_diag = Diagonal
end

t₀ = 0
Δt = 0.01
n_time_steps = 1000
H = Hamiltonian{Float64}()
U = CrankNicholson(Δt, H)