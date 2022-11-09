Base.@kwdef struct SIRSolution
	sol
	tspan :: Tuple{Float64, Float64}
	parameters :: SIRParameters
	K :: Vector{Matrix{Float64}}
end

function (sol::SIRSolution)(t::Float64)
	@assert t <= sol.tspan[2] "Time span of solution $(sol.tspan) - got $t"
	return SIRState(sol.sol(t), sol.K, sol.parameters)
end

function Base.getindex(sol::SIRSolution, tvals::AbstractVector, i, g, a)
	idx = SIR_partial_sweep.sir_index(i, g, a, sol.parameters)
	return map(t -> sol.sol(t)[idx], tvals)
end

function simulate(X::SIRState, tspan)
	u0 = vec(X)
	p = (sir = X.parameters, K = [r.K for r in regions(X)])
	return SIRSolution(;
		sol = solve(ODEProblem(SIR!, u0, tspan, p), Tsit5()),
		tspan = tspan,
		parameters = parameters(X),
		K = p.K,
	)
end


function SIR!(du, u, p, t)
	# u is a concatenation of SIR of all regions
	# it is organized as [S1 I1 R1, S2 I2 R2, ..., SM IM RM]
	# - i indexes regions, in 1:p.M
	# - g indexes SIR compartment, in 1:3 or (:S, :I, :R)
	# - a indexes viruses, in 1:N
	# u[sir_index(i, g, a, parameters)] gives the right position in u

	@unpack M, N, α, γ, δ, C = p.sir # K is in p.K
	du .= 0.

	for i in 1:M, a in 1:N
		# S
		idx = sir_index(i, :S, a, N)
		for j in 1:M, b in 1:N
			# region j infecting region i
			# and virus b generating immunity to a
			du[idx] -=
				α * C[i,j] * p.K[i][a,b] * u[sir_index(i,:S,b,N)] * u[sir_index(j,:I,b,N)]
		end
		du[idx] += γ * u[sir_index(i, :R, a, N)]

		# I
		idx = sir_index(i, :I, a, N)
		for j in 1:M
			# region j infecting region i
			du[idx] += α * C[i,j] * u[sir_index(i, :S, a, N)] * u[sir_index(j, :I, a, N)]
		end
		du[idx] -= δ * u[idx]

		# R
		idx = sir_index(i, :R, a, N)
		for b in 1:N
			du[idx] += δ * p.K[i][a,b] * u[sir_index(i, :I, b, N)]
		end
		du[idx] -= γ * u[idx]
	end
	# @assert isapprox(sum(du), 0, atol=1e-15) "$(sum(du))"
	@debug _conserved(du, p), u, p
	return nothing
end

function _conserved(du, parameters)
	return map(1:parameters.sir.M) do i
		map(1:parameters.sir.N) do a
			z = du[sir_index(i, :S, a, parameters.sir)]
			for b in 1:parameters.sir.N
				z += parameters.K[i][a,b] * du[sir_index(i, :I, b, parameters.sir)]
			end
			z += du[sir_index(i, :R, a, parameters.sir)]
		end
	end
end




function _SIR!(du, u, p::SIRParameters, t)
	# u ~ [S1, S2, I1, I2, R1, R2]
	# g is about S/I/R and a is about wt/m
	# --> u[(2*g + a)*M + i]: immune type `i`, virus `a`, compartment `g`
	# the function sir_index(g, a, i, p) takes care of this

	# S
	for a in 1:2, i in 1:p.M
		du[sir_index(:S, a, i, p)] = 0

		# Increase at rate γ*∑_i R_i^a
		du[sir_index(:S, a, i, p)] += p.γ * u[sir_index(:R, a, i, p)]

		# Infections by other regions at rate p.M * α * c_ij S_i^a I_j^a
		for b in 1:2, j in 1:p.M
			du[sir_index(:S, a, i, p)] -=
				p.α * p.M * p.C[i,j] * p.K[i][a,b] * u[sir_index(:S,b,i,p)] * u[sir_index(:I,b,j,p)]
		end
	end

	# I
	for a in 1:2, i in 1:p.M
		du[sir_index(:I, a, i, p)] = 0

		# Exponential decrease at rate δ
		du[sir_index(:I,a,i,p)] = -p.δ * u[sir_index(:I,a,i,p)]


		for j in 1:p.M
			du[sir_index(:I,a,i,p)] +=
				p.α * p.M * p.C[i,j] * u[sir_index(:S,a,i,p)] * u[sir_index(:I,a,j,p)]
		end
	end

	# R
	for a in 1:2, i in 1:p.M
		du[sir_index(:R, a, i, p)] = 0

		# Exponential decrease at rate γ
		du[sir_index(:R,a,i,p)] -= p.γ * u[sir_index(:R,a,i,p)]

		# Increase at rate δ * ∑_b K_i[a,b] I_i^b
		for b in 1:2
			du[sir_index(:R, a, i, p)] +=
				p.δ * p.K[i][a,b] * u[sir_index(:I, b, i, p)]
		end
	end
end
