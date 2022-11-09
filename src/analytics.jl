function equilibrium_1(X::SIRState)
	@assert parameters(X).M == 1 "Only for one region"

	@unpack N, α, γ, δ = parameters(X)
	K = regions(X)[1].K

	S = δ/α * ones(N)
	I = γ * (1-δ/α) / (δ + γ) * (K \ ones(N))
	R = 1 .- S - K*I
	region = SIRRegion(; S, I, R, K)

	return SIRState(;regions = [region], parameters = deepcopy(parameters(X)))
end

function equilibrium_M(X::SIRState)
	return nothing
end

function equilibrium(X::SIRState)
	return if length(regions(X)) == 1
		equilibrium_1(X)
	else
		equilibrium_M(X)
	end
end
