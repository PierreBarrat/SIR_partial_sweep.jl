function equilibrium_1(X::SIRState)
	@assert X.parameters.M == 1 "Only for one region"
	region = _equilibrium(X.regions[1], X.parameters)
	return SIRState(; regions = [region], parameters = deepcopy(X.parameters))
end

function _equilibrium(r::SIRRegion, p::SIRParameters)
	S, I, R = _equilibrium(r.K, p)
	return SIRRegion(; S, I, R, K=copy(r.K))
end
function _equilibrium(K::Matrix, p::SIRParameters)
	@unpack N, α, γ, δ = p

	S = δ/α * ones(N)
	I = try
		γ * (1-δ/α) / (δ + γ) * (K \ ones(N))
	catch err
		@error "Singular `K`: rank = $(rank(K)) < $N. Cannot compute equilibrium." K
		error(err)
	end
	R = 1 .- S - K*I

	return S, I, R
end

"""
!! THIS IS FALSE FOR S and R (but should be good for I)
"""
function equilibrium_M(X::SIRState)
	@unpack M, N, α, γ, δ, C = X.parameters

	eq_regions = []
	for (i, region) in enumerate(X.regions)
		Z = sum(C[i,:])
		K_av = sum(C[i,j] * r.K for (j,r) in enumerate(X.regions)) / Z

		S = δ/α * ones(N)
		I = try
			γ * (1-δ/α) / (δ + γ) * (K_av \ ones(N))
		catch err
			@error "Singular `K_av`: rank = $(rank(K_av)) < $N. Cannot compute equilibrium." K_av
			error(err)
		end
		R = 1 .- S - region.K*I

		push!(eq_regions, SIRRegion(; S, I, R, K=copy(region.K)))
	end

	return SIRState(; regions = eq_regions, parameters = deepcopy(X.parameters))
end

"""
	equilibrium(X::SIRState)

Return an `SIRState` object giving the equilibrium reached by the system if initiated
from `X`.
"""
function equilibrium(X::SIRState)
	return if length(regions(X)) == 1
		equilibrium_1(X)
	else
		equilibrium_M(X)
	end
end
