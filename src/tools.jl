# # NOT EASY TO IMPLEMENT
# # because I need the initial value of R, computed from equilibrium.
# function add_virus(X::SIRState, b::Float64, f::Float64; kwargs...)
# 	bvec = b * ones(X.parameters.N)
# 	fvec = f * ones(X.parameters.N)
# 	return add_virus(X, bvec, fvec; kwargs...)
# end

# function add_virus(X::SIRState, b::AbstractVector, f::AbstractVector; I0 = 0.)
# 	# Update parameters
# 	@unpack N, M, α, γ, δ, c, C = X.parameters
# 	new_parameters = SIRParameters(; N+1, α, γ, δ, M, c, copy(C))

# 	# Add new virus in each region
# 	new_regions = SIRRegion[]
# 	for r in X.regions
# 		new_r = add_virus(r, b, f, I0)
# 		push!(new_regions, new_r)
# 	end

# 	#
# 	return SIRState(; regions = new_regions, parameters = new_parameters)
# end

# function add_virus(region::SIRRegion, b, f, I0)
# 	N = size(region)

# 	# new cross-immunity
# 	new_K = zeros(Float64, N+1, N+1)
# 	new_K[1:N, 1:N] .= region.K
# 	new_K[1:N, N+1] .= b
# 	new_K[N+1, 1:N] .= f
# 	new_K[N+1, N+1] = 1

# 	#
# end

"""
	cross_immunity(N::Int, variant_list)

Cross-immunity matrix. Each element in `variant_list` is interpreted as a pair of the form
`(b,f)` giving the cross-immunity of this variant to all others.
The value `variant_list[1]` for the first variant is not used: there are no previous
variants, making them useless.

### Example

```
cross_immunity(2, [(b=0.8, f=0.5)])
cross_immunity(2, [(0., 0.), (b=0.8, f=0.5)]) # first value of variant_list ignored
cross_immunity(3, [(0.6, 0.5), (0.3, 0.2)])
"""
function cross_immunity(N::Int, variant_list)
	K = diagm(ones(N))
	# custom function for indexing into variant_list
	idx(b) = length(variant_list) == N ? b : b-1
	for b in 1:N, a in 1:(b-1)
		# b > a --> a ancient & b recent
		β, ϕ = variant_list[idx(b)]
		K[a,b] = β # immunity from b infection to variant a
		K[b,a] = ϕ # immunity from a infection to variant b
	end
	return K
end

"""
	cross_immunity(N::Int, backward::Real, forward::Real)

Cross-immunity matrix with the same backward and forward effect for all variants.
"""
function cross_immunity(N::Int, backward::Real, forward::Real)
	K = diagm(ones(N))
	for a in 1:N, b in (a+1):N
		K[a,b] = backward # cross-immunity from next virus to previous virus (a < b)
		K[b,a] = forward
	end
	return K
end

"""
	cross_immunity(N)

Diagonal cross-immunity matrix (*i.e.* no interaction between variants).
"""
cross_immunity(N::Int) = diagm(ones(N))

"""
	set_infected(r::SIRRegion, a::Int, val)

Set `r.I[a]` to `val`. Useful to introduce a new virus.
"""
function _set_infected!(r::SIRRegion, a::Int, val)
	N = size(r)
	newI = copy(r.I)
	newI[a] = val
	newS = ones(N) - r.K*newI - r.R
	@info r.I newI
	check_population_size(newS, newI, r.R, r.K)

	r.S .= newS
	r.I .= newI

	return nothing
end
"""
	set_infected(X::SIRState, a::Int, val)

Return a copy of `X` with `I[a]` to `val` in all regions.
"""
function set_infected(X::SIRState, a::Int, val)
	X_copy = deepcopy(X)
	for r in regions(X_copy)
		_set_infected!(r, a, val)
	end
	return X_copy
end

"""
	frequency(X::SIRState, region::Int, virus::Int)

Fraction of infections caused by `virus` in `region`.
"""
function frequency(X::SIRState, region, virus)
	return _frequency(X, region, virus)
end
"""
	frequency(X::SIRState, virus::Int)

Frequency of `virus` accross regions.
"""
function frequency(X::SIRState, virus::Int)
	return mean(frequency(X, i, virus) for i in 1:length(regions(X)))
end

_frequency(X::SIRState, i::Int, a) = X.regions[i].I[a] / sum(X.regions[i].I)
function _frequency(X::SIRState, i::Colon, a)
	return [_frequency(X, ii, a) for ii in 1:length(regions(X))]
end

"""
	frequency(sol::SIRSolution, tvals, i, a)

Frequency of infections by virus `a` in region `i`.
"""
function frequency(sol::SIRSolution, tvals, i, a)
	I = sol[tvals, i, :I, :]
	return map(x -> x[a], I) ./ map(sum, I)
end
"""
	frequency(sol::SIRSolution, tvals, a)

Frequency of infections by virus `a` accross regions.
"""
function frequency(sol::SIRSolution, tvals, a)
	mean(frequency(sol, tvals, i, a) for i in 1:sol.parameters.M)
end
