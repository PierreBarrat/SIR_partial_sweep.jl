"""
	set_infected(r::SIRRegion, a::Int, val)

Set `r.I[a]` to `val`. Useful to introduce a new virus.
"""
function set_infected!(r::SIRRegion, a::Int, val)
	N = size(r)
	newI = copy(r.I)
	newI[a] = val
	newS = ones(N) - r.K*newI - r.R

	check_population_size(newS, newI, r.R, r.K)

	r.S .= newS
	r.I .= newI

	return nothing
end
"""
	set_infected!(X::SIRState, a::Int, val)

Set `I[a]` to `val` in all regions.
"""
function set_infected!(X::SIRState, a::Int, val)
	for r in regions(X)
		set_infected!(r, a, val)
	end
	return nothing
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
	frequency(SIRSolution, i, a)

Frequency of infections by virus `a` in region `i`.
"""
function frequency(sol::SIRSolution, tvals, i, a)
	I = sol[tvals, i, :I, :]
	return map(x -> x[a], I) ./ map(sum, I)
end
function frequency(sol::SIRSolution, tvals, a)
	mean(frequency(sol, tvals, i, a) for i in 1:sol.parameters.M)
end
