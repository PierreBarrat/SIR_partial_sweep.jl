################## SIRParameters ##################

"""
	struct SIRParameters

Parameters of the SIR.

## Arguments
- `N`: number of virus variants.
- `α`, `γ`, `δ`: "base" parameters of the SIR (simplified equations below):
```math
dS/dt = -αSI + γR

dI/dt = αSI - δI
```

- `M`: number of geographic regons
- `c`: strength of connection between region, in `[0,1]`
- `C`: connectivity matrix set from `c`: `c` off-diagonal and `1 - M*c` on diagonal
"""
@with_kw struct SIRParameters
	# General
	N :: Int = 2 # Number of viruses

	# SIR
	α :: Float64 = 10
	γ :: Float64 = .01
	δ :: Float64 = 1

	# Geographic structure
	M :: Int = 1 # Number of geographic regions
	c :: Float64 = .1 # Strength of connection between region, in `[0,1]`
	C :: Matrix{Float64} = c * ones(M, M) .+ diagm(ones(M)*(1 - M*c))
end


################## SIRRegion ##################

Base.@kwdef struct SIRRegion
	S :: Vector{Float64}
	I :: Vector{Float64}
	R :: Vector{Float64}
	K :: Matrix{Float64}

	function SIRRegion(S, I, R, K)
		@assert length(S) == length(I) == length(R) == size(K,1) == size(K,2) "Size mismatch"
		@assert check_population_size(S, I, R, K) "Size of host pop. is not one"
		return new(S, I, R, K)
	end
end

"""
	SIRRegion(
		N::Int, I::AbstractVector, R::AbstractVector = zeros(N);
		b=0., f=0., K = cross_immunity(N, b, f)
	)
	SIRRegion(N::Int; I0 = 0., R0 = 0., b = 0., f = 0., K = cross_immunity(N, b, f))

`N`: number of viruses.
"""
function SIRRegion(
	N::Int, I::AbstractVector, R::AbstractVector = zeros(N);
	b=0., f=0., K = cross_immunity(N, b, f)
)
	S = ones(N) - K*I - R

	for (a, (s, i, r)) in enumerate(zip(S, I, R))
		c = 0 <= s <= 1 && 0 <= i <= 1 && 0 <= r <= 1
		@assert c "Problem with pop. size - virus $a has (S, I, R) = $((s, i, r))"
	end

	return SIRRegion(; S, I, R, K)
end
function SIRRegion(N::Int; I0 = 0., R0 = 0., b = 0., f = 0., K = cross_immunity(N, b, f))
	return SIRRegion(N, I0 * ones(N), R0 * ones(N); b, f, K)
end


Base.size(X::SIRRegion) = length(X.S)
Base.vec(X::SIRRegion) = vcat(X.S, X.I, X.R)

function check_population_size(region::SIRRegion)
	mapreduce(x -> isapprox(x, 1, rtol=1e-10), *, population(region), init=true)
end
function check_population_size(S, I, R, K)
	mapreduce(*, 1:length(S), init=true) do a
		isapprox(population(S, I, R, K, a), 1, rtol=1e-10)
	end
end
"""
	population(region::SIRRegion, a::Int)

Size of host population when divided in compartment regarding
"""
population(S, I, R, K, a) = S[a] + sum(K[a,b] * I[b] for b in 1:length(S)) + R[a]
population(r::SIRRegion, a) = population(r.S, r.I, r.R, r.K, a)
population(r::SIRRegion) = map(a -> population(r, a), 1:size(r))


################## SIRState ##################

"""
	SIRState(; regions, parameters)
"""
Base.@kwdef mutable struct SIRState
	parameters :: SIRParameters
	regions :: Vector{SIRRegion}

	function SIRState(parameters, regions)
		@assert !isempty(regions) "Must have at least one region"
		@assert size(regions[1]) == parameters.N
		@assert length(regions) == parameters.M
		@assert mapreduce(*, regions, init=true) do r
			size(r) == size(regions[1])
		end "Regions do not match in size."

		return new(parameters, regions)
	end
end

"""
	SIRState(parameters::SIRParameters; I0=0., R0=0., b=0., f=0.)
"""
function SIRState(parameters::SIRParameters; I0=0., R0=0., b=0., f=0.)
	regions = [SIRRegion(parameters.N; I0, R0, b, f) for m in 1:parameters.M]
	return SIRState(parameters, regions)
end
"""
	SIRState(dat::Vector{Float64}, Ks, parameters::SIRParameters)

Create `SIRState` from a vector and an array of cross-immunity matrices `Ks`.
"""
function SIRState(dat::Vector{Float64}, Ks, parameters::SIRParameters)
	@assert length(dat) == parameters.M * parameters.N * 3
	regions = SIRRegion[]
	for i in 1:parameters.M
		r = SIRRegion(;
			S = dat[sir_index(i, :S, :, parameters)],
			I = dat[sir_index(i, :I, :, parameters)],
			R = dat[sir_index(i, :R, :, parameters)],
			K = Ks[i],
		)
		push!(regions, r)
	end
	return SIRState(; parameters, regions)
end

function Base.getindex(X::SIRState, i::Int, g, a)
	g = compartment_to_symb(g)
	return getfield(X.regions[i], g)[a]
end

"""
	vec(X::SIRState)

Transform `X::SIRState` into a vector for use in differential equation solvers.
"""
Base.vec(X::SIRState) = vcat(map(vec, X.regions)...)

regions(X::SIRState) = X.regions
parameters(X::SIRState) = X.parameters



################## Utils ##################

function symb_to_virus(a::Symbol)
	return if a == :wt
		1
	elseif a == :m
		2
	else
		error("Unknown symbol $a")
	end
end
symb_to_virus(a::Number) = a

function symb_to_compartment(g::Symbol)
	if g == :S
		return 1
	elseif g == :I
		return 2
	elseif g == :R
		return 3
	else
		error()
	end
end
symb_to_compartment(g::Number) = g
compartment_to_symb(g::Symbol) = g
function compartment_to_symb(g::Number)
	return if g == 1
		:S
	elseif g == 2
		:I
	elseif g == 3
		:R
	else
		error()
	end
end


"""
	sir_index(i, g, a, N::Int)

Index of strain `a` in SIR class `g` in region `i`.

- `a` is a number in `[1, N]` (or `:wt`/`:m` if `N==2`)
- `g` can be any of `(:S, :I, :R)`
- `i` is a number in `[1, M]`: it's the region

## Drawing
```
[
	// i --> Region i
	[S_1, S_2, ..., S_N] // g == 1
	[I_1, I_2, ..., I_N] // g == 2
	[R_1, R_2, ..., R_N] // g == 3

	etc...
]
```

- Linear size of each region: `3*N`
- Total size of array `3*N*M` if `M` is the number of regions
"""
function sir_index(i, g, a, N::Int)
	# @assert 1 <= i <= M "Error accessing region $i in size $M SIR model"
	# @assert 1 <= a <= N "Error accessing region $i in size $M SIR model"

	a = symb_to_virus.(a)
	g = symb_to_compartment.(g)

	region_offset = 3*N * (i .- 1)
	compartment_offset = (g .- 1) * N

	# @info "region_offset: $(region_offset), compartment_offset = $(compartment_offset)"

	return region_offset .+ compartment_offset .+ (a)
end
sir_index(i, g, ::Colon, N::Int) = sir_index(i, g, 1:N, N)

# function sir_index(i, g, ::Colon, N::Int)
# 	a = symb_to_virus(a)
# 	g = symb_to_compartment(g)

# 	region_offset = 3*N * (i .- 1)
# 	compartment_offset = (g .- 1) * N

# 	return region_offset .+ compartment_offset .+ (1:N)
# end

sir_index(i, g, a, p::SIRParameters) = sir_index(i, g, a, p.N)
sir_index(i, g, a, s::SIRState) = sir_index(i, g, a, s.parameters)
sir_index(::Colon, g, a, p::SIRParameters) = sir_index(1:p.M, g, a, p)

