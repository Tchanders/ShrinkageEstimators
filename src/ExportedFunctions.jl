# Exported functions

# TODO:
# Offer alternative target
# Look into standardising the documentation
# Offer alternative discretisation methods

export get_shrinkage_entropy, get_shrinkage_mi, get_shrinkage_cmi


"""
Calculates a James-Stein shrinkage estimate for the entropy of a set
of observed values, via the following steps:
(1) convert the observed values to frequencies of discrete bins
(2) convert the frequencies to probablities using shrinkage estimation
(3) run the probabilities through the entropy formula.
Can calculate joint entropy for up to 3 sets of observed values.

Parameters:

valuesX - dxn Array{Float64,2} - The observed values, where d is
the number of dimensions and n is the number of values.

[valuesY] - dxn Array{Float64,2} - The observed values, where d is
the number of dimensions and n is the number of values.

[valuesZ] - dxn Array{Float64,2} - The observed values, where d is
the number of dimensions and n is the number of values.

[base=2] - Int - The base of the logarithm (this determines the units).

[lambda=false] - Int - The shrinkage intensity, between 0 and 1
inclusive. If omitted, the optimal lambda it will be calculated.
"""
function get_shrinkage_entropy(valuesX::Array{Float64,2}, base=2, lambda=false)
	probabilities = getprobabilitiesshrinkage(getfrequencies(valuesX), lambda)
	return applyentropyformula(probabilities, base)
end
function get_shrinkage_entropy(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, base=2, lambda=false)
	jointprobabilities = getprobabilitiesshrinkage(getjointfrequencies(valuesX, valuesY), lambda)
	return applyentropyformula(jointprobabilities, base)
end
function get_shrinkage_entropy(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, valuesZ::Array{Float64,2}, base=2, lambda=false)
	jointprobabilities = getprobabilitiesshrinkage(getjointfrequencies(valuesX, valuesY, valuesZ), lambda)
	return applyentropyformula(jointprobabilities, base)
end

"""
Calculates a James-Stein shrinkage estimate for the mutual information
between two sets of observed values.

Parameters:

valuesX - dxn Array{Float64,2} - A set of observed values, where d is
the number of dimensions and n is the number of values.

valuesY - dxn Array{Float64,2} - Another set of observed values, where
d is the number of dimensions and n is the number of values.

[base=2] - Int - The base of the logarithm (this determines the units).

[lambda=false] - Int - The shrinkage intensity, between 0 and 1
inclusive. If omitted, the optimal lambda it will be calculated.
"""
function get_shrinkage_mi(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, base=2, lambda=false)
	entropyX = get_shrinkage_entropy(valuesX, base, lambda)
	entropyY = get_shrinkage_entropy(valuesY, base, lambda)
	entropyXY = get_shrinkage_entropy(valuesX, valuesY, base, lambda)
	return applymutualinformationformula(entropyX, entropyY, entropyXY)
end

"""
Calculates a James-Stein shrinkage estimate for the conditional mutual
information between two sets of observed values, given a third set of
observed values.

Parameters:

valuesX - dxn Array{Float64,2} - A set of observed values, where d is
the number of dimensions and n is the number of values.

valuesY - dxn Array{Float64,2} - Another set of observed values, where
d is the number of dimensions and n is the number of values.

valuesZ - dxn Array{Float64,2} - The set of observed values, that will
be conditioned on, where d is the number of dimensions and n is the
number of values.

[base=2] - Int - The base of the logarithm (this determines the units).

[lambda=false] - Int - The shrinkage intensity, between 0 and 1
inclusive. If omitted, the optimal lambda it will be calculated.
"""
function get_shrinkage_cmi(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, valuesZ::Array{Float64,2}, base=2, lambda=false)
	entropyZ = get_shrinkage_entropy(valuesZ, base, lambda)
	entropyXZ = get_shrinkage_entropy(valuesX, valuesZ, base, lambda)
	entropyYZ = get_shrinkage_entropy(valuesY, valuesZ, base, lambda)
	entropyXYZ = get_shrinkage_entropy(valuesX, valuesY, valuesZ, base, lambda)
	return applyconditionalmutualinformationformula(entropyZ, entropyXZ, entropyYZ, entropyXYZ)
end
