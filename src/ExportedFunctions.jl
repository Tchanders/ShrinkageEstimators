# Exported functions

export getentropy, getmutualinformation, getconditionalmutualinformation

function getentropy(values::Array{Float64,2}, base=2, lambda=false)
	probabilities = getprobabilitiesshrinkage(getfrequencies(values), lambda)
	return applyentropyformula(probabilities, base)
end
function getentropy(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, base=2, lambda=false)
	jointprobabilities = getprobabilitiesshrinkage(getjointfrequencies(valuesX, valuesY), lambda)
	return applyentropyformula(jointprobabilities, base)
end
function getentropy(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, valuesZ::Array{Float64,2}, base=2, lambda=false)
	jointprobabilities = getprobabilitiesshrinkage(getjointfrequencies(valuesX, valuesY, valuesZ), lambda)
	return applyentropyformula(jointprobabilities, base)
end

function getmutualinformation(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, base=2, lambda=false)
	entropyX = getentropy(valuesX, base, lambda)
	entropyY = getentropy(valuesY, base, lambda)
	entropyXY = getentropy(valuesX, valuesY, base, lambda)
	return applymutualinformationformula(entropyX, entropyY, entropyXY)
end

function getconditionalmutualinformation(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, valuesZ::Array{Float64,2}, base=2, lambda=false)
	entropyZ = getentropy(valuesZ, base, lambda)
	entropyXZ = getentropy(valuesX, valuesZ, base, lambda)
	entropyYZ = getentropy(valuesY, valuesZ, base, lambda)
	entropyXYZ = getentropy(valuesX, valuesY, valuesZ, base, lambda)
	return applyconditionalmutualinformationformula(entropyZ, entropyXZ, entropyYZ, entropyXYZ)
end
