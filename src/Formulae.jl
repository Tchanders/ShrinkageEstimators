# Functions that apply the entropy and mutual information formulae

export applyconditionalmutualinformationformula

function applyentropyformula(probabilities, base)
	return -sum([p == 0 ? 0 : p .* log(p) for p in probabilities]) / log(base)
end

function applymutualinformationformula(entropyX, entropyY, jointentropy)
	return entropyX + entropyY - jointentropy
end

function applyconditionalmutualinformationformula(entropyZ, entropyXZ, entropyYZ, entropyXYZ)
	return entropyXZ + entropyYZ - entropyXYZ - entropyZ
end
