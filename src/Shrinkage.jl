# Functions for performing shrinkage

function getprobabilitiesshrinkage(frequencies::Array{Int}, lambda::Void)
	function getlambda(n::Int, normalizedvalues, target)
		# Unbiased estimator of variance of u
		varu = normalizedvalues .* (1 - normalizedvalues) / (n - 1)
		msp = sum((normalizedvalues - target).^2) # misspecification ???

		# Estimate shrinkage intensity
		lambda = msp == 0 ? 1 : sum(varu) / msp
		
		# Make lambda be between 0 and 1 inclusive
		return lambda > 1 ? 1 : (lambda < 0 ? 0 : lambda)
	end

	target = 1 / length(frequencies) # Target is uniform distribution
	n = sum(frequencies)
	normalizedvalues = frequencies / n
	lambda = n == 1 || n == 0 ? 1 : getlambda(n, normalizedvalues, frequencies)
	
	return lambda * target + (1 - lambda) * normalizedvalues
end
function getprobabilitiesshrinkage(frequencies::Array{Int}, lambda::Real)
	target = 1 / length(frequencies) # Target is uniform distribution
	n = sum(frequencies)
	normalizedcounts = frequencies / n

	return lambda * target + (1 - lambda) * normalizedcounts
end
