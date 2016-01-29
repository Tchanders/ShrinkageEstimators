# Functions for performing discretization

function getfrequencies(values::Array{Float64,2})
	binids, numberofbins = getbinids(values)
	frequencies = zeros(Int, (1, numberofbins))
	for binid in binids
		frequencies[binid] += 1
	end
	return frequencies
end

function getjointfrequencies(valuesX::Array{Float64,2}, valuesY::Array{Float64,2})
	n = size(valuesX)[2]
	binidsX, numberofbins = getbinids(valuesX)
	binidsY = getbinids(valuesY)[1]
	frequencies = zeros(Int, (numberofbins, numberofbins))
	for i in 1:n
		frequencies[binidsY[i], binidsX[i]] += 1
	end
	return frequencies
end
function getjointfrequencies(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, valuesZ::Array{Float64,2})
	n = size(valuesX)[2]
	binidsX, numberofbins = getbinids(valuesX)
	binidsY = getbinids(valuesY)[1]
	binidsZ = getbinids(valuesZ)[1]
	frequencies = zeros(Int, (numberofbins, numberofbins, numberofbins))
	for i in 1:n
		frequencies[binidsZ[i], binidsY[i], binidsX[i]] += 1
	end
	return frequencies
end

function getbinids(values::Array{Float64,2})
	# This is a separate function because might want to offer different methods in the future
	function getnumberofbins(values)
		# size(values)[2] is n
		return round(Int, sqrt(size(values)[2]))
	end

	# This is a separate function because might want to offer different methods in the future
	function getvalueadjustmentparameters(values, numberofbins)
		min, max = extrema(values)
		binwidth = (max - min) / (numberofbins - 1)
		return min, binwidth
	end

	function adjustvalues(values, numberofbins)
		min, binwidth = getvalueadjustmentparameters(values, numberofbins)
		# If all the values were 0, they don't need adjusting
		# For now, test this by testing whether binwidth is 0
		if binwidth == 0
			return convert(Array{Int}, values) + 1
		end
		return floor(Int, (values - min) * (1 / binwidth)) + 1 # This should be ceil?
	end

	numberofbins = getnumberofbins(values)
	numberofdimensions, n = size(values)
	binids = zeros(Int, size(values))
	for i in 1:numberofdimensions
		binids[i:i, 1:end] += adjustvalues(values, numberofbins)
	end
	return binids, numberofbins
end
