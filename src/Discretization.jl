# Functions for performing discretization

using Discretizers

function getfrequencies(values::Array{Float64,2}, mode)
	binids, numberofbins = getbinids(values, mode)
	frequencies = zeros(Int, (1, numberofbins))
	for binid in binids
		frequencies[binid] += 1
	end
	return frequencies
end

function getjointfrequencies(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, mode)
	n = size(valuesX)[2]
	binidsX, numberofbins = getbinids(valuesX, mode)
	binidsY = getbinids(valuesY, mode)[1]
	frequencies = zeros(Int, (numberofbins, numberofbins))
	for i in 1:n
		frequencies[binidsY[i], binidsX[i]] += 1
	end
	return frequencies
end
function getjointfrequencies(valuesX::Array{Float64,2}, valuesY::Array{Float64,2}, valuesZ::Array{Float64,2}, mode)
	n = size(valuesX)[2]
	binidsX, numberofbins = getbinids(valuesX, mode)
	binidsY = getbinids(valuesY, mode)[1]
	binidsZ = getbinids(valuesZ, mode)[1]
	frequencies = zeros(Int, (numberofbins, numberofbins, numberofbins))
	for i in 1:n
		frequencies[binidsZ[i], binidsY[i], binidsX[i]] += 1
	end
	return frequencies
end

function getbinids(values::Array{Float64,2}, mode)
	# This is a separate function because might want to offer different methods in the future
	function getnumberofbins(values)
		# size(values)[2] is n
		return round(Int, sqrt(size(values)[2]))
	end

	numberofbins = getnumberofbins(values)
	numberofdimensions, n = size(values)
	binids = zeros(Int, size(values))
	for i in 1:numberofdimensions
		# If values are all the same, assign them all to bin 1
		min, max = extrema(values)
		if min == max
			binids[i:i, 1:end] += convert(Array{Int}, values) + 1
		elseif mode == "uniformwidth"
			numberofbins = getnumberofbins(values)
			binids[i:i, 1:end] += encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(numberofbins), values)), values)
		elseif mode == "uniformcount"
			numberofbins = getnumberofbins(values)
			binids[i:i, 1:end] += encode(LinearDiscretizer(binedges(DiscretizeUniformCount(numberofbins), reshape(values, length(values)))), values)
		end
	end
	return binids, numberofbins
end
