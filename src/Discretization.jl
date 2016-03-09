# Functions for performing discretization

using Discretizers
# using PyCall

# @pyimport astroML.density_estimation as de

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
	binidsX, numberofbinsX = getbinids(valuesX, mode)
	binidsY, numberofbinsY = getbinids(valuesY, mode)
	frequencies = zeros(Int, (numberofbinsY, numberofbinsX))
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

	numberofbins = 0 # So can be accessed outside the loops
	numberofdimensions, n = size(values)
	binids = zeros(Int, size(values))
	for i in 1:numberofdimensions
		# If values are all the same, assign them all to bin 1
		min, max = extrema(values)
		if min == max
			numberofbins += 1
			binids[i:i, 1:end] += convert(Array{Int}, values) + 1
		elseif mode == "uniformwidth"
			numberofbins += getnumberofbins(values)
			binids[i:i, 1:end] += encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(numberofbins), values)), values)
		elseif mode == "uniformwidth2"
			numberofbins += getnumberofbins(values) * 2
			binids[i:i, 1:end] += encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(numberofbins), values)), values)
		elseif mode == "uniformcount"
			numberofbins = getnumberofbins(values)
			try
				binids[i:i, 1:end] += encode(LinearDiscretizer(binedges(DiscretizeUniformCount(numberofbins), reshape(values, length(values)))), values)
			catch
				binids[i:i, 1:end] += encode(LinearDiscretizer(binedges(DiscretizeUniformWidth(numberofbins), values)), values)
			end
		elseif mode == "bayesianblocks"
			# edges = de.bayesian_blocks(reshape(values, length(values)))
			# binsizes = de.histogram(values, edges)[1]
			# binsizes = length(binsizes) == 0 ? [n] : binsizes
			# # The highest index in the current bin
			# highestindex = 0
			# # j is the index of the current bin
			# for j in 1:length(binsizes)
			# 	highestindex += binsizes[j]
			# 	binids[i:i, 1:highestindex] += 1
			# end
			# numberofbins = length(edges) - 1 <= 1 ? 1 : length(edges) - 1
		# TODO: Handle "mode doesn't exist" error
		end
	end
	return binids, numberofbins
end
