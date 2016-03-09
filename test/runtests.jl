using ShrinkageEstimators
using Base.Test

d = 1
n = 10000
b = sqrt(n) # getfrequencies function uses sqrt(n) bins

arr = rand(d, n)
arrX = arr
arrY = rand(1, 10000)
arrZ = rand(1, 10000)

# Since the estimator is only an estimator, the tests are only approximate.
# In general a tolerance of 0.5 is used, which is relatively large. This
# should pass even if the estimate is slightly off, but fail if something
# is going drastically wrong.

# Check entropy is roughly log_base(b) for uniform distribution with b bins
@test_approx_eq_eps get_shrinkage_entropy(arr) log2(b) 0.5
println("Entropy passed")

# Check change of base works
@test_approx_eq_eps get_shrinkage_entropy(arr, base=e) log(b) 0.5
println("Entropy with change of base passed")

# Check specifying lambda works
@test_approx_eq_eps get_shrinkage_entropy(arr, lambda=0) log2(b) 0.5
println("Entropy with with lambda 0 passed")

# Check uniform width discretization works
@test_approx_eq_eps get_shrinkage_entropy(arr, mode="uniformcount") log2(b) 0.5
println("Entropy with with uniform count passed")

# TODO: Work out how to test "bayesianblocks" mode

# Check joint entropy between multiple identical distributions is roughly
# the same as entropy for one of them
@test_approx_eq_eps get_shrinkage_entropy(arr) get_shrinkage_entropy(arr, arr) 0.5
println("Joint entropy between 2 variables passed")
@test_approx_eq_eps get_shrinkage_entropy(arr) get_shrinkage_entropy(arr, arr, arr) 0.5
println("Joint entropy between 3 variables passed")

# Check I(X; Y) is H(X) + H(Y) - H(X, Y)
a = get_shrinkage_entropy(arrX) + get_shrinkage_entropy(arrY) - get_shrinkage_entropy(arrX, arrY)
@test a == get_shrinkage_mi(arrX, arrY)
println("Mutual information passed")

# Check I(X; Y | Z) is H(X, Z) + H(Y, Z) - H(X, Y, Z) - H(Z)
a = get_shrinkage_entropy(arrX, arrZ) + get_shrinkage_entropy(arrY, arrZ) - get_shrinkage_entropy(arrX, arrY, arrZ) - get_shrinkage_entropy(arrZ)
@test a == get_shrinkage_cmi(arrX, arrY, arrZ)
println("Conditional mutual information passed")
