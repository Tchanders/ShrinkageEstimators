using ShrinkageEstimators
using Base.Test

# write your own tests here
@test 1 == 1

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
@test_approx_eq_eps getentropy(arr) log2(b) 0.5
println(getentropy(arr) - log2(b))
println("Entropy passed")

# Check change of base works
@test_approx_eq_eps getentropy(arr, e) log(b) 0.5
println("Entropy with change of base passed")

# Check specifying lambda works
@test_approx_eq_eps getentropy(arr, 2, 0) log2(b) 0.5
println("Entropy with with lambda 0 passed")

# Check joint entropy between multiple identical distributions is roughly
# the same as entropy for one of them
@test_approx_eq_eps getentropy(arr) getentropy(arr, arr) 0.5
println("Joint entropy between 2 variables passed")
@test_approx_eq_eps getentropy(arr) getentropy(arr, arr, arr) 0.5
println("Joint entropy between 3 variables passed")

# Check I(X; Y) is H(X) + H(Y) - H(X, Y)
a = getentropy(arrX) + getentropy(arrY) - getentropy(arrX, arrY)
@test a == getmutualinformation(arrX, arrY)
println("Mutual information passed")

# Check I(X; Y | Z) is H(X, Z) + H(Y, Z) - H(X, Y, Z) - H(Z)
a = getentropy(arrX, arrZ) + getentropy(arrY, arrZ) - getentropy(arrX, arrY, arrZ) - getentropy(arrZ)
@test a == getconditionalmutualinformation(arrX, arrY, arrZ)
println("Conditional mutual information passed")
