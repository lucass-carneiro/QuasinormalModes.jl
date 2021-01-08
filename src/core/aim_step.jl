"""
    AIMStep!(p::AIMProblem{N,T}, c::AIMCache{N,U}) where {N <: Unsigned, T <: Number, U <: Any}

Performs a single step of the AIM algorithm:
1. The initial data arrays are not altered.
2. The previous arrays receive the values of the current arrays.
3. The results of the next step computed using the initial and current arrays. Results are stored in the buffer arrays.
4. The current arrays receive the contents of the buffer arrays.

# Input
- `p::AIMProblem`: The problem data to use in the computation.
- `c::AIMCache`: The cache of arrays that corresponds to the problem p.

# Output
    nothing
"""
function AIMStep!(p::AIMProblem{N,T}, c::AIMCache{N,U}) where {N <: Unsigned, T <: Number, U <: Any}

    zeroN = zero(N)
    oneN = one(N)
    twoN = oneN + oneN
    
    # The currently computed coefficients will be the previous iteration coefficients
    # after the iteration is done, so we store the current values to the previous array
    copy!(c.pcda, c.ccda)
    copy!(c.pdda, c.cdda)
    
    # Applay the AIM formula to compute the coeficients writing the results to the
    # buffer arrays and reading data from all the other arrays.
    # The use of a buffer array as temporary storage to the coefficients allow for the
    # computation of the coefficients to be parallelized
    Threads.@threads for i in zeroN:(get_niter(p) - oneN)
        @inbounds sumc = sum(k -> c.icda[k + oneN] * c.ccda[i - k + oneN], zeroN:i)
        @inbounds sumd = sum(k -> c.idda[k + oneN] * c.ccda[i - k + oneN], zeroN:i)

        @inbounds c.bcda[i + oneN] = (i + oneN) * c.ccda[i + twoN] + c.cdda[i + oneN] + sumc
        @inbounds c.bdda[i + oneN] = (i + oneN) * c.cdda[i + twoN] + sumd
    end

    # Copy the data in the buffer arrays to the current arrays
    copy!(c.ccda, c.bcda)
    copy!(c.cdda, c.bdda)

    return nothing
end
