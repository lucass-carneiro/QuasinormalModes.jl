"""
    function AIMStep!(data::QNMData, icda, ccda, pcda, bcda, idda, cdda, pdda, bdda)

Performs a single step of the AIM algorithm:
1. The initial data arrays are not altered.
2. The previous arrays receive the values of the current arrays.
3. The results of the next step computed using the initial and current arrays. Results are stored in the buffer arrays.
4. The current arrays receive the contents of the buffer arrays.

# Input
- `data::QNMData`: The spacetime data to use in the computation.
- `icda`: Initial c data array.
- `ccda`: Current c data array.
- `pcda`: Previous c data array.
- `bcda`: Buffer c data array.
- `idda`: Initial d data array.
- `cdda`: Current d data array.
- `pdda`: Previous d data array.
- `bdda`: Buffer d data array.

# Output
    nothing
"""
function AIMStep!(data::QNMData, icda, ccda, pcda, bcda, idda, cdda, pdda, bdda)
    #= The currently computed coefficients will be the previous iteration coefficients
     = after the iteration is done, so we store the current values to the previous array
    =#
    copy!(pcda, ccda)
    copy!(pdda, cdda)
    
    #= Applay the AIM formula to compute the coeficients writing the results to the
     = buffer arrays and reading data from all the other arrays.
     = The use of a buffer array as temporary storage to the coefficients allow for the
     = computation of the coefficients to be parallelized
    =#
    Threads.@threads for i in 0x00000:(data.nIter - 0x00001)
        sumc = sum(k -> icda[k + 0x00001]*ccda[i - k + 0x00001], 0x00000:i)
        sumd = sum(k -> idda[k + 0x00001]*ccda[i - k + 0x00001], 0x00000:i)

        @inbounds bcda[i + 0x00001] = (i + 0x00001)*ccda[i + 0x00002] + cdda[i + 0x00001] + sumc
        @inbounds bdda[i + 0x00001] = (i + 0x00001)*cdda[i + 0x00002] + sumd
    end

    # Copy the data in the buffer arrays to the current arrays
    copy!(ccda, bcda)
    copy!(cdda, bdda)

    return nothing
end
