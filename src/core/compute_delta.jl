"""
    computeDelta(data::QuadFreqData)
    
Compute and return the AIM "quantization condition".

# Input
- `data::QuadFreqData`: A quadratic frequency space-time data structure.

# Output
An object of type `Polynomial{Complex{BigFloat}}` whose roots are the quasinormal modes.
"""
function computeDelta(data::QuadFreqData)
    # Pre-allocate the coeficient arrays. We need 4 arrays for each AIM iteration equation, thus 4 arrays of the same size.
    icda = zeros(Polynomial{Complex{BigFloat}}, data.nIter + 0x00001) # hold the initial c data, i.e., c^i_0
    ccda = zeros(Polynomial{Complex{BigFloat}}, data.nIter + 0x00001) # hold the coefficients for the current aim step, c^i_n
    pcda = zeros(Polynomial{Complex{BigFloat}}, data.nIter + 0x00001) # hold the coefficients for the previous aim step, c^i_{n-1}
    bcda = zeros(Polynomial{Complex{BigFloat}}, data.nIter + 0x00001) # The work buffer used to actually compute the c coefficients in parallel

    idda = zeros(Polynomial{Complex{BigFloat}}, data.nIter + 0x00001) # hold the initial d data, i.e., c^i_0
    cdda = zeros(Polynomial{Complex{BigFloat}}, data.nIter + 0x00001) # hold the coefficients for the current aim step, d^i_n
    pdda = zeros(Polynomial{Complex{BigFloat}}, data.nIter + 0x00001) # hold the coefficients for the previous aim step, d^i_{n-1}
    bdda = zeros(Polynomial{Complex{BigFloat}}, data.nIter + 0x00001) # The work buffer used to actually compute the d coefficients in parallel


    # The initial data arrays are initialized using the Fourier coefficient formula f^(k)(x0)/k!
    for i in 0x00000:data.nIter
        icda[i + 0x00001] = createPoly(data, i, 0x00001)/factorial(big(i))
        idda[i + 0x00001] = createPoly(data, i, 0x00002)/factorial(big(i))
    end

    # Initialize the current data arryas with the initial data for the first step
    copy!(ccda, icda)
    copy!(cdda, idda)

    # Perform the aim steps
    for i in 0x00001:data.nIter
        AIMStep!(data, icda, ccda, pcda, bcda, idda, cdda, pdda, bdda)
    end

    # Compute and return the AIM "quantization condition"
    return cdda[0x00001]*pcda[0x00001] - pdda[0x00001]*ccda[0x00001]
end

"""
    computeDelta(data::GenFreqData, ω::Complex{BigFloat})

Compute and return the AIM "quantization condition".

# Input
- `data::GenFreqData`: A generic frequency space-time data structure.
- `ω::Complex{BigFloat}`: The value of `ω` at which the quantization condition is to be evaluated.

# Output
An object of type `Complex{BigFloat}` that is zero if `ω` is a quasinormal mode.
"""
function computeDelta(data::GenFreqData, ω::Complex{BigFloat})
    # Pre-allocate the coeficient arrays. We need 4 arrays for each AIM iteration equation, thus 4 arrays of the same size.
    ccda = zeros(Complex{BigFloat}, data.nIter + 0x00001) # hold the coefficients for the current aim step, c^i_n
    pcda = zeros(Complex{BigFloat}, data.nIter + 0x00001) # hold the coefficients for the previous aim step, c^i_{n-1}
    bcda = zeros(Complex{BigFloat}, data.nIter + 0x00001) # The work buffer used to actually compute the c coefficients in parallel
    
    cdda = zeros(Complex{BigFloat}, data.nIter + 0x00001) # hold the coefficients for the current aim step, d^i_n
    pdda = zeros(Complex{BigFloat}, data.nIter + 0x00001) # hold the coefficients for the previous aim step, d^i_{n-1}
    bdda = zeros(Complex{BigFloat}, data.nIter + 0x00001) # The work buffer used to actually compute the d coefficients in parallel
    
    # The initial data array elements are the fourier expansion coefficients of λ0 and S0 in x around x0  of order nIter given ω
    t = Complex{BigFloat}(BigFloat(data.x0), BigFloat(0.0)) + Taylor1(Complex{BigFloat}, Int64(data.nIter))
    icda::Array{Complex{BigFloat},1} = data(0x00001,t,ω)
    idda::Array{Complex{BigFloat},1} = data(0x00002,t,ω)

    # Initialize the current data arryas with the initial data for the first step
    copy!(ccda, icda)
    copy!(cdda, idda)

    # Perform the aim steps
    for i in 0x00001:data.nIter
        AIMStep!(data, icda, ccda, pcda, bcda, idda, cdda, pdda, bdda)
    end

    # Compute and return the AIM "quantization condition"
    return cdda[0x00001]*pcda[0x00001] - pdda[0x00001]*ccda[0x00001]
end
