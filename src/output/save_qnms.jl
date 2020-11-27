"""
    saveQNMs(qnms::Array{Complex{BigFloat}, 1};
             filename::String = "qnms.dat",
             cutoff::BigFloat = BigFloat(1.0e-5),
             instab::Bool = false
             )

Save the quasinormal modes previously computed using computeQNMs.

# Input
- `qnms::Array{Complex{BigFloat}, 1}`: Array of QNMs.
- `filename::String`: The name of the file with the computed QNMs.
- `cutoff::BigFloat`: All QNMs whose real parts are smaller than this value are not saved in the output file.
- `instab::Bool`: If this flag is set to true, the output file will also contain instabilities (QNMs with positive imaginary part).

# Output
    nothing
"""
function saveQNMs(qnms::Array{Complex{BigFloat}, 1};
                  filename::String = "qnms.dat",
                  cutoff::BigFloat = BigFloat(1.0e-5),
                  instab::Bool = false
                  )
    
    file = open(filename, "w")

    println(file, "# AIM Iterations = ", data.nIter)
    println(file, "# Re(omega)    Im(omega)")

    for qnm in qnms
        if real(qnm) > cutoff && ( instab ? true : imag(qnm) < big"0.0" )
        println(file, real(qnm), "    ", imag(qnm))
        end
    end

    close(file)

    return nothing
end
