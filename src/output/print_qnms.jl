"""
    function printQNMs(qnms::Array{Complex{BigFloat}, 1};
                       cutoff::BigFloat = BigFloat(1.0e-5),
                       instab::Bool = false
                       )

Prints the quasinormal modes previously computed using computeQNMs to stdout.

# Input
- `qnms::Array{Complex{BigFloat}, 1}`: Array of QNMs.
- `cutoff::BigFloat`: All QNMs whose real parts are smaller than this value are not printed.
- `instab::Bool`: If this flag is set to true, the output file will also contain instabilities (QNMs with positive imaginary part).

# Output
    nothing
"""
function printQNMs(qnms::Array{Complex{BigFloat}, 1};
                   cutoff::BigFloat = BigFloat(1.0e-5),
                   instab::Bool = false
                   )
    println("# Re(omega)    Im(omega)")

    for qnm in qnms
        if real(qnm) > cutoff && ( instab ? true : imag(qnm) < big"0.0" )
        println(real(qnm), "    ", imag(qnm))
        end
    end

    return nothing
end
