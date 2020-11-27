"""
    function filterQNMs!(qnms::Array{Complex{BigFloat}, 1},
                         cutoff::BigFloat = BigFloat(1.0e-5),
                         instab::Bool = false
                         )

Filter the array of quasinormal modes previously computed using computeQNMs **in-place**.

# Input
- `qnms::Array{Complex{BigFloat}, 1}`: Array of QNMs to be filtered.
- `cutoff::BigFloat`: All QNMs whose real parts are smaller than this value are excluded.
- `instab::Bool`: If this flag is set to true, the output file will also contain instabilities (QNMs with positive imaginary part).

# Output
    nothing
"""
function filterQNMs!(qnms::Array{Complex{BigFloat}, 1},
                     cutoff::BigFloat = BigFloat(1.0e-5),
                     instab::Bool = false
                     )
    
    function qnm_filter(x::Complex{BigFloat})
        if real(x) > cutoff && ( instab ? true : imag(x) < big"0.0" )
            return true
        else
            return false
        end
    end
    
    filter!(qnm_filter, qnms)
    
    return nothing
end
