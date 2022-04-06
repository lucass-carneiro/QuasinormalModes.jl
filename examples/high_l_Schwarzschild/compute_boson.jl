# --------------------------------------------------------------------
# Compute Schwarzschild QNMs perturbed by spin 0,1 or 2 field and save
# --------------------------------------------------------------------

function compute_boson(iter, s)
    x0 = Complex(big"0.43", big"0.0")
    M = Complex(big"1.0", big"0.0")

    p_num = Schwarzschild_boson(iter, x0, M, 0x00000, convert(UInt32, s))
    c_num = AIMCache(p_num)

    file = open("qnm_" * Dates.format(now(), "yyyy-mm-dd HH:MM:SS") * "_s_$(s).dat", "w")
    println(file, "# 1:Iterations 2:l 3:n 4:real(guess) 5:imag(guess) 6:re(omega) 7:im(omega)")

    if s == 0
       refArray = spin_0
    elseif s == 1
       refArray = spin_1
    elseif s == 2
       refArray = spin_2
    else
        error("Bosonic spin must be either 0, 1 or 2")
    end

    println("Starting computation")

    for reference in refArray
        l = UInt32(reference[1])
        n = reference[2]

        println("Computing l = ", l, " n = ", n)

        guess = Complex(BigFloat(reference[3]), BigFloat(reference[4]))

        print(file, iter, "    ", l, "    ", n, "    ", reference[3], "    ", reference[4],"    ")

        p_num = Schwarzschild_boson(iter, x0, M, l, convert(UInt32, s))
        qnm = computeEigenvalues(Serial(), p_num, c_num, guess, nls_xtol = big"1.0e-20", nls_ftol = big"1.0e-20", nls_iterations = 10000000)

        if(qnm.x_converged || qnm.f_converged)
            print(file, qnm.zero[1], "    ", qnm.zero[2], "\n")
        else
            print(file, "Did not converge    Did not converge\n")
        end

        flush(file)
        flush(stdout)
    end

    close(file)
end
