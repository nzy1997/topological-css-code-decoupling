# Color Code Example
# ====================
# This example demonstrates how to use ToricBuilder to construct and analyze
# a color code using toric building blocks.

using ToricBuilder
using Oscar

function main(; show_progress::Bool=false)
    # Step 1: Setup the Laurent polynomial ring over GF(2)
    # ----------------------------------------------------
    # We work over GF(2) (binary field) with variables x and y representing
    # translations in two spatial directions
    Rxy, (x, y) = laurent_polynomial_ring(GF(2), ["x", "y"])

    # Step 2: Define the polynomial vector for the color code
    # --------------------------------------------------------
    # These polynomials define the stabilizer structure of the quantum code
    # Each polynomial represents a constraint on the qubits
    poly_vector = [1 + x + x*y, 1 + y + x*y]

    # Step 3: Find the lattice parameter L
    # -------------------------------------
    # L determines the size of the finite lattice (torus) on which the code lives
    # find_L finds the minimal L such that the ideal relations are satisfied
    l = find_L([x, y], ideal(Rxy, poly_vector))

    # Step 4: Generate the L-table (optional, for analysis)
    # ------------------------------------------------------
    # This creates a table showing how the ideal relations behave on the L×L torus
    L_table = make_L_table([x, y], ideal(Rxy, poly_vector), l)

    # Step 5: Create the excitation matrix
    # -------------------------------------
    # The excitation matrix encodes the structure of the code's check operators
    A = excitation_matrix(poly_vector)

    # Step 6: Apply coarse graining
    # ------------------------------
    # Coarse graining reduces the system to an l×l periodic lattice
    # This step converts the infinite lattice to a finite torus
    A_cg = coarse_grain(A, [x, y], [l, l])

    # Step 7: Convert to toric form with the high-level workflow
    # -----------------------------------------------------------
    result = build_toric_form(poly_vector; show_progress=show_progress, compute_inverse=true)

    # Step 8: Verify the result
    # --------------------------
    check_result(result, result.input_matrix; require_inverse=true)

    size(A_cg) == (2 * l^2, 4 * l^2) || error("unexpected coarse_grain output shape")
    result.l == l || error("build_toric_form returned an unexpected period bound")

    return (; poly_vector, l, L_table, A, A_cg, result)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
