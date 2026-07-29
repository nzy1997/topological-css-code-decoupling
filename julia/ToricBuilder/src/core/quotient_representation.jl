function _monomial_exponent_vector(monomial)
    exp_vecs = collect(AbstractAlgebra.exponent_vectors(monomial))
    if length(exp_vecs) != 1
        error("Expected a monomial, got $(length(exp_vecs)) terms")
    end
    return exp_vecs[1]
end

function _leading_support_and_monomial(term, F)
    for i in 1:rank(F)
        coeff = term[i]
        if !iszero(coeff)
            return i, coeff
        end
    end
    error("Unable to locate a nonzero leading support in module term")
end

function _degree_bounds_from_leading_monomials(leading_monomials, F)
    bounds = Dict{Int, Vector{Int}}()
    for lm in leading_monomials
        support, monom = _leading_support_and_monomial(lm, F)
        exps = _monomial_exponent_vector(monom)
        if !haskey(bounds, support)
            bounds[support] = collect(exps)
        else
            for i in 1:length(exps)
                bounds[support][i] = max(bounds[support][i], exps[i])
            end
        end
    end
    return bounds
end

function _monomial_from_exps(P, vars, exps)
    m = one(P)
    for (var, exp) in zip(vars, exps)
        if exp != 0
            m *= var^exp
        end
    end
    return m
end

function _monomial_basis_from_lmm(leading_monomials, F, P, gb_lmm)
    vars = gens(P)
    bounds = _degree_bounds_from_leading_monomials(leading_monomials, F)
    for i in 1:rank(F)
        if !haskey(bounds, i)
            error("Missing leading monomial for basis vector e[$i]; quotient is likely infinite")
        end
    end

    basis_vecs = basis(F)
    mb = FreeModElem[]
    for support in 1:rank(F)
        deg_bounds = bounds[support]
        ranges = ntuple(i -> 0:deg_bounds[i], length(deg_bounds))
        for exps in Iterators.product(ranges...)
            monom = _monomial_from_exps(P, vars, exps)
            elem = monom * basis_vecs[support]
            if !iszero(normal_form(elem, gb_lmm))
                push!(mb, elem)
            end
        end
    end
    return mb
end

function _build_monomial_index(mb, F, P)
    index = Dict{Tuple{Int, Tuple}, Int}()
    for (i, elem) in enumerate(mb)
        support, monom = _leading_support_and_monomial(elem, F)
        exps = _monomial_exponent_vector(monom)
        index[(support, Tuple(exps))] = i
    end
    return index
end

function _action_matrix_from_basis(mb, gb, F, P, var)
    dim = length(mb)
    F2 = base_ring(P)
    A = zero_matrix(F2, dim, dim)
    index = _build_monomial_index(mb, F, P)

    for (col, elem) in enumerate(mb)
        prod = var * elem
        reduced = normal_form(prod, gb)
        for support in 1:rank(F)
            poly = reduced[support]
            if iszero(poly)
                continue
            end
            for (coeff, monom) in zip(Oscar.coefficients(poly), Oscar.monomials(poly))
                if iszero(coeff)
                    continue
                end
                exps = _monomial_exponent_vector(monom)
                row = index[(support, Tuple(exps))]
                A[row, col] += coeff
            end
        end
    end
    return A
end

"""
    construct_quotient_representation(mat)

Construct the representation matrices of x and y on the quotient module
defined by the excitation matrix `mat`.
"""
function construct_quotient_representation(mat::MatElem)
    R = base_ring(mat)
    f = Oscar._polyringquo(R)
    module_data = _setup_module(mat, f)
    gb = groebner_basis(module_data.M)

    leading_monomials = [Oscar.leading_monomial(g, ordering=gb.ordering) for g in gb]
    LMM, _ = sub(module_data.F, leading_monomials)
    gb_lmm = groebner_basis(LMM; ordering=gb.ordering)

    mb = _monomial_basis_from_lmm(leading_monomials, module_data.F, module_data.P, gb_lmm)
    vars = gens(module_data.P)

    Ax = _action_matrix_from_basis(mb, gb, module_data.F, module_data.P, vars[1])
    Ay = _action_matrix_from_basis(mb, gb, module_data.F, module_data.P, vars[2])
    return Ax, Ay
end

function _validate_anyon_polynomial_vector(poly_vec::Vector)
    isempty(poly_vec) && throw(ArgumentError("poly_vec must be nonempty"))

    R = parent(poly_vec[1])
    for (i, poly) in enumerate(poly_vec)
        if parent(poly) != R
            throw(ArgumentError("poly_vec[$i] has parent $(parent(poly)); expected $R"))
        end
    end

    if n_variables(R) != 2
        throw(ArgumentError("construct_anyon_translation_representation expects a two-variable ring"))
    end

    if base_ring(R) != Oscar.GF(2)
        throw(ArgumentError("construct_anyon_translation_representation expects coefficient field GF(2)"))
    end

    return R
end

function _polynomial_generator_from_terms(P, vars, coeffs, exp_vecs)
    generator = zero(P)
    offsets = [max(0, -minimum(exp_vec[i] for exp_vec in exp_vecs)) for i in 1:length(vars)]

    for (coeff, exp_vec) in zip(coeffs, exp_vecs)
        monom = one(P)
        for i in 1:length(vars)
            exp = exp_vec[i] + offsets[i]
            if exp != 0
                monom *= vars[i]^exp
            end
        end
        generator += P(coeff) * monom
    end

    return generator
end

function _ordinary_polynomial_generators(poly_vec::Vector)
    R = _validate_anyon_polynomial_vector(poly_vec)
    P, vars = polynomial_ring(base_ring(R), string.(symbols(R)))
    is_laurent_input = R isa AbstractAlgebra.Generic.LaurentMPolyWrapRing

    generators = elem_type(P)[]
    for poly in poly_vec
        coeffs = collect(Oscar.coefficients(poly))
        exp_vecs = collect(AbstractAlgebra.exponent_vectors(poly))
        push!(generators, _polynomial_generator_from_terms(P, vars, coeffs, exp_vecs))
    end

    return P, vars, generators, is_laurent_input
end

function _finite_quotient_degree_bounds(leading_monomials)
    n = length(_monomial_exponent_vector(leading_monomials[1]))
    bounds = fill(-1, n)

    for lm in leading_monomials
        exps = collect(_monomial_exponent_vector(lm))
        nonzero_vars = findall(!iszero, exps)
        if length(nonzero_vars) == 1
            i = only(nonzero_vars)
            bounds[i] = max(bounds[i], exps[i])
        end
    end

    if any(<(0), bounds)
        error("Quotient is likely infinite-dimensional; leading monomials do not bound every variable")
    end

    return bounds
end

function _normal_form(poly, reduction_data)
    return normal_form(poly, reduction_data.ideal; ordering=reduction_data.ordering)
end

function _quotient_ring_monomial_basis(P, vars, gb_vec)
    leading_monomials = [Oscar.leading_monomial(g) for g in gb_vec]
    leading_exps = [collect(_monomial_exponent_vector(lm)) for lm in leading_monomials]
    bounds = _finite_quotient_degree_bounds(leading_monomials)

    basis = elem_type(P)[]
    ranges = ntuple(i -> 0:bounds[i], length(bounds))
    for exps in Iterators.product(ranges...)
        if !any(lm_exps -> all(i -> exps[i] >= lm_exps[i], eachindex(exps)), leading_exps)
            push!(basis, _monomial_from_exps(P, vars, exps))
        end
    end
    return basis
end

function _build_quotient_basis_index(basis)
    index = Dict{Tuple, Int}()
    for (i, monom) in enumerate(basis)
        exps = _monomial_exponent_vector(monom)
        index[Tuple(exps)] = i
    end
    return index
end

function _quotient_ring_action_matrix(basis, reduction_data, P, var)
    dim = length(basis)
    F2 = base_ring(P)
    A = zero_matrix(F2, dim, dim)
    index = _build_quotient_basis_index(basis)

    for (col, monom) in enumerate(basis)
        reduced = _normal_form(var * monom, reduction_data)
        if iszero(reduced)
            continue
        end

        for (coeff, term) in zip(Oscar.coefficients(reduced), Oscar.monomials(reduced))
            if iszero(coeff)
                continue
            end
            exps = _monomial_exponent_vector(term)
            row = index[Tuple(exps)]
            A[row, col] += coeff
        end
    end

    return A
end

"""
    construct_anyon_translation_representation(poly_vec)

Compute the finite quotient-ring basis of `F2[x, y] / ideal(poly_vec)` and the
translation action matrices for multiplication by `x` and `y`.

Inputs may be ordinary two-variable polynomials or Laurent polynomials over
`GF(2)`. Laurent generators are first multiplied by monomial units to clear
negative exponents.
"""
function construct_anyon_translation_representation(poly_vec::Vector)
    P, vars, generators, is_laurent_input = _ordinary_polynomial_generators(poly_vec)
    I = ideal(P, generators)
    if is_laurent_input
        I = saturation(I, ideal(P, [vars[1] * vars[2]]))
    end
    groebner_basis_data = groebner_basis(I)
    gb_vec = collect(groebner_basis_data)
    reduction_data = (; ideal=I, ordering=getfield(groebner_basis_data, :ord))

    basis = _quotient_ring_monomial_basis(P, vars, gb_vec)
    Ax = _quotient_ring_action_matrix(basis, reduction_data, P, vars[1])
    Ay = _quotient_ring_action_matrix(basis, reduction_data, P, vars[2])

    return (
        ;
        groebner_basis=groebner_basis_data,
        basis=basis,
        anyon_count=length(basis),
        Ax=Ax,
        Ay=Ay,
    )
end

"""
    construct_period_table(Ax, Ay)

Compute the period table and period bound from representation matrices Ax and Ay.
"""
function construct_period_table(Ax::MatElem, Ay::MatElem)
    F = base_ring(Ax)
    if base_ring(Ay) != F
        throw(ArgumentError("Ax and Ay must be over the same base ring"))
    end

    dim = nrows(Ax)
    I = identity_matrix(F, dim)

    lx = 1
    while Ax^lx != I
        lx += 1
    end
    ly = 1
    while Ay^ly != I
        ly += 1
    end
    l = lcm(lx, ly)

    period_table = Tuple{Int, Int}[]
    for ix in 0:l
        for iy in 0:l
            if ix == 0 && iy == 0
                continue
            end
            if Ax^ix * Ay^iy == I
                push!(period_table, (ix, iy))
            end
        end
    end
    return period_table, l
end
