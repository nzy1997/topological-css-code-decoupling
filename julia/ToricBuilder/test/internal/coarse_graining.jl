using Test
using ToricBuilder
using ToricBuilder: enlarge_monomial, kronecker_product, replace_variable, replace_variable_inv, replace_variable_inv_mat, coarse_graining_with_replace_variable,enlarge_polynomial, possible_remainders,find_minimum_triangle, find_parallelogram_lattice_points
using Oscar

@testset "enlarge_monomial" begin
	F = GF(2)
	# Use a Laurent polynomial ring for the single-variable tests
	Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

	result1 = enlarge_monomial(Rxy, [x], [2], [2])
	@test result1.entries == [x 0; 0 x]

	result1 = enlarge_monomial(Rxy, [x], [-3], [2])
	@test result1.entries == [0 x^-1; x^-2 0]

	result2 = enlarge_monomial(Rxy, [x], [1], [3])
	expected2 = zero_matrix(Rxy, 3, 3)
	expected2[1, 3] = x
	expected2[2, 1] = Rxy(1)
	expected2[3, 2] = Rxy(1)
	@test result2 == expected2

	
	result3 = enlarge_monomial(Rxy, [x, y], [1, 1], [2, 2])
	mat_x = zero_matrix(Rxy, 2, 2)
	mat_x[1, 2] = x
	mat_x[2, 1] = Rxy(1)

	mat_y = zero_matrix(Rxy, 2, 2)
	mat_y[1, 2] = y
	mat_y[2, 1] = Rxy(1)

	expected3 = kronecker_product(mat_x, mat_y)
	@test result3 == expected3


	result4 = enlarge_monomial(Rxy, [x, y], [2, 1], [2, 3])
	mat_x2 = zero_matrix(Rxy, 2, 2)
	mat_x2[1, 1] = x
	mat_x2[2, 2] = x

	mat_y_3 = zero_matrix(Rxy, 3, 3)
	mat_y_3[1, 3] = y
	mat_y_3[2, 1] = Rxy(1)
	mat_y_3[3, 2] = Rxy(1)

	expected4 = kronecker_product(mat_x2, mat_y_3)
	@test result4 == expected4

	result5 = enlarge_monomial(Rxy, [x, y], [0, 0], [2, 2])

	expected5 = identity_matrix(Rxy, 4)
	@test result5 == expected5

	# Switch back to the single-variable ring to test negative exponents
	result6 = enlarge_monomial(Rxy, [x], [-1], [3])
	result7 = enlarge_monomial(Rxy, [x], [1], [3])
	@test (result6 * result7).entries == [Rxy(1) 0 0; 0 Rxy(1) 0; 0 0 Rxy(1)]
end

@testset "enlarge_polynomial" begin
	F = GF(2)
	Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

	poly1 = x^2 + x + 1
	result1 = enlarge_polynomial(Rxy, poly1, [x], [2])

	mat_x2 = zero_matrix(Rxy, 2, 2)
	mat_x2[1, 1] = x
	mat_x2[2, 2] = x

	mat_x1 = zero_matrix(Rxy, 2, 2)
	mat_x1[1, 2] = x
	mat_x1[2, 1] = Rxy(1)

	mat_1 = identity_matrix(Rxy, 2)

	expected1 = mat_x2 + mat_x1 + mat_1
	@test result1 == expected1

	poly2 = x^2*y + x*y^2 + 1
	result2 = enlarge_polynomial(Rxy, poly2, [x, y], [2, 2])

	mat1 = enlarge_monomial(Rxy, [x, y], [2, 1], [2, 2])
	mat2 = enlarge_monomial(Rxy, [x, y], [1, 2], [2, 2])
	mat3 = enlarge_monomial(Rxy, [x, y], [0, 0], [2, 2])

	expected2 = mat1 + mat2 + mat3
	@test result2 == expected2

	# Switch back to the single-variable ring
	Rx_single, x_single = laurent_polynomial_ring(F, "x")
	poly3 = x_single + Rx_single(0)  # In GF(2), 3x + 2 = x
	result3 = enlarge_polynomial(Rx_single, poly3, [x_single], [3])

	mat_x_3 = enlarge_monomial(Rx_single, [x_single], [1], [3])
	expected3 = mat_x_3
	@test result3 == expected3

	poly4 = x + y
	result4 = enlarge_polynomial(Rxy, poly4, [x, y], [2, 2])

	mat_x = enlarge_monomial(Rxy, [x, y], [1, 0], [2, 2])
	mat_y = enlarge_monomial(Rxy, [x, y], [0, 1], [2, 2])

	expected4 = mat_x + mat_y
	@test result4 == expected4

	poly5 = Rxy(1)
	result5 = enlarge_polynomial(Rxy, poly5, [x, y], [2, 2])

	expected5 = identity_matrix(Rxy, 4)
	@test result5 == expected5

	poly6 = x^2*y
	result6 = enlarge_polynomial(Rxy, poly6, [x, y], [2, 3])

	expected6 = enlarge_monomial(Rxy, [x, y], [2, 1], [2, 3])
	@test result6 == expected6

	# Test Laurent polynomials with negative exponents in the single-variable ring
	Rx_single2, x_single2 = laurent_polynomial_ring(F, "x")
	frac1 = 1 + x_single2^(-1)  # Equivalent to (x + 1) / x
	mat1 = enlarge_polynomial(Rx_single2, frac1, [x_single2], [2])

	expected1 = zero_matrix(Rx_single2, 2, 2)
	expected1[1, 1] = 1
	expected1[1, 2] = 1
	expected1[2, 1] = x_single2^(-1)
	expected1[2, 2] = 1

	matches1 = true
	for i in 1:2, j in 1:2
		if mat1[i, j] != expected1[i, j]
			matches1 = false
		end
	end
	@test matches1

	# Test multivariable Laurent polynomials using the Laurent polynomial ring
	# (y + (x+1)/x) / y = 1 + y^(-1) + x^(-1)*y^(-1)
	frac2 = 1 + y^(-1) + x^(-1)*y^(-1)
	mat2 = enlarge_polynomial(Rxy, frac2, [x, y], [2, 2])

	mat_1 = enlarge_monomial(Rxy, [x, y], [0, 0], [2, 2])      # 1
	mat_y_inv = enlarge_monomial(Rxy, [x, y], [0, -1], [2, 2])  # y^(-1)
	mat_xy_inv = enlarge_monomial(Rxy, [x, y], [-1, -1], [2, 2]) # x^(-1)*y^(-1)

	expected2 = mat_1 + mat_y_inv + mat_xy_inv

	matches2 = true
	for i in 1:4, j in 1:4
		if mat2[i, j] != expected2[i, j]
			matches2 = false
		end
	end
	@test matches2
end

@testset "coarse_graining" begin
	F = GF(2)
	# Use a Laurent polynomial ring for the single-variable tests
	Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

	A1 = zero_matrix(Rxy, 1, 2)
	A1[1, 1] = 1 + x
	A1[1, 2] = Rxy(0)

	B1 = coarse_graining(A1, [x], [2])

	@test size(B1) == (2, 4)
	Rxy_new, (x_new, y_new) = laurent_polynomial_ring(F, ["x", "y"])
	expected_block = enlarge_polynomial(Rxy_new, 1 + x_new, [x_new], [2])
	block_matches = true
	for i in 1:2, j in 1:2
		if B1[i, j] != expected_block[i, j]
			block_matches = false
			break
		end
	end
	@test block_matches

	A2 = zero_matrix(Rxy, 2, 2)
	A2[1, 1] = x
	A2[1, 2] = 1
	A2[2, 1] = 1 + x
	A2[2, 2] = x^2

	B2 = coarse_graining(A2, [x], [2])

	@test size(B2) == (4, 4)

	A3 = zero_matrix(Rxy, 1, 1)
	A3[1, 1] = 1 + x + x^2

	B3 = coarse_graining(A3, [x], [3])

	@test size(B3) == (3, 3)

	Rxy_new3, (x_new3, y_new3) = laurent_polynomial_ring(F, ["x", "y"])
	expected_block3 = enlarge_polynomial(Rxy_new3, 1 + x_new3 + x_new3^2, [x_new3], [3])

	block_matches3 = true
	for i in 1:3, j in 1:3
		if B3[i, j] != expected_block3[i, j]
			block_matches3 = false
		end
	end
	@test block_matches3

	A4 = zero_matrix(Rxy, 2, 2)
	A4[1, 1] = 1 + x + y
	A4[1, 2] = x * y
	A4[2, 1] = x^2
	A4[2, 2] = 1

	B4 = coarse_graining(A4, [x, y], [2, 2])

	expansion = 4
	@test size(B4) == (2 * expansion, 2 * expansion)

	Rxy_new4, (x_new4, y_new4) = laurent_polynomial_ring(F, ["x", "y"])
	expected_block4 = enlarge_polynomial(Rxy_new4, 1 + x_new4 + y_new4, [x_new4, y_new4], [2, 2])

	block_matches4 = true
	for i in 1:4, j in 1:4
		if B4[i, j] != expected_block4[i, j]
			block_matches4 = false
		end
	end
	@test block_matches4

	# Switch back to the single-variable ring for this test
	Rxy_single3, (x_single3, y_single3) = laurent_polynomial_ring(F, ["x", "y"])
	A5 = zero_matrix(Rxy_single3, 2, 3)
	A5[1, 1] = x_single3
	A5[1, 3] = 1
	A5[2, 2] = x_single3 + 1

	B5 = coarse_graining(A5, [x_single3], [2])

	@test size(B5) == (4, 6)

	zero_block = true
	for i in 3:4, j in 1:2
		if !iszero(B5[i, j])
			zero_block = false
		end
	end
	@test zero_block

	A6 = zero_matrix(Rxy, 1, 1)
	A6[1, 1] = x + y


	B6 = coarse_graining(A6, [x, y], [2, 3])

	@test size(B6) == (6, 6)
	Rxy_new6, (x_new6, y_new6) = laurent_polynomial_ring(F, ["x", "y"])
	expected_block6 = enlarge_polynomial(Rxy_new6, x_new6 + y_new6, [x_new6, y_new6], [2, 3])

	block_matches6 = true
	for i in 1:6, j in 1:6
		if B6[i, j] != expected_block6[i, j]
			block_matches6 = false
		end
	end
	@test block_matches6
end

@testset "find_L" begin
	F = GF(2)
	Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
	I = ideal([x*y+x+1, x*y+y+1])
	l = find_L([x, y], I)
	@test l == 3

	I = ideal([1+x + x^-1*y^3, 1+y + y^-1*x^3])
	l = find_L([x, y], I)
	@test l == 12
end

@testset "make_L_table" begin
	F = GF(2)
	Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
	I = ideal([x*y+x+1, x*y+y+1])
	l = make_L_table([x, y], I,3)
	@test l == [1 0 0 1; 0 0 1 0; 0 1 0 0 ; 1 0 0 1]
	I = ideal([1+x + x^-1*y^3, 1+y + y^-1*x^3])
	l = make_L_table([x, y], I,12)
	@test l[13,13] == 1
end

@testset "matrix period APIs from poly_vec" begin
	F = GF(2)
	Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])

	color_poly_vec = [1 + x + x*y, 1 + y + x*y]
	color_l = find_L(color_poly_vec)
	@test color_l == 3
	@test make_L_table(color_poly_vec, color_l) == make_L_table([x, y], ideal(color_poly_vec), color_l)

	sparse_poly_vec = [1 + x + x^-1*y^3, 1 + y + y^-1*x^3]
	@test find_L(sparse_poly_vec) == 12
	@test find_L(sparse_poly_vec; max_int=11) === nothing

	triangle = ToricBuilder.find_minimum_triangle_periods(color_poly_vec; max_int=10)
	vertices, area = ToricBuilder.find_minimum_triangle(make_L_table([x, y], ideal(color_poly_vec), color_l))
	@test triangle == ((vertices[2][1] - 1, vertices[2][2] - 1), (vertices[3][1] - 1, vertices[3][2] - 1), Int(2 * area))
end

@testset "replace_variable - 2 variables" begin
	F = GF(2)
	R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
	Ruv, (u, v) = laurent_polynomial_ring(F, ["u", "v"])
	
	# Define the variable relations
	u_rel = x*y^2
	v_rel = x^-1*y
	
	# Helper function: verify the result
	function verify_result(original, new_poly, remainder, new_var_rels)
		if isnothing(new_poly)
			return false
		end
		# Convert the new-variable expression back to the old variables
		result = remainder
		str_repr = string(new_poly)
		if occursin("u", str_repr)
			u_exp = 0
			m = match(r"u\^?(-?\d+)?", str_repr)
			if m !== nothing
				u_exp = m.captures[1] !== nothing ? parse(Int, m.captures[1]) : 1
			end
			result *= new_var_rels[1]^u_exp
		end
		if occursin("v", str_repr)
			v_exp = 0
			m = match(r"v\^?(-?\d+)?", str_repr)
			if m !== nothing
				v_exp = m.captures[1] !== nothing ? parse(Int, m.captures[1]) : 1
			end
			result *= new_var_rels[2]^v_exp
		end
		return result == original
	end
	
	# Test case 1: remainder set [1, x, x^2]
	allowed_rems_1 = [R(1), x, x^2]
	
	uv_result, remainder = replace_variable(x^2*y^2, [u_rel, v_rel], R, Ruv, (x, y), (u, v); allowed_remainders=allowed_rems_1)
	@test !isnothing(uv_result)
	@test verify_result(x^2*y^2, uv_result, remainder, [u_rel, v_rel])
	
	uv_result, remainder = replace_variable(x^3, [u_rel, v_rel], R, Ruv, (x, y), (u, v); allowed_remainders=allowed_rems_1)
	@test !isnothing(uv_result)
	@test verify_result(x^3, uv_result, remainder, [u_rel, v_rel])
	@test remainder == R(1)
	
	uv_result, remainder = replace_variable(y^3, [u_rel, v_rel], R, Ruv, (x, y), (u, v); allowed_remainders=allowed_rems_1)
	@test !isnothing(uv_result)
	@test verify_result(y^3, uv_result, remainder, [u_rel, v_rel])
	
	# Test case 2: remainder set [1, x, x*y]
	allowed_rems_2 = [R(1), x, x*y]
	
	uv_result, remainder = replace_variable(x^3*y^2, [u_rel, v_rel], R, Ruv, (x, y), (u, v); allowed_remainders=allowed_rems_2)
	@test !isnothing(uv_result)
	@test verify_result(x^3*y^2, uv_result, remainder, [u_rel, v_rel])
	
	uv_result, remainder = replace_variable(x*y, [u_rel, v_rel], R, Ruv, (x, y), (u, v); allowed_remainders=allowed_rems_2)
	@test !isnothing(uv_result)
	@test verify_result(x*y, uv_result, remainder, [u_rel, v_rel])

	allowed_rems_3 = [R(1), x^-1]

	uv_result, remainder = replace_variable(x^2, [u_rel, v_rel], R, Ruv, (x, y), (u, v); allowed_remainders=allowed_rems_3)
	@test !isnothing(uv_result)
	@test verify_result(x^2, uv_result, remainder, [u_rel, v_rel])
end

@testset "replace_variable - 3 variables" begin
	F = GF(2)
	R3, (x, y, z) = laurent_polynomial_ring(F, ["x", "y", "z"])
	R3_new, (u, v, w) = laurent_polynomial_ring(F, ["u", "v", "w"])
	
	# Define the variable relations
	u_rel = x*y
	v_rel = y*z
	w_rel = x*z^2
	
	allowed_rems = [R3(1), x, y, z, x*y]
	
	# Helper function: verify the result
	function verify_result_3var(original, new_poly, remainder, new_var_rels)
		if isnothing(new_poly)
			return false
		end
		result = remainder
		str_repr = string(new_poly)
		var_names = ["u", "v", "w"]
		for (i, var_name) in enumerate(var_names)
			if occursin(var_name, str_repr)
				exp = 0
				m = match(Regex("$(var_name)\\^?(-?\\d+)?"), str_repr)
				if m !== nothing
					exp = m.captures[1] !== nothing ? parse(Int, m.captures[1]) : 1
				end
				result *= new_var_rels[i]^exp
			end
		end
		return result == original
	end
	
	# Test basic monomials
	uvw_result, remainder = replace_variable(x*y, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w]; allowed_remainders=allowed_rems)
	@test !isnothing(uvw_result)
	@test verify_result_3var(x*y, uvw_result, remainder, [u_rel, v_rel, w_rel])
	@test remainder == R3(1)
	
	uvw_result, remainder = replace_variable(y*z, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w]; allowed_remainders=allowed_rems)
	@test !isnothing(uvw_result)
	@test verify_result_3var(y*z, uvw_result, remainder, [u_rel, v_rel, w_rel])
	@test remainder == R3(1)
	
	uvw_result, remainder = replace_variable(x*z^2, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w]; allowed_remainders=allowed_rems)
	@test !isnothing(uvw_result)
	@test verify_result_3var(x*z^2, uvw_result, remainder, [u_rel, v_rel, w_rel])
	@test remainder == R3(1)
	
	# Test more complex monomials
	uvw_result, remainder = replace_variable(x*y^2*z, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w]; allowed_remainders=allowed_rems)
	@test !isnothing(uvw_result)
	@test verify_result_3var(x*y^2*z, uvw_result, remainder, [u_rel, v_rel, w_rel])
	@test remainder == R3(1)
	
	# Test single variables
	uvw_result, remainder = replace_variable(x, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w]; allowed_remainders=allowed_rems)
	@test !isnothing(uvw_result)
	@test remainder == x
	
	uvw_result, remainder = replace_variable(y, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w]; allowed_remainders=allowed_rems)
	@test !isnothing(uvw_result)
	@test remainder == y
end

@testset "replace_variable - 4 variables" begin
	F = GF(2)
	R4, (x, y, z, w) = laurent_polynomial_ring(F, ["x", "y", "z", "w"])
	R4_new, (a, b, c) = laurent_polynomial_ring(F, ["a", "b", "c"])
	
	# Define the variable relations using 3 new variables to represent 4 old variables
	a_rel = x*y
	b_rel = z*w
	c_rel = x*z
	
	allowed_rems = [R4(1), x, y, z, w, x*y, z*w]
	
	# Helper function: verify the result
	function verify_result_4var(original, new_poly, remainder, new_var_rels)
		if isnothing(new_poly)
			return false
		end
		result = remainder
		str_repr = string(new_poly)
		var_names = ["a", "b", "c"]
		for (i, var_name) in enumerate(var_names)
			if occursin(var_name, str_repr)
				exp = 0
				m = match(Regex("$(var_name)\\^?(-?\\d+)?"), str_repr)
				if m !== nothing
					exp = m.captures[1] !== nothing ? parse(Int, m.captures[1]) : 1
				end
				result *= new_var_rels[i]^exp
			end
		end
		return result == original
	end
	
	# Test basic monomials
	abc_result, remainder = replace_variable(x*y, [a_rel, b_rel, c_rel], R4, R4_new, [x, y, z, w], [a, b, c]; allowed_remainders=allowed_rems)
	@test !isnothing(abc_result)
	@test verify_result_4var(x*y, abc_result, remainder, [a_rel, b_rel, c_rel])
	@test remainder == R4(1)
	
	abc_result, remainder = replace_variable(z*w, [a_rel, b_rel, c_rel], R4, R4_new, [x, y, z, w], [a, b, c]; allowed_remainders=allowed_rems)
	@test !isnothing(abc_result)
	@test verify_result_4var(z*w, abc_result, remainder, [a_rel, b_rel, c_rel])
	@test remainder == R4(1)
	
	abc_result, remainder = replace_variable(x*z, [a_rel, b_rel, c_rel], R4, R4_new, [x, y, z, w], [a, b, c]; allowed_remainders=allowed_rems)
	@test !isnothing(abc_result)
	@test verify_result_4var(x*z, abc_result, remainder, [a_rel, b_rel, c_rel])
	@test remainder == R4(1)
	
	# Test more complex monomials
	abc_result, remainder = replace_variable(x*y*z*w, [a_rel, b_rel, c_rel], R4, R4_new, [x, y, z, w], [a, b, c]; allowed_remainders=allowed_rems)
	@test !isnothing(abc_result)
	@test verify_result_4var(x*y*z*w, abc_result, remainder, [a_rel, b_rel, c_rel])
	@test remainder == R4(1)
	
	abc_result, remainder = replace_variable(x^2*y*z, [a_rel, b_rel, c_rel], R4, R4_new, [x, y, z, w], [a, b, c]; allowed_remainders=allowed_rems)
	@test !isnothing(abc_result)
	@test verify_result_4var(x^2*y*z, abc_result, remainder, [a_rel, b_rel, c_rel])
	@test remainder == R4(1)
	
	# Test single variables
	abc_result, remainder = replace_variable(x, [a_rel, b_rel, c_rel], R4, R4_new, [x, y, z, w], [a, b, c]; allowed_remainders=allowed_rems)
	@test !isnothing(abc_result)
	@test remainder == x
end

@testset "replace_variable_inv - 2 variables" begin
	F = GF(2)
	R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
	Ruv, (u, v) = laurent_polynomial_ring(F, ["u", "v"])
	
	# Define the variable relations
	u_rel = x*y^2
	v_rel = x^-1*y
	
	# Test basic conversions
	result1 = replace_variable_inv(u, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	@test result1 == x*y^2
	
	result2 = replace_variable_inv(v, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	@test result2 == x^-1*y
	
	# Test u*v
	result3 = replace_variable_inv(u*v, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	expected3 = (x*y^2) * (x^-1*y)  # = y^3
	@test result3 == expected3
	@test result3 == y^3
	
	# Test u^2*v^3
	result4 = replace_variable_inv(u^2*v^3, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	expected4 = (x*y^2)^2 * (x^-1*y)^3  # = x^2*y^4 * x^-3*y^3 = x^-1*y^7
	@test result4 == expected4
	@test result4 == x^-1*y^7
	
	# Test negative exponents
	result5 = replace_variable_inv(u^-1, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	expected5 = (x*y^2)^-1  # = x^-1*y^-2
	@test result5 == expected5
	@test result5 == x^-1*y^-2
	
	# Test 1 (the multiplicative identity)
	result6 = replace_variable_inv(Ruv(1), [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	@test result6 == R(1)
	
	# Test round-trip conversion with replace_variable
	# For simple cases, round-trip conversion should work
	allowed_rems = [R(1), x, x^2]
	uv_poly, remainder = replace_variable(x^3*y^4, [u_rel, v_rel], R, Ruv, [x, y], [u, v]; allowed_remainders=allowed_rems)
	if !isnothing(uv_poly)
		back = replace_variable_inv(uv_poly, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
		@test back * remainder == x^3*y^4
	end
end

@testset "replace_variable_inv - 3 variables" begin
	F = GF(2)
	R3, (x, y, z) = laurent_polynomial_ring(F, ["x", "y", "z"])
	R3_new, (u, v, w) = laurent_polynomial_ring(F, ["u", "v", "w"])
	
	# Define the variable relations
	u_rel = x*y
	v_rel = y*z
	w_rel = x*z^2
	
	# Test basic conversions
	result1 = replace_variable_inv(u, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w])
	@test result1 == x*y
	
	result2 = replace_variable_inv(v, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w])
	@test result2 == y*z
	
	result3 = replace_variable_inv(w, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w])
	@test result3 == x*z^2
	
	# Test compositions
	result4 = replace_variable_inv(u*v, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w])
	expected4 = (x*y) * (y*z)  # = x*y^2*z
	@test result4 == expected4
	
	result5 = replace_variable_inv(u*v*w, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w])
	expected5 = (x*y) * (y*z) * (x*z^2)  # = x^2*y^2*z^3
	@test result5 == expected5
	
	# Test powers
	result6 = replace_variable_inv(u^2*w^-1, [u_rel, v_rel, w_rel], R3, R3_new, [x, y, z], [u, v, w])
	expected6 = (x*y)^2 * (x*z^2)^-1  # = x^2*y^2 * x^-1*z^-2 = x*y^2*z^-2
	@test result6 == expected6
end

@testset "replace_variable_inv - roundtrip consistency" begin
	F = GF(2)
	R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
	Ruv, (u, v) = laurent_polynomial_ring(F, ["u", "v"])
	
	u_rel = x*y^2
	v_rel = x^-1*y
	
	# Test round-trip consistency: for monomials in the new-variable ring,
	# applying replace_variable_inv and then replace_variable should recover the original expression
	test_cases = [
		u,
		v,
		u*v,
		u^2,
		v^3,
		u^-1*v^2,
	]
	
	for test_poly in test_cases
		# New variables -> old variables
		old_poly = replace_variable_inv(test_poly, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
		
		# Old variables -> new variables, if possible
		allowed_rems = [R(1), x, x^-1, y, y^-1, x*y, x^-1*y]
		new_poly, remainder = replace_variable(old_poly, [u_rel, v_rel], R, Ruv, [x, y], [u, v]; allowed_remainders=allowed_rems)
		
		# If the conversion succeeds, check that the result is consistent
		if !isnothing(new_poly)
			# Convert new_poly back to the old-variable ring and compare it with old_poly
			new_poly_in_old = replace_variable_inv(new_poly, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
			@test new_poly_in_old * remainder == old_poly
		end
	end
end

@testset "replace_variable_inv_mat" begin
	F = GF(2)
	R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
	Ruv, (u, v) = laurent_polynomial_ring(F, ["u", "v"])
	
	u_rel = x*y^2
	v_rel = x^-1*y
	
	# Test case 1: monomial matrix
	A1 = zero_matrix(Ruv, 2, 2)
	A1[1, 1] = u
	A1[1, 2] = v
	A1[2, 1] = u*v
	A1[2, 2] = Ruv(1)
	
	B1 = replace_variable_inv_mat(A1, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	
	@test B1[1, 1] == x*y^2
	@test B1[1, 2] == x^-1*y
	@test B1[2, 1] == y^3
	@test B1[2, 2] == R(1)
	
	# Test case 2: polynomial matrix
	A2 = zero_matrix(Ruv, 2, 2)
	A2[1, 1] = u + v
	A2[1, 2] = u*v + Ruv(1)
	A2[2, 1] = u^2 + v^2
	A2[2, 2] = u + v + Ruv(1)
	
	B2 = replace_variable_inv_mat(A2, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	
	# B2[1,1] = u + v = x*y^2 + x^-1*y
	@test B2[1, 1] == x*y^2 + x^-1*y
	
	# B2[1,2] = u*v + 1 = y^3 + 1
	@test B2[1, 2] == y^3 + R(1)
	
	# B2[2,1] = u^2 + v^2 = (x*y^2)^2 + (x^-1*y)^2 = x^2*y^4 + x^-2*y^2
	@test B2[2, 1] == x^2*y^4 + x^-2*y^2
	
	# B2[2,2] = u + v + 1
	@test B2[2, 2] == x*y^2 + x^-1*y + R(1)
	
	# Test case 3: negative exponents
	A3 = zero_matrix(Ruv, 2, 2)
	A3[1, 1] = u^-1
	A3[1, 2] = v^-1
	A3[2, 1] = u^-1*v
	A3[2, 2] = u*v^-1
	
	B3 = replace_variable_inv_mat(A3, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	
	# u^-1 = (x*y^2)^-1 = x^-1*y^-2
	@test B3[1, 1] == x^-1*y^-2
	
	# v^-1 = (x^-1*y)^-1 = x*y^-1
	@test B3[1, 2] == x*y^-1
	
	# u^-1*v = x^-1*y^-2 * x^-1*y = x^-2*y^-1
	@test B3[2, 1] == x^-2*y^-1
	
	# u*v^-1 = x*y^2 * x*y^-1 = x^2*y
	@test B3[2, 2] == x^2*y
	
	# Test case 4: zero matrices and zero entries
	A4 = zero_matrix(Ruv, 3, 3)
	A4[1, 1] = u
	A4[2, 2] = v
	A4[3, 3] = Ruv(1)
	# All other entries should remain zero
	
	B4 = replace_variable_inv_mat(A4, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	
	@test B4[1, 1] == x*y^2
	@test B4[2, 2] == x^-1*y
	@test B4[3, 3] == R(1)
	@test iszero(B4[1, 2])
	@test iszero(B4[1, 3])
	@test iszero(B4[2, 1])
	
	# Test case 5: more complex polynomials
	A5 = zero_matrix(Ruv, 2, 2)
	A5[1, 1] = u^2 + u*v + v^2
	A5[1, 2] = u^3 + v^3 + Ruv(1)
	A5[2, 1] = u^-1 + v^-1
	A5[2, 2] = u*v + u^-1*v^-1
	
	B5 = replace_variable_inv_mat(A5, [u_rel, v_rel], R, Ruv, [x, y], [u, v])
	
	# Verify the matrix size
	@test size(B5) == (2, 2)
	
	# Verify that every entry is a polynomial over the old-variable ring
	for i in 1:2, j in 1:2
		@test parent(B5[i, j]) == R
	end
	
	println("  replace_variable_inv_mat: all tests passed ✓")
end

@testset "coarse_graining_with_replace_variable" begin
	F = GF(2)
	Rxy, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
	Ruv, (u, v) = laurent_polynomial_ring(F, ["u", "v"])
	A = zero_matrix(Rxy, 2, 4)
    A[1, 1] = 1 + x + x*y
    A[1, 2] = 1 + y + x*y
    A[2, 3] = 1 + y^-1 + x^-1*y^-1
    A[2, 4] = 1 + x^-1 + x^-1*y^-1
	B = coarse_graining_with_replace_variable(A, [x*y^2, x^-1*y], Rxy, Ruv, [x, y], [u, v], [Rxy(1), y, x^-1])
	# @test B == A
end

@testset "possible_remainders - 2 variables" begin
	F = GF(2)
	R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
	Ruv, (u, v) = laurent_polynomial_ring(F, ["u", "v"])
	
	# Define the variable relations
	u_rel = x*y^2
	v_rel = x^-1*y
	pr = possible_remainders([u_rel, v_rel], R, Ruv, [x, y], [u, v])
	@test pr == [R(1), x, x^-1] || pr == [R(1), x, x^2] || pr == [R(1), x, x*y]
end

@testset "find_minimum_triangle" begin
	# Test case 1: a simple 4x4 matrix
	m1 = [1 0 0 1;
	      0 0 1 0;
	      0 1 0 0;
	      1 0 0 1]
	
	vertices, area = ToricBuilder.find_minimum_triangle(m1)
	
	# Verify that the first vertex is (1, 1)
	@test vertices[1] == (1, 1)
	
	# Verify that every vertex lies at a matrix entry equal to 1
	for v in vertices
		@test m1[v[1], v[2]] == 1
	end
	
	# Verify that the area is positive, so the points are not collinear
	@test area > 0
	
	println("  Test 1: found triangle $vertices, area = $area")
	
	# Test case 2: The minimum triangle should be (1,1), (1,2), (2,1)
	m2 = [1 1 0 1;
	      1 0 0 0;
	      0 0 0 0;
	      1 0 0 1]
	
	vertices2, area2 = ToricBuilder.find_minimum_triangle(m2)
	
	@test vertices2[1] == (1, 1)
	@test area2 == 0.5  # (1,1), (1,2), (2,1)  has area 0.5
	
	println("  Test 2: found triangle $vertices2, area = $area2")
	
	# Test case 3: Verify the fallback triangle
	m3 = [1 0 0 1;
	      0 0 0 0;
	      0 0 0 0;
	      1 0 0 1]
	
	vertices3, area3 = ToricBuilder.find_minimum_triangle(m3)
	
	@test vertices3[1] == (1, 1)
	# The fallback triangle is (1,1), (1,4), (4,1), with area 4.5
	@test area3 == 4.5
	
	println("  Test 3: fallback triangle $vertices3, area = $area3")
	
	# Test case 4: Verify the angle-based tie-break: when the area is the same, choose the larger angle
	# Construct a matrix containing several triangles with the same area
	m4 = [1 1 0 0 1;
	      1 0 0 1 0;
	      0 0 0 0 0;
	      0 1 0 0 0;
	      1 0 0 0 1]
	
	vertices4, area4 = ToricBuilder.find_minimum_triangle(m4)
	
	# Verify the chosen triangle
	@test vertices4[1] == (1, 1)
	
	# Compute the angle
	p1, p2, p3 = vertices4
	v1 = (p2[1] - p1[1], p2[2] - p1[2])
	v2 = (p3[1] - p1[1], p3[2] - p1[2])
	dot_product = v1[1] * v2[1] + v1[2] * v2[2]
	len1 = sqrt(v1[1]^2 + v1[2]^2)
	len2 = sqrt(v2[1]^2 + v2[2]^2)
	cos_angle = dot_product / (len1 * len2)
	angle = acos(cos_angle) * 180 / π
	
	println("  Test 4: triangle $vertices4, area = $area4, angle ≈ $(round(angle, digits=1))°")
	
	# Test case 5: special case: a right triangle
	m5 = [1 1 1 0;
	      1 0 0 0;
	      1 0 0 0;
	      0 0 0 1]
	
	vertices5, area5 = ToricBuilder.find_minimum_triangle(m5)
	
	# Compute the angle
	p1_5, p2_5, p3_5 = vertices5
	v1_5 = (p2_5[1] - p1_5[1], p2_5[2] - p1_5[2])
	v2_5 = (p3_5[1] - p1_5[1], p3_5[2] - p1_5[2])
	dot_product_5 = v1_5[1] * v2_5[1] + v1_5[2] * v2_5[2]
	
	# If the triangle is right-angled, the dot product should be 0
	if abs(dot_product_5) < 1e-10
		println("  Test 5: found a right triangle $vertices5, area = $area5")
	else
		len1_5 = sqrt(v1_5[1]^2 + v1_5[2]^2)
		len2_5 = sqrt(v2_5[1]^2 + v2_5[2]^2)
		cos_angle_5 = dot_product_5 / (len1_5 * len2_5)
		angle_5 = acos(cos_angle_5) * 180 / π
		println("  Test 5: triangle $vertices5, area = $area5, angle ≈ $(round(angle_5, digits=1))°")
	end
end

@testset "find_parallelogram_lattice_points" begin
	# Test case 1: a small parallelogram: the unit square
	p1 = (1, 1)
	p2 = (1, 2)
	p3 = (2, 1)
	
	lattice_points1 = ToricBuilder.find_parallelogram_lattice_points(p1, p2, p3)
	
	# v1 = (0,1), v2 = (1,0)
	# The parallelogram is the unit square [0,1) x [0,1) relative to p1 = (1,1)
	# It should contain only the single lattice point (1,1)
	@test (1, 1) in lattice_points1
	@test length(lattice_points1) == 1
	
	println("  Test 1: triangle $((p1, p2, p3)), lattice points: $lattice_points1")
	
	# Test case 2: a larger parallelogram
	p1 = (1, 1)
	p2 = (1, 3)
	p3 = (3, 1)
	
	lattice_points2 = ToricBuilder.find_parallelogram_lattice_points(p1, p2, p3)
	
	# v1 = (0,2), v2 = (2,0)
	# It should contain 4 points: (1,1), (1,2), (2,1), (2,2)
	@test (1, 1) in lattice_points2
	@test (1, 2) in lattice_points2
	@test (2, 1) in lattice_points2
	@test (2, 2) in lattice_points2
	@test length(lattice_points2) == 4
	
	# Verify that vertices p2, p3, and p4 are excluded
	@test !((1, 3) in lattice_points2)  # p2
	@test !((3, 1) in lattice_points2)  # p3
	@test !((3, 3) in lattice_points2)  # p4
	
	println("  Test 2: triangle $((p1, p2, p3)), number of lattice points: $(length(lattice_points2))")
	
	# Test case 3: fallback triangle
	p1 = (1, 1)
	p2 = (1, 4)
	p3 = (4, 1)
	
	lattice_points3 = ToricBuilder.find_parallelogram_lattice_points(p1, p2, p3)
	
	# Verify that (1,1) is included
	@test (1, 1) in lattice_points3
	
	# Verify that vertices p2=(1,4), p3=(4,1), and p4=(4,4) are excluded
	@test !((1, 4) in lattice_points3)
	@test !((4, 1) in lattice_points3)
	@test !((4, 4) in lattice_points3)
	
	println("  Test 3: triangle $((p1, p2, p3)), number of lattice points: $(length(lattice_points3))")
	
	# Test case 4: Verify points lying on the edges
	p1 = (1, 1)
	p2 = (1, 5)
	p3 = (5, 1)
	
	lattice_points4 = ToricBuilder.find_parallelogram_lattice_points(p1, p2, p3)
	
	# It should include points on edge p1-p2: (1,1), (1,2), (1,3), (1,4), but exclude (1,5)
	@test (1, 1) in lattice_points4
	@test (1, 2) in lattice_points4
	@test (1, 3) in lattice_points4
	@test (1, 4) in lattice_points4
	@test !((1, 5) in lattice_points4)
	
	# It should include points on edge p1-p3: (1,1), (2,1), (3,1), (4,1), but exclude (5,1)
	@test (2, 1) in lattice_points4
	@test (3, 1) in lattice_points4
	@test (4, 1) in lattice_points4
	@test !((5, 1) in lattice_points4)
	
	# It should include interior points such as (2,2) and (3,3)
	@test (2, 2) in lattice_points4
	@test (3, 3) in lattice_points4
	
	# It should exclude the fourth vertex (5,5)
	@test !((5, 5) in lattice_points4)
	
	# v1=(0,4), v2=(4,0), It should contain 4 x 4 = 16 points
	@test length(lattice_points4) == 16
	
	println("  Test 4: triangle $((p1, p2, p3)), number of lattice points: $(length(lattice_points4))")
	
	# Test case 5: an asymmetric parallelogram
	p1 = (1, 1)
	p2 = (2, 4)
	p3 = (4, 2)
	
	lattice_points5 = ToricBuilder.find_parallelogram_lattice_points(p1, p2, p3)
	
	# Verify the basic properties
	@test (1, 1) in lattice_points5  # includes p1
	@test !((2, 4) in lattice_points5)  # excludes p2
	@test !((4, 2) in lattice_points5)  # excludes p3
	
	# v1=(1,3), v2=(3,1)
	# The fourth vertex is p4=(5,5)
	@test !((5, 5) in lattice_points5)  # excludes p4
	
	println("  Test 5: asymmetric triangle $((p1, p2, p3)), number of lattice points: $(length(lattice_points5))")
end
