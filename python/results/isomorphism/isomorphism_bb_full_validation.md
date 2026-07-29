# Isomorphism BB Full Validation

- Source snapshot: pre-release paper-companion implementation; numerical rows
  were preserved unchanged during the clean-history export
- Sage version: `10.9`
- Rows checked: 28
- Rows passed: 28
- Row timeout seconds: 1800.0
- QCA check: `product`

| Row | L | Status | Total s | Period s | Coarse s | Decouple s | QCA s | Cell | phi1_inverse shape | Forward maps | Raw QCA=standard | Redefined QCA=standard | Verification method | Error |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | ok | 0.092 | 0.009 | 0.047 | 0.023 | 0.003 | [[1, 2], [2, 1]] | [6, 6] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 2 | 7 | ok | 0.097 | 0.007 | 0.041 | 0.031 | 0.007 | [[2, 3], [3, 1]] | [14, 14] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 3 | 3 | ok | 0.076 | 0.006 | 0.040 | 0.020 | 0.002 | [[1, 2], [2, 1]] | [6, 6] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 4 | 3 | ok | 0.097 | 0.007 | 0.040 | 0.030 | 0.007 | [[0, 3], [3, 0]] | [18, 18] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 5 | 7 | ok | 0.092 | 0.008 | 0.039 | 0.029 | 0.009 | [[1, 3], [3, 2]] | [14, 14] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 6 | 3 | ok | 0.095 | 0.008 | 0.041 | 0.030 | 0.007 | [[0, 3], [3, 0]] | [18, 18] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 7 | 21 | ok | 0.299 | 0.010 | 0.048 | 0.188 | 0.030 | [[2, 5], [5, 2]] | [42, 42] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 8 | 15 | ok | 0.151 | 0.008 | 0.044 | 0.068 | 0.017 | [[0, 5], [3, 3]] | [30, 30] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 9 | 31 | ok | 0.508 | 0.009 | 0.058 | 0.338 | 0.064 | [[3, 8], [5, 3]] | [62, 62] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 10 | 3 | ok | 0.088 | 0.008 | 0.045 | 0.026 | 0.004 | [[0, 3], [1, 1]] | [6, 6] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 11 | 12 | ok | 13.683 | 0.012 | 0.295 | 11.347 | 0.910 | [[0, 12], [12, 0]] | [288, 288] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 12 | 7 | ok | 0.970 | 0.009 | 0.074 | 0.624 | 0.178 | [[0, 7], [7, 0]] | [98, 98] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 13 | 63 | ok | 37.062 | 0.015 | 0.427 | 31.852 | 2.600 | [[6, 15], [15, 6]] | [378, 378] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 14 | 217 | ok | 56.807 | 0.021 | 0.594 | 49.496 | 3.612 | [[3, 13], [19, 10]] | [434, 434] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 15 | 21 | ok | 0.498 | 0.009 | 0.062 | 0.338 | 0.050 | [[0, 21], [1, 1]] | [42, 42] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 16 | 105 | ok | 7.589 | 0.011 | 0.163 | 6.016 | 0.937 | [[3, 12], [9, 1]] | [210, 210] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 17 | 63 | ok | 2.676 | 0.011 | 0.093 | 1.973 | 0.410 | [[0, 9], [7, 3]] | [126, 126] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 18 | 73 | ok | 3.548 | 0.012 | 0.112 | 2.595 | 0.487 | [[5, 7], [14, 5]] | [146, 146] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 19 | 255 | ok | 59.417 | 0.025 | 0.748 | 47.859 | 6.161 | [[0, 17], [15, 13]] | [510, 510] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 20 | 21 | ok | 131.865 | 0.022 | 2.163 | 113.138 | 8.714 | [[0, 21], [21, 0]] | [882, 882] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 21 | 42 | ok | 3.897 | 0.012 | 0.129 | 2.947 | 0.392 | [[2, 16], [6, 6]] | [168, 168] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 22 | 21 | ok | 0.358 | 0.007 | 0.049 | 0.227 | 0.042 | [[3, 6], [4, 1]] | [42, 42] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 23 | 186 | ok | 201.156 | 0.464 | 1.617 | 176.704 | 9.156 | [[26, 40], [34, 38]] | [744, 744] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 24 | 105 | ok | 7.451 | 0.012 | 0.171 | 5.766 | 0.673 | [[10, 11], [15, 6]] | [210, 210] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 25 | 63 | ok | 419.552 | 0.029 | 2.273 | 384.657 | 21.859 | [[7, 14], [35, 7]] | [882, 882] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 26 | 217 | ok | 50.715 | 0.050 | 0.582 | 43.003 | 4.742 | [[10, 19], [13, 3]] | [434, 434] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 27 | 186 | ok | 305.637 | 0.028 | 1.569 | 272.513 | 18.009 | [[14, 32], [16, 10]] | [744, 744] | False | False | True | inverse_degree_two_triangular_certificate |  |
| 28 | 217 | ok | 53.046 | 0.021 | 0.612 | 45.639 | 3.371 | [[8, 17], [17, 9]] | [434, 434] | False | False | True | inverse_degree_two_triangular_certificate |  |
