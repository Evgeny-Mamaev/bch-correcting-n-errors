import random
import unittest

from finitefield import get_primitive_polynomial, \
    trim_polynomial, \
    build_logarithmic_table, \
    multiply_polynomials, \
    divide_polynomials, \
    get_cyclotomic_cosets, \
    get_positions_of_binary_ones, \
    get_polynomial_from_roots, polynomial_of_argument_to_power


class FinateFieldTest(unittest.TestCase):
    n = 21
    k = 6
    r = 15
    t = 3
    d = 2 * t + 1

    def test_get_primitive_polynomial(self):
        assert get_primitive_polynomial(10, 3) == 0b10011100111

    def test_trim_polynomial(self):
        assert (trim_polynomial(polynomial=15, length=3)) == 7

    def test_build_logarithmic_table(self):
        print()
        primitive_polynomial = 0b10011100111
        power = 10
        number_of_elements = 2 ** power
        table = build_logarithmic_table(power=power, primitive_polynomial=primitive_polynomial)
        for i in range(number_of_elements ** 2):
            a = random.randrange(number_of_elements - 1)
            b = random.randrange(number_of_elements - 1)
            expected = divide_polynomials(
                polynomial1=multiply_polynomials(table[a], table[b]),
                polynomial2=primitive_polynomial)[1]
            actual = table[(a + b) % (number_of_elements - 1)]
            assert actual == expected

    def test_build_logarithmic_table_exception(self):
        with self.assertRaises(ValueError):
            build_logarithmic_table(power=10, primitive_polynomial=19)

    def test_multiply_polynomials(self):
        polynomial1 = 0b10
        polynomial2 = 0b11
        assert multiply_polynomials(polynomial1=polynomial1, polynomial2=polynomial1) == 0b100
        assert multiply_polynomials(polynomial1=polynomial2, polynomial2=polynomial2) == 0b101
        polynomial3 = 0b101
        polynomial4 = 0b1110
        assert multiply_polynomials(polynomial1=polynomial3, polynomial2=polynomial4) == 0b110110

    def test_divide_polynomials(self):
        polynomial1 = 0b1010011
        polynomial2 = 0b10011
        assert divide_polynomials(polynomial1=polynomial1, polynomial2=polynomial2) == (0b101, 0b1100)
        polynomial3 = 0b101
        polynomial4 = 0b11
        assert divide_polynomials(polynomial1=polynomial3, polynomial2=polynomial4) == (polynomial4, 0)
        polynomial5 = 0b110110
        polynomial6 = 0b10011
        assert divide_polynomials(polynomial1=polynomial5, polynomial2=polynomial6) == (0b11, 0b11)
        assert divide_polynomials(polynomial1=0b100100100000001, polynomial2=0b11001) == (1863, 14)
        assert multiply_polynomials(polynomial1=1863, polynomial2=0b11001) ^ 14 == 18689

    def test_polynomial_of_argument_to_power(self):
        assert polynomial_of_argument_to_power(polynomial=0b11001, power=3) == 0b1001000000001
        assert polynomial_of_argument_to_power(polynomial=0b1, power=4) == 0b1
        assert polynomial_of_argument_to_power(polynomial=0b100, power=3) == 0b1000000
        assert polynomial_of_argument_to_power(polynomial=0b100, power=1) == 0b100

    def test_get_cyclotomic_cosets(self):
        power = 8
        result = 0
        counter = 0
        for i in get_cyclotomic_cosets(power=power):
            result ^= i
            counter += 1
        result ^= 1
        assert result == 2 ** (2 ** power - 1) - 1
        assert counter >= power

    def test_get_positions_of_binary_ones(self):
        assert get_positions_of_binary_ones(0b110110) == [1, 2, 4, 5]

    def test_get_polynomial_from_roots(self):
        power = 4
        assert get_polynomial_from_roots(
            roots=get_cyclotomic_cosets(power)[1],
            power=power,
            logarithmic_table=build_logarithmic_table(power, 0b10011)) == 0b11111


if __name__ == '__main__':
    unittest.main()
