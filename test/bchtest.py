import unittest

from bch import BCH, calculate_generator_polynomial, encode, reverse_int, get_nth_bit, decode, get_syndromes
from finatefield import get_primitive_polynomial, get_cyclotomic_cosets, build_logarithmic_table


class BchTest(unittest.TestCase):
    n = 4
    k = 1
    t = 2
    d = 2 * t + 1

    def test_init(self):
        BCH(self.n, self.k, self.t)

    def test_calculate_generator_polynomial(self):
        primitive_polynomial = get_primitive_polynomial(n=self.n, k=1)
        print("{0:b}".format(primitive_polynomial))
        cyclotomic_cosets = get_cyclotomic_cosets(n=self.n)
        for coset in cyclotomic_cosets:
            print("{0:0>{width}b}".format(coset, width=2 ** self.n))
        logarithmic_table = build_logarithmic_table(n=self.n, primitive_polynomial=primitive_polynomial)
        generator_polynomial = calculate_generator_polynomial(
            primitive_polynomial=primitive_polynomial,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithmic_table=logarithmic_table,
            n=self.n,
            t=self.t
        )
        print("{0:b}".format(generator_polynomial))

    def test_encode(self):
        assert encode(0b111010001, 0b110110100000000) == 0b110110110110110

    def test_reverse_int(self):
        a = 0b1101000000
        b = 0b0101000001
        c = 0b1101000001
        assert reverse_int(a, a.bit_length()) == 0b0000001011
        assert reverse_int(b, 10) == 0b1000001010
        assert reverse_int(c, c.bit_length()) == 0b1000001011

    def test_get_nth_bit(self):
        assert get_nth_bit(0b11001001, 0) == 1

    def test_get_syndromes(self):
        n = 4
        t = 2
        primitive_polynomial = 0b11001
        received_message = 0b100100100000001
        cyclotomic_cosets = get_cyclotomic_cosets(n)
        logarithmic_table = build_logarithmic_table(n, primitive_polynomial)
        assert get_syndromes(
            primitive_polynomial=primitive_polynomial,
            received_message=received_message,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithmic_table=logarithmic_table,
            n=n,
            t=t
        ) == [8, 1, 6, 2]

    def test_decode(self):
        assert decode(0b110110110110110, self.n) == 0b110110100000000


if __name__ == '__main__':
    unittest.main()
