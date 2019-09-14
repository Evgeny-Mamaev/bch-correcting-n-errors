import random
import unittest

from bch import calculate_generator_polynomial, encode, reverse_int, get_nth_bit, get_syndromes, \
    get_order_of_sigma, find_roots_of_sigma, get_error_positions, decode, \
    get_random_number_of_hamming_weight, berlekamp_massey_decode, get_hamming_weight, text_to_bits, \
    translate_message_to_bits_and_split_on_blocks_of_length_k, initiate, \
    translate_bits_to_message_and_glue_blocks_of_length_k
from finatefield import get_primitive_polynomial, get_cyclotomic_cosets, build_logarithmic_table


class BchTest(unittest.TestCase):
    n = 4
    power = 4
    k = 1
    t = 2
    d = 2 * t + 1
    p = 0.1

    def test_calculate_generator_polynomial(self):
        primitive_polynomial = get_primitive_polynomial(power=self.n, k=1)
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
        n = 4
        t = 2
        assert encode(0b111010001, 0b1101101, n, t) == 0b110110110110110

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
        n = 5
        t = 3
        primitive_polynomial = 0b100101
        received_message = 0b100101000000001
        cyclotomic_cosets = get_cyclotomic_cosets(n)
        logarithmic_table = build_logarithmic_table(n, primitive_polynomial)
        syndromes = get_syndromes(
            primitive_polynomial=primitive_polynomial,
            received_message=received_message,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithmic_table=logarithmic_table,
            n=n,
            t=t
        )
        assert syndromes == [0, 0, 29, 0, 23, 27]

    def test_berlekamp_massey_decode(self):
        n = 4
        t = 2
        primitive_polynomial = 0b10011
        received_message = 0b000011001100011
        cyclotomic_cosets = get_cyclotomic_cosets(n)
        logarithmic_table = build_logarithmic_table(n, primitive_polynomial)
        syndromes = get_syndromes(
            primitive_polynomial=primitive_polynomial,
            received_message=received_message,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithmic_table=logarithmic_table,
            n=n,
            t=t
        )
        assert berlekamp_massey_decode(
            syndromes=syndromes,
            logarithmic_table=logarithmic_table,
            n=n,
            t=t) == [1, 4, 9, 0]

    def test_find_roots_of_sigma(self):
        n = 4
        primitive_polynomial = 0b11001
        logarithmic_table = build_logarithmic_table(n, primitive_polynomial)
        assert find_roots_of_sigma([0, 11, 2, 3, -1, -1], n, logarithmic_table) == [0, 3, 9]
        pass

    def test_get_order_of_sigma(self):
        assert get_order_of_sigma([-1, 0, 0, -1, -1]) == 2
        assert get_order_of_sigma([1, 0, 0, 0, 1]) == 4
        assert get_order_of_sigma([-1, -1, -1, -1, -1]) == 0

    def test_get_error_positions(self):
        assert sorted(get_error_positions([0, 3, 9], 4)) == [1, 6, 12]

    def test_decode(self):
        n = 4
        t = 2
        primitive_polynomial = 0b10011
        received_message = 16896
        cyclotomic_cosets = get_cyclotomic_cosets(n)
        logarithmic_table = build_logarithmic_table(n, primitive_polynomial)
        assert decode(primitive_polynomial, received_message, cyclotomic_cosets, logarithmic_table, n, t) >> n * t == 0

    def test_encode_decode(self):
        n = 2
        t = 1
        cyclotomic_cosets = get_cyclotomic_cosets(n)
        primitive_polynomial = get_primitive_polynomial(power=n, k=1)
        logarithmic_table = build_logarithmic_table(n, primitive_polynomial)
        generator_polynomal = calculate_generator_polynomial(
            primitive_polynomial=primitive_polynomial,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithmic_table=logarithmic_table,
            n=n,
            t=t
        )
        count = 0
        for i in range(1000):
            error_vector = get_random_number_of_hamming_weight(2 ** n - 1, t)
            for message in range(0, 2 ** (2 ** n - 1 - n * t)):
                codeword = encode(generator_polynomal, message, n, t)
                received_message = codeword ^ error_vector
                decoded = decode(primitive_polynomial, received_message, cyclotomic_cosets, logarithmic_table, n, t)
                if decoded != message:
                    count += 1
        assert count == 0

    def test_get_random_number_of_hamming_weight(self):
        for i in range(100):
            length = random.randrange(100)
            if length == 0:
                continue
            weight = random.randrange(length)
            probe = get_random_number_of_hamming_weight(length=length, weight=weight)
            assert get_hamming_weight(num=probe) == weight

    def test_get_hamming_weight(self):
        assert get_hamming_weight(num=7) == 3
        assert get_hamming_weight(num=8) == 1
        assert get_hamming_weight(num=1) == 1

    def test_text_to_bits(self):
        assert text_to_bits("hello") == 0b0110100001100101011011000110110001101111

    def test_translate_message_to_bit_and_split_on_blocks_of_length_k(self):
        blocks = [0b1101111, 0b1011000,
                  0b0110001,
                  0b0101011, 0b0000110, 0b01101]
        k = 7
        bits = 0b0110100001100101011011000110110001101111
        assert translate_message_to_bits_and_split_on_blocks_of_length_k("hello", k) == blocks
        assert translate_bits_to_message_and_glue_blocks_of_length_k(blocks, k) == (bits, "hello")

    def test_initiate(self):
        assert initiate(1 / 3, 15) == (1, 3, 1, 2)
        assert initiate(1 / 3, 10) == (1, 3, 1, 2)
        assert initiate(1 / 10, 15) == (2, 15, 7, 4)
        assert initiate(1 / 31, 31) == (1, 31, 26, 5)
        assert initiate(1 / 100, 100) == (1, 63, 57, 6)


if __name__ == '__main__':
    unittest.main()
