import random
import unittest

from bch import calculate_generator_polynomial, encode, reverse_int, get_nth_bit, get_syndromes, \
    get_order_of_sigma, find_roots_of_sigma, get_error_positions, decode, \
    get_random_number_of_hamming_weight, berlekamp_massey_decode, get_hamming_weight, text_to_bits, \
    translate_message_to_bits_and_split_on_blocks_of_length_k, initiate, \
    translate_bits_to_message_and_glue_blocks_of_length_k
from finitefield import get_primitive_polynomial, get_cyclotomic_cosets, build_logarithmic_table


class BchTest(unittest.TestCase):
    power = 4
    k = 1
    t = 2
    d = 2 * t + 1
    p = 0.1

    def test_calculate_generator_polynomial(self):
        primitive_polynomial = get_primitive_polynomial(power=self.power, k=1)
        cyclotomic_cosets = get_cyclotomic_cosets(power=self.power)
        logarithmic_table = build_logarithmic_table(power=self.power, primitive_polynomial=primitive_polynomial)
        generator_polynomial = calculate_generator_polynomial(
            primitive_polynomial=primitive_polynomial,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithmic_table=logarithmic_table,
            power=self.power,
            t=self.t
        )
        assert generator_polynomial == 0b111010001

    def test_encode(self):
        power = 4
        t = 2
        assert encode(
            generator_polynomial=0b111010001,
            message=0b1101101,
            power=power,
            t=t) == 0b110110110110110

    def test_reverse_int(self):
        a = 0b1101000000
        b = 0b0101000001
        c = 0b1101000001
        assert reverse_int(number=a, width=a.bit_length()) == 0b0000001011
        assert reverse_int(number=b, width=10) == 0b1000001010
        assert reverse_int(number=c, width=c.bit_length()) == 0b1000001011

    def test_get_nth_bit(self):
        assert get_nth_bit(number=0b11001001, n=0) == 1

    def test_get_syndromes(self):
        power = 5
        t = 3
        primitive_polynomial = 0b100101
        received_message = 0b100101000000001
        cyclotomic_cosets = get_cyclotomic_cosets(power=power)
        logarithmic_table = build_logarithmic_table(power=power, primitive_polynomial=primitive_polynomial)
        syndromes = get_syndromes(
            primitive_polynomial=primitive_polynomial,
            received_message=received_message,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithmic_table=logarithmic_table,
            power=power,
            t=t
        )
        assert syndromes == [0, 0, 29, 0, 23, 27]

    def test_berlekamp_massey_decode(self):
        power = 4
        t = 2
        primitive_polynomial = 0b10011
        received_message = 0b000011001100011
        cyclotomic_cosets = get_cyclotomic_cosets(power=power)
        logarithmic_table = build_logarithmic_table(power=power, primitive_polynomial=primitive_polynomial)
        syndromes = get_syndromes(
            primitive_polynomial=primitive_polynomial,
            received_message=received_message,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithmic_table=logarithmic_table,
            power=power,
            t=t
        )
        assert berlekamp_massey_decode(
            syndromes=syndromes,
            logarithmic_table=logarithmic_table,
            power=power,
            t=t) == [0, 2, 14, -1, -1]

    def test_find_roots_of_sigma(self):
        power = 4
        primitive_polynomial = 0b11001
        logarithmic_table = build_logarithmic_table(power=power, primitive_polynomial=primitive_polynomial)
        assert find_roots_of_sigma([0, 11, 2, 3, -1, -1], power, logarithmic_table) == [0, 3, 9]
        pass

    def test_get_order_of_sigma(self):
        assert get_order_of_sigma(sigma=[-1, 0, 0, -1, -1]) == 2
        assert get_order_of_sigma(sigma=[1, 0, 0, 0, 1]) == 4
        assert get_order_of_sigma(sigma=[-1, -1, -1, -1, -1]) == 0

    def test_get_error_positions(self):
        assert sorted(get_error_positions(roots=[0, 3, 9], power=4)) == [0, 6, 12]

    def test_decode(self):
        power = 4
        t = 2
        primitive_polynomial = 0b10011
        cyclotomic_cosets = get_cyclotomic_cosets(power=power)
        logarithmic_table = build_logarithmic_table(power=power, primitive_polynomial=primitive_polynomial)
        generator_polynomial = calculate_generator_polynomial(
            primitive_polynomial=primitive_polynomial,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithmic_table=logarithmic_table,
            power=power,
            t=t
        )
        count = 0
        for i in range(10):
            error_vector = get_random_number_of_hamming_weight(length=2 ** power - 1, weight=t)
            for message in range(2 ** (2 ** power - 1 - power * t)):
                codeword = encode(generator_polynomial=generator_polynomial, message=message, power=power, t=t)
                received_message = codeword ^ error_vector
                decoded = decode(
                    primitive_polynomial=primitive_polynomial,
                    received_message=received_message,
                    cyclotomic_cosets=cyclotomic_cosets,
                    logarithmic_table=logarithmic_table,
                    power=power,
                    t=t)
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
        assert translate_message_to_bits_and_split_on_blocks_of_length_k(message="hello", k=k) == blocks
        assert translate_bits_to_message_and_glue_blocks_of_length_k(blocks=blocks, k=k) == (bits, "hello")

    def test_initiate(self):
        assert initiate(p=1 / 3, n=15) == (1, 3, 1, 2)
        assert initiate(p=1 / 3, n=10) == (1, 3, 1, 2)
        assert initiate(p=1 / 10, n=15) == (2, 15, 7, 4)
        assert initiate(p=1 / 31, n=31) == (1, 31, 26, 5)
        assert initiate(p=1 / 100, n=100) == (1, 63, 57, 6)


if __name__ == '__main__':
    unittest.main()
