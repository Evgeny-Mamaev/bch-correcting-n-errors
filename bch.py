from finatefield import get_primitive_polynomial, build_logarithmic_table, get_cyclotomic_cosets, multiply_polynomials, \
    get_polynomial_from_roots, divide_polynomials, get_positions_of_binary_ones, polynomial_of_argument_to_power


class BCH(object):
    def __init__(self, n, k, t):
        """
        Constructs a BCH code with the
        specified parameters.
        :param n: length of a code.
        :param k: length of a message.
        :param t: number of errors to
        be corrected.
        """
        primitive_polynomial = get_primitive_polynomial(n=n, k=1)
        logarithmic_table = build_logarithmic_table(n=n, primitive_polynomial=primitive_polynomial)
        cyclotomic_cosets = get_cyclotomic_cosets(n=n)
        self.generator_polynomial = calculate_generator_polynomial(
            primitive_polynomial=primitive_polynomial,
            cyclotomic_cosets=cyclotomic_cosets,
            logarithmic_table=logarithmic_table,
            n=n,
            t=t
        )
        print()
        print("{0:b}".format(self.generator_polynomial))


def calculate_generator_polynomial(primitive_polynomial, cyclotomic_cosets, logarithmic_table, n, t):
    generator_polynomial = primitive_polynomial
    for i in range(1, t):
        generator_polynomial = multiply_polynomials(
            polynomial1=generator_polynomial,
            polynomial2=get_polynomial_from_roots(
                roots=cyclotomic_cosets[i],
                n=n,
                logarithmic_table=logarithmic_table)
        )
    return generator_polynomial


def encode(generator_polynomial, message):
    return message ^ divide_polynomials(polynomial1=message, polynomial2=generator_polynomial)[1]


def get_syndromes(primitive_polynomial, received_message, cyclotomic_cosets, logarithmic_table, n, t):
    """
    Calculates syndromes based on a particular
    received message, in number of 2 * t.
    E.g. received message is 100100100000001
    should be divided by each of the minimal
    polynomials. The reminder is 1110 if it is
    divided by M1 = 11001. Using roots of M1
        2    4
    a, a  , a  three syndromes are obtained:
    S1 = 1110, due to x = a;
                                2
    S2 = 101010100, due to x = a for 1110;
                                        3
    S3 = 10001000100010000, due to x = a for
    1110.
                          3    2    1
    Note that 1110 means x  + x  + x .

    :param primitive_polynomial: a primitive
    polynomial, primitive in sense of the
              n
    given GF(2 ).
    :param received_message: a message which
    was received by the decoder.
    :param cyclotomic_cosets: cyclotomic cosets
    for building minimal polynomials.
    :param logarithmic_table: a convenient form
    of the field elements for multiplication.
                               n
    :param n: the power in GF(2 )
    :param t: a number of errors to correct.
    :return: a list of powers of a primitive
    element a as shortcuts for polynomials.
    """
    length = t * 2
    syndromes = [0] * length
    flipped_logarithmic_table = flip_dictionary(logarithmic_table)
    for i in cyclotomic_cosets:
        for position in get_positions_of_binary_ones(i):
            if position - 1 < length:
                syndrome_polynomial = divide_polynomials(
                    polynomial1=received_message,
                    polynomial2=get_polynomial_from_roots(
                        roots=i,
                        n=n,
                        logarithmic_table=logarithmic_table)
                )[1]
                syndrome_polynomial_of_argument_to_power = polynomial_of_argument_to_power(
                    polynomial=syndrome_polynomial,
                    power=position)
                if syndrome_polynomial_of_argument_to_power >= primitive_polynomial:
                    syndrome_polynomial_of_argument_to_power = divide_polynomials(
                        polynomial1=syndrome_polynomial_of_argument_to_power,
                        polynomial2=primitive_polynomial
                    )[1]
                syndrome = flipped_logarithmic_table[syndrome_polynomial_of_argument_to_power]
                syndromes[position - 1] = syndrome
    return syndromes


def flip_dictionary(dictionary):
    return dict((v, k) for k, v in dictionary.items())


def decode(codeword, n):
    number_elements_in_field = 2 ** n - 1
    b = 1
    c = 1
    l = 0
    m = -1
    for i in range(number_elements_in_field):
        d = 0
        for j in range(l):
            d ^= get_nth_bit(c, j + 1) & get_nth_bit(codeword, i - j - 1)
        d ^= get_nth_bit(codeword, i)
        if d == 0:
            annihilator = 2 ** (codeword.bit_length() + 1) - 1
            for k in range(i - l, i + 1):
                annihilator ^= 1 << k
            codeword ^= annihilator
        t = c
        annihilator = 0
        for k in range(m):
            annihilator |= 1 << k
        bcopy = b & annihilator
        bcopy <<= i - m
        c ^= bcopy
        if l <= i / 2:
            l = i + 1 - l
            m = i
            b = t
    return c


def reverse_int(number, width):
    return int('{:0{width}b}'.format(number, width=width)[::-1], 2)


def get_nth_bit(number, n):
    return 1 & number >> n