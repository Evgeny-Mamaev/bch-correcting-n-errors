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


def berlekamp_massey_decode(syndromes, primitive_polynomial, logarithmic_table, n, t):
    flipped_logarithmic_table = flip_dictionary(logarithmic_table)
    length = len(syndromes)
    discrepancy = -1
    l = 0
    discrepancy_last_not_null = -1
    sigma = [-1] * length
    sigma[0] = 0
    sigma_order = 0
    i = 0
    l_ro = 1
    ro_last_not_null_discrepancy = -1
    ro_last_not_null_discrepancy_previous = -1
    # due to sigma represents a polynomial with
    # coefficients if form of powers of a primitive
    # element of the field the special case is
    # considered on the 0th place:
    # the -1 means the place is empty
    # 0 means the place is occupied by 1,
    #            0
    # since alpha  = 1.
    sigma_last_not_null_discrepancy = [-1] * length
    sigma_last_not_null_discrepancy[0] = 0
    sigma_last_not_null_discrepancy_previous = [-1] * length
    sigma_last_not_null_discrepancy_previous[0] = 0
    discrepancy_last_not_null_reciprocal = 0
    while i < l + t:
        d_intermediate = 0
        for j in range(i + 1):
            if j == 0:
                d_intermediate |= 1 << syndromes[i - j]
            else:
                if sigma[j] > 0:
                    d_intermediate ^= 1 << syndromes[i - j] + sigma[j] if sigma[j] >= 0 else 0
        discrepancy = flipped_logarithmic_table[divide_polynomials(d_intermediate, primitive_polynomial)[1]]
        if discrepancy != -1:
            for k in range(length):
                shift = k - i + ro_last_not_null_discrepancy
                if (shift < 0) | (sigma_last_not_null_discrepancy_previous[shift] < 0):
                    continue
                intermediate = discrepancy + discrepancy_last_not_null_reciprocal + \
                               sigma_last_not_null_discrepancy_previous[shift]
                if intermediate == sigma[k]:
                    sigma[k] = -1
                    continue
                intermediate = flipped_logarithmic_table[divide_polynomials(
                    logarithmic_table[sigma[k]] ^ logarithmic_table[
                        intermediate % (2 ** n - 1)], primitive_polynomial)[1]]
                sigma[k] = intermediate % (2 ** n - 1)
            sigma_last_not_null_discrepancy_previous = sigma_last_not_null_discrepancy.copy()
            sigma_last_not_null_discrepancy = sigma.copy()
            discrepancy_last_not_null = discrepancy
            discrepancy_last_not_null_reciprocal = 2 ** n - 1 - (
                discrepancy_last_not_null if discrepancy_last_not_null >= 0 else 0)
            l = max(l, l_ro + i - (ro_last_not_null_discrepancy if ro_last_not_null_discrepancy >= 0 else 0))
            l_ro = get_order_of_sigma(sigma_last_not_null_discrepancy_previous)
            ro_last_not_null_discrepancy = i
        i += 1
    return sigma


def find_roots_of_sigma(sigma, n, logarithmic_table):
    roots = []
    for candidate in range(2 ** n - 1):
        result = logarithmic_table[sigma[0]]
        for power in range(1, get_order_of_sigma(sigma) + 1):
            result ^= logarithmic_table[(sigma[power] + candidate * power) % (2 ** n - 1)]
        if result == 0:
            roots.append(candidate)
    return roots


def get_error_positions(roots, n):
    positions = []
    for root in roots:
        positions.append((15 - root) % (2 ** n - 2))
    return positions


def get_order_of_sigma(sigma):
    counter = 0
    for i in reversed(sigma):
        if i > -1:
            return len(sigma) - 1 - counter
        counter += 1
    return 0


def reverse_int(number, width):
    return int('{:0{width}b}'.format(number, width=width)[::-1], 2)


def get_nth_bit(number, n):
    return 1 & number >> n
