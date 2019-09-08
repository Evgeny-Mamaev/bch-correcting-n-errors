import random

import deprecation

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
    """
    Calculates a generator polynomial of
    a particular BCH code.
    :param primitive_polynomial: a primitive
    polynomial of a particular Galois field.
    :param cyclotomic_cosets: cyclotomic
    cosets of primitive roots of a Galois
    field.
    :param logarithmic_table: a convenient form
    of the field elements for multiplication.
                               n
    :param n: the power in GF(2 ).
    :param t: a number of errors to correct.
    :return: a generator polynomial which is
    a product of several polynomails including
    a primitive polynomial obtained from
    cyclotomic cosets.
    """
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


def encode(generator_polynomial, message, n, t):
    """
    Encodes a message by shifting
    it to the highest power and
    adding a divided by a generator
    polynomial message to it.
    :param generator_polynomial:
    a generator polynomial of a
    certain finite field GF.
    :param message: a message to encode.
                               n
    :param n: the power in GF(2 ).
    :param t: a number of errors to correct.
    :return: an encoded message.
    """
    message = message << n * t
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
                if len(bin(syndrome_polynomial_of_argument_to_power)) >= len(bin(primitive_polynomial)):
                    syndrome_polynomial_of_argument_to_power = divide_polynomials(
                        polynomial1=syndrome_polynomial_of_argument_to_power,
                        polynomial2=primitive_polynomial
                    )[1]
                syndrome = flipped_logarithmic_table[syndrome_polynomial_of_argument_to_power]
                syndromes[position - 1] = syndrome
    return syndromes


def flip_dictionary(dictionary):
    """
    Flips a dictionary: key <-> value.
    :param dictionary: a dictionary where
    to flip values and keys.
    :return: a flipped dictionary.
    """
    return dict((v, k) for k, v in dictionary.items())


@deprecation.deprecated(deprecated_in="1.0",
                        details="Use the berlekamp_massey_decode function instead")
def berlekamp_massey_decode_2(syndromes, primitive_polynomial, logarithmic_table, n, t):
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
            if (syndromes[i - j] < 0) | (sigma[j] < 0):
                continue
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


def berlekamp_massey_decode(syndromes, logarithmic_table, n, t):
    """
    Calculates an error locator polynomial using
    the Berlekamp-Massey algorithm.
     [−1]          [0]
    Λ    (x) ← 1, Λ   (x) ← 1
    L(-1) ← 0, L(0) ← 0
    ∆(1) = 1
    for k ← 1 to 2t do
                                                     [k - 1]
        ∆(k) ← sum(from i=0 to L(k-1) including) of Λ       (i) * S(k−i)
        m ← min {i ∈ {0, 1, . . . , k − 1} | L(i) = L(k−1) }
        if ∆(k) = 0 then
             [k]       [k−1]
            Λ   (x) = Λ     (x)
            L(k−1) = L(k)
        else
             [k]       [k−1]                       k−m   [m−1]
            Λ   (x) ← Λ     (x) − ∆(k) * 1/∆(m) * x   * Λ     (x)
            L(k) ← max (L(k−1) , L(m−1) + k − m)
        end if
    end for
            [2t]
    Λ(x) ← Λ    (x)
    :param syndromes: a previously calculated
    array of syndromes corresponding to the
    received distorted message.
    :param logarithmic_table: a mapping of a
    power to a corresponding element of a particular
    Galois field.
                                              n
    :param n: the power in a Galois field GF(2 )
    :param t: a number of error to be corrected.
    :return: an array representing an error locator
    polynomial.
                               2        14    2
    The result looks like 1 + α  * x + α   * x .
    The array represents the powers of x accordingly
                                  0               1
    to the position. array[0] is x , array[1] is x .
    In these cells the powers of α are stored.

    """
    flipped_logarithmic_table = flip_dictionary(logarithmic_table)
    number_of_elements = 2 ** n - 1
    capital_lambdas = [[0] * (2 * t + 1)] * (2 * t + 1)
    capital_lambdas[-1][0] = 1
    capital_lambdas[0][0] = 1
    capital_ls = [0] * number_of_elements
    deltas = [0] * (2 * t + 1)
    deltas[0] = 1
    for k in range(1, 2 * t + 1):
        for i in range(capital_ls[k - 1] + 1):
            if (k - 1 - i < len(syndromes)) and (syndromes[k - 1 - i] >= 0) and (capital_lambdas[k - 1][i] > 0):
                deltas[k] ^= logarithmic_table[
                    (flipped_logarithmic_table[capital_lambdas[k - 1][i]] + syndromes[k - 1 - i]) % number_of_elements
                    ]
        m_bag = []
        for i in range(k):
            if capital_ls[i] == capital_ls[k - 1]:
                m_bag.append(i)
        m = sorted(m_bag)[0]
        if deltas[k] == 0:
            capital_lambdas[k] = capital_lambdas[k - 1].copy()
            capital_ls[k] = capital_ls[k - 1]
        else:
            capital_lambdas[k] = capital_lambdas[k - 1].copy()
            multiplied_c_l = [0] * (len(capital_lambdas[m - 1]) + k - m)
            for i in range(len(capital_lambdas[m - 1])):
                if capital_lambdas[m - 1][i] > 0:
                    multiplied_c_l[i + k - m] = logarithmic_table[
                        (flipped_logarithmic_table[capital_lambdas[m - 1][i]] +
                         flipped_logarithmic_table[deltas[k]] +
                         number_of_elements - flipped_logarithmic_table[deltas[m]]
                         ) % number_of_elements
                        ]
            who = len(capital_lambdas[k]) - len(multiplied_c_l)
            if who <= 0:
                for i in range(who):
                    capital_lambdas[k][i].append(0)
            else:
                for i in range(who):
                    multiplied_c_l.append(0)
            for i in range(len(capital_lambdas[k])):
                capital_lambdas[k][i] = capital_lambdas[k][i] ^ multiplied_c_l[i]
            capital_ls[k] = max(capital_ls[k - 1], capital_ls[m - 1] + k - m)
    result = []
    for i in capital_lambdas[2 * t]:
        result.append(flipped_logarithmic_table[i])
    return result


def find_roots_of_sigma(sigma, n, logarithmic_table):
    """
    Tries all the elements of a field to solve an
    error polynomial equls to zero.
    :param sigma: a sigma is an array representing
    an error locator polynomial.
                               2        14    2
    The result looks like 1 + α  * x + α   * x .
    The array represents the powers of x accordingly
                                  0               1
    to the position. array[0] is x , array[1] is x .
    In these cells the powers of α are stored.
                                              n
    :param n: the power in a Galois field GF(2 )
    :param logarithmic_table: a mapping of a power
    to a corresponding element of a particular
    Galois field.
    :return: an array of roots of a sigma.
    """
    roots = []
    for candidate in range(2 ** n - 1):
        result = logarithmic_table[sigma[0]]
        for power in range(1, get_order_of_sigma(sigma) + 1):
            if sigma[power] >= 0:
                result ^= logarithmic_table[(sigma[power] + candidate * power) % (2 ** n - 1)]
        if result == 0:
            roots.append(candidate)
    return roots


def get_error_positions(roots, n):
    """
    Flips the established roots of
    an error locator polynomial to
    the positions of errors.
    :param roots: an array of roots
    of an error locator polynomial.
    :param n: the power in a Galois
              n
    field GF(2 )
    :return: an array of error
    positions.
    """
    positions = []
    for root in roots:
        positions.append((2 ** n - 1 - root) % (2 ** n - 1))
    return positions


def get_order_of_sigma(sigma):
    """
    Returns an order of an
    arror locator polynomial.
    :param sigma: an error
    locator polynomial.
    :return: an order of
    sigma as an integer.
    """
    counter = 0
    for i in reversed(sigma):
        if i > -1:
            return len(sigma) - 1 - counter
        counter += 1
    return 0


def decode(primitive_polynomial, received_message, cyclotomic_cosets, logarithmic_table, n, t):
    """
    Performs decoding of a received message.
    :param primitive_polynomial: a primitive
    polynomial over a Galois field.
    :param received_message: a message is
    to be decoded.
    :param cyclotomic_cosets: cyclotomic
    cosets of primitive roots of a Galois
    field.
    :param logarithmic_table: a mapping from
    the power of a primitive element of a Galois
    field to this element.
    :param n: the power in the Galois field
       n
    G(2 ).
    :param t: a number of errors to be corrected.
    :return: a decoded message.
    """
    syndromes = get_syndromes(primitive_polynomial, received_message, cyclotomic_cosets, logarithmic_table, n, t)
    sigma = berlekamp_massey_decode(syndromes, logarithmic_table, n, t)
    roots = find_roots_of_sigma(sigma, n, logarithmic_table)
    error_positions = get_error_positions(roots, n)
    for position in error_positions:
        received_message ^= 1 << position
    received_message >>= n * t
    return received_message


def reverse_int(number, width):
    """
    Reverses bit of an integer.
    :param number: an integer to
    reverse bits positions.
    :param width: the width of
    a number.
    :return: an integer with
    the reverted order of bits.
    """
    return int('{:0{width}b}'.format(number, width=width)[::-1], 2)


def get_nth_bit(number, n):
    """
    Returns nth bit of an
    integer
    :param number: an integer
    to get a bit from.
    :param n: a bit position.
    :return: a bit on the
    specified position.
    """
    return 1 & number >> n


def get_random_number_of_hamming_weight(length, weight):
    """
    Returns a random number of a particular Hamming
    weight.
    :param length: a length of a number.
    :param weight: a weight of a number.
    :return: the number meets the requirements of
    length and Hamming weight.
    :raises: ValueError if the weight is greater
    than the length.
    """
    if weight > length:
        raise ValueError('The weight shouldn\'t be greater'
                         ' than the length: {0} > {1}'
                         .format(weight, length))
    i = 0
    result = 0
    while True:
        if i == weight:
            return result
        shift = random.randrange(length)
        power_of_two = 1 << shift
        if power_of_two & result == power_of_two:
            continue
        result |= power_of_two
        i += 1


def get_hamming_weight(num):
    """
    Calculates the Hamming weight.
    [10001] the weight is 2,
    [10] the weight is 1.
    :param num: a number to calculate
    the weight of.
    :return: the Hamming weight of
    the binary representation of
    the given number.
    """
    length = len(bin(num)) - 2
    weight = 0
    for i in range(length):
        if num & 1 == 1:
            weight += 1
        num >>= 1
    return weight
