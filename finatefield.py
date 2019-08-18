import csv

import yaml


def get_primitive_polynomial(n, k):
    if n < 1 | n > 51:
        raise ValueError('The parameter n should be 1 <= n <= 51, '
                         'n is {0}'.
                         format(n))
    if k < 1 | k > 3:
        raise ValueError('The parameter k should be 1 <= k <= 3, '
                         'k is {0}'.
                         format(k))
    config = yaml.safe_load(open("config.yml"))
    dictionary = {1: 1, 2: 2, 3: 5, 4: 10}
    with open(config['primitive-polynomials'], newline='') as csvfile:
        primitive_polynomials = csv.reader(csvfile, delimiter=',')
        for row in primitive_polynomials:
            ints = list(map(lambda x: 0 if x == '' else int(x), row))
            if ints[0] == n:
                polynomial_powers = ints[dictionary[k]:dictionary[k + 1]]
                polynomial_binary = 1 << (n - 1)
                polynomial_binary |= 1
                for i in polynomial_powers:
                    polynomial_binary |= 1 << i
                return polynomial_binary


def build_logarithmic_table(n, primitive_polynomial):
    """
    Builds a logarithmic table where a
    logarithm is mapped to the binary
    representation of the corresponding
    polynomial.
    The field is generated by the
    :param primitive_polynomial.
    :param n: the power in size of the
              n
    field GF(2 ).
    :return: the logarithmic table of
    the field.
    """
    if len(bin(primitive_polynomial)) - 3 != n:
        raise ValueError('The primitive polynomial {0:b} '
                         'is not of the specified power n = {1}'.
                         format(primitive_polynomial, n))
    logarithmic_table = {-1: 0}
    for i in range(n):
        logarithmic_table[i] = 1 << i
    logarithmic_table[n] = trim_polynomial(primitive_polynomial, n)
    for i in range(n + 1, 2 ** n - 1):
        multiplied_by_x_polynomial = logarithmic_table[i - 1] << 1
        if multiplied_by_x_polynomial & (2 ** n):
            multiplied_by_x_polynomial ^= logarithmic_table[n]
        logarithmic_table[i] = trim_polynomial(multiplied_by_x_polynomial, n)
    return logarithmic_table


def multiply_polynomials(polynomial1, polynomial2):
    result = 0
    for i in range(len(bin(polynomial2)) - 2):
        if polynomial2 & (1 << i):
            result ^= polynomial1 << i
    return result


def divide_polynomials(polynomial1, polynomial2):
    quotient = 0
    reminder = polynomial1
    while len(bin(reminder)) >= len(bin(polynomial2)):
        shift = len(bin(reminder)) - len(bin(polynomial2))
        reminder ^= polynomial2 << shift
        quotient ^= 1 << shift
    return quotient, reminder


def get_cyclotomic_cosets(n):
    """
    Fills a list of cyclotomic cosets.
                     4
    For GF(16) = GF(2 ) based on the
                4
    polynomial x  + x + 1 with the
    primitive root α=x+0 cyclotomic
    cosets are:
                                     4
    m1(x) = m2(x) = m4(x) = m8(x) = x  + x + 1,
                                      4    3    2
    m3(x) = m6(x) = m9(x) = m12(x) = x  + x  + x  + x + 1,
                      2
    m5(x) = m10(x) = x  + x + 1,
                                        4    3
    m7(x) = m11(x) = m13(x) = m14(x) = x  + x  + 1.
    :param n: the power in size of the
              n
    field GF(2 ).
    :return: a list of cyclotomic cosets
    in the binary representation.
    """
    cyclotomic_sets = []
    all_cyclotomic_members = 1
    i = 0
    while all_cyclotomic_members < 2 ** (2 ** n - 2) - 1:
        cyclotomic_sets.append(0)
        k = 0
        while True:
            if not 1 & (all_cyclotomic_members >> k):
                break
            k += 1
        while True:
            k = k % (2 ** n - 1)
            if 1 & (cyclotomic_sets[i] >> k):
                break
            cyclotomic_sets[i] ^= 1 << k
            k *= 2
        all_cyclotomic_members ^= cyclotomic_sets[i]
        i += 1
    return cyclotomic_sets


def trim_polynomial(polynomial, length):
    return polynomial & ((2 ** length) - 1)


def is_power_of_two(num):
    """
    Determines the binary length of the number and
    checks whether it has only one 1 => the power of 2.
    :param num: a number to check.
    :return: true if the number has only one 1
    and all the remaining 0, false otherwise.
    """
    counter = 0
    for i in (range(len(bin(num)) - 2)):
        if num & 1:
            counter += 1
        num >>= 1
        i += 1
    return counter == 1