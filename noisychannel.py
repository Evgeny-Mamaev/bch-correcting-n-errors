from bch import translate_message_to_bits_and_split_on_blocks_of_length_k, encode, get_random_number_of_hamming_weight, \
    decode


def stage_encoding(code, envelope):
    codewords = []
    print("Stage: encoding...")
    print("__________________")
    for block in translate_message_to_bits_and_split_on_blocks_of_length_k(envelope, code.k):
        codeword = encode(code.generator_polynomial, block, code.power, code.t)
        print("Block: {0:b}, codeword: {1:b}".format(block, codeword))
        codewords.append(codeword)
    return codewords


def stage_distorting(code, codewords):
    print("Stage: distorting...")
    print("__________________")
    for codeword in codewords:
        error_vector = get_random_number_of_hamming_weight(2 ** code.power - 1, code.t)
        print("Codeword: {0:b}, error vector: {1:b}, distorted codeword: {2:b}".
              format(codeword, error_vector, codeword ^ error_vector))
        codeword ^= error_vector
    return codewords


def stage_decoding(code, distorted_codewords):
    decoded = []
    print("Stage: correcting...")
    print("__________________")
    for codeword in distorted_codewords:
        message = decode(
            code.primitive_polynomial, codeword, code.cyclotomic_cosets, code.logarithmic_table, code.power, code.t)
        print("Codeword: {0:b}, decoded message: {1:b}".
              format(codeword, message))
        decoded.append(message)
    return decoded
