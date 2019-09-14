import sys

from bch import translate_message_to_bits_and_split_on_blocks_of_length_k, encode, get_random_number_of_hamming_weight, \
    decode, translate_bits_to_message_and_glue_blocks_of_length_k, BCH


def transmit_envelope_through_noisy_channel(code, envelope):
    codewords = stage_encoding(code, envelope)
    distorted_codewords = stage_distorting(code, codewords)
    decoded = stage_correcting(code, distorted_codewords)
    return stage_recovering(decoded, code.k)


def simulate():
    try:
        while True:
            p = float(input("Enter the probability of error of a noisy channel p: "))
            n = int(input("Enter the desired length of a code word n: "))
            code = BCH(p, n)
            envelope = input("Enter a message you want to transmit through a noisy channel: ")
            transmit_envelope_through_noisy_channel(code, envelope)
    except KeyboardInterrupt:
        print()
        print("Shutting down...")
        sys.exit()


def stage_encoding(code, envelope):
    codewords = []
    print("Stage: encoding...")
    print("__________________")
    for block in translate_message_to_bits_and_split_on_blocks_of_length_k(envelope, code.k):
        codeword = encode(code.generator_polynomial, block, code.power, code.t)
        print(
            "Block: {0:>0{a}b}, codeword: {1:>0{b}b}".format(block, codeword, a=code.n - code.t * code.power, b=code.n))
        codewords.append(codeword)
    print()
    return codewords


def stage_distorting(code, codewords):
    print("Stage: distorting...")
    print("__________________")
    distorted_codewords = []
    for codeword in codewords:
        error_vector = get_random_number_of_hamming_weight(2 ** code.power - 1, code.t)
        print("Codeword: {0:0>{a}b}, error vector: {1:0>{a}b}, distorted codeword: {2:0>{a}b}".
              format(codeword, error_vector, codeword ^ error_vector, a=code.n))
        distorted_codewords.append(codeword ^ error_vector)
    print()
    return distorted_codewords


def stage_correcting(code, distorted_codewords):
    decoded = []
    print("Stage: correcting...")
    print("__________________")
    for codeword in distorted_codewords:
        message = decode(
            code.primitive_polynomial, codeword, code.cyclotomic_cosets, code.logarithmic_table, code.power, code.t)
        print("Codeword: {0:0>{a}b}, decoded message: {1:0>{b}b}".
              format(codeword, message, a=code.n, b=code.n - code.t * code.power))
        decoded.append(message)
    print()
    return decoded


def stage_recovering(decoded, k):
    print("Stage: recovering...")
    print("__________________")
    bits, text = translate_bits_to_message_and_glue_blocks_of_length_k(decoded, k)
    print("Bits: {0:b}".format(bits))
    print()
    print("Text: {0}".format(text))
    print()
    return text
