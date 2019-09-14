import sys

from bch import translate_message_to_bits_and_split_on_blocks_of_length_k, encode, get_random_number_of_hamming_weight, \
    decode, translate_bits_to_message_and_glue_blocks_of_length_k, BCH


def transmit_envelope_through_noisy_channel(code, envelope):
    """
    Sends an envelope through a noisy channel, encodes
    a split onto messages envelope, distorts messages,
    corrects random errors, and recovers an envelope.
    :param code: a code is used in the transmission.
    :param envelope: a string to transmit.
    :return: a recovered string.
    """
    codewords = stage_encoding(code=code, envelope=envelope)
    distorted_codewords = stage_distorting(code=code, codewords=codewords)
    decoded = stage_correcting(code=code, distorted_codewords=distorted_codewords)
    return stage_recovering(decoded=decoded, k=code.k)


def simulate():
    """
    A main method to use from the command line.
    Simulates a noisy channel and performs coding-
    decoding operations on a input string.
    """
    try:
        while True:
            p = float(input("Enter the probability of error of a noisy channel p: "))
            n = int(input("Enter the desired length of a code word n: "))
            code = BCH(p=p, n=n)
            envelope = input("Enter a message you want to transmit through a noisy channel: ")
            transmit_envelope_through_noisy_channel(code=code, envelope=envelope)
    except KeyboardInterrupt:
        print()
        print("Shutting down...")
        sys.exit()


def stage_encoding(code, envelope):
    """
    Simulates a encoding stage.
    :param code: a BCH code to use.
    :param envelope: an envelope
    to decode.
    :return: code words which are
    to be directed to the noisy
    channel (distoring stage).
    """
    codewords = []
    print("Stage: encoding...")
    print("__________________")
    for block in translate_message_to_bits_and_split_on_blocks_of_length_k(message=envelope, k=code.k):
        codeword = encode(
            generator_polynomial=code.generator_polynomial,
            message=block,
            power=code.power,
            t=code.t)
        print(
            "Block: {0:>0{a}b}, codeword: {1:>0{b}b}".format(block, codeword, a=code.n - code.t * code.power, b=code.n))
        codewords.append(codeword)
    print()
    return codewords


def stage_distorting(code, codewords):
    """
    Distorts given codewords.
    :param code: a BCH code to use.
    :param codewords: code words to
    distort with random errors of
    length of number of correctable
    errors of the code.
    :return: distorted code words.
    """
    print("Stage: distorting...")
    print("__________________")
    distorted_codewords = []
    for codeword in codewords:
        error_vector = get_random_number_of_hamming_weight(length=2 ** code.power - 1, weight=code.t)
        distorted_codeword = codeword ^ error_vector
        print("Codeword: {0:0>{a}b}, error vector: {1:0>{a}b}, distorted codeword: {2:0>{a}b}".
              format(codeword, error_vector, distorted_codeword, a=code.n))
        distorted_codewords.append(distorted_codeword)
    print()
    return distorted_codewords


def stage_correcting(code, distorted_codewords):
    """
    Corrects all the errors in the distorted
    code words.
    :param code: a BCH code to use.
    :param distorted_codewords: an array of
    distorted code words.
    :return: an array of decoded messages.
    """
    decoded = []
    print("Stage: correcting...")
    print("__________________")
    for codeword in distorted_codewords:
        message = decode(
            primitive_polynomial=code.primitive_polynomial,
            received_message=codeword,
            cyclotomic_cosets=code.cyclotomic_cosets,
            logarithmic_table=code.logarithmic_table,
            power=code.power,
            t=code.t)
        print("Codeword: {0:0>{a}b}, decoded message: {1:0>{b}b}".
              format(codeword, message, a=code.n, b=code.n - code.t * code.power))
        decoded.append(message)
    print()
    return decoded


def stage_recovering(decoded, k):
    """
    Recovers an original envelope.
    :param decoded: a list of decoded
    messages.
    :param k: length of a block of a
    split message.
    :return: a recovered integer string
    of a previously sent envelope.
    """
    print("Stage: recovering...")
    print("__________________")
    bits, text = translate_bits_to_message_and_glue_blocks_of_length_k(blocks=decoded, k=k)
    print("Bits: {0:b}".format(bits))
    print()
    print("Text: {0}".format(text))
    print()
    return text
