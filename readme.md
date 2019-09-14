Refer to MacWilliams and Sloane "The theory of Error-Correcting Codes" and Berlekamp "Algebraic coding theory" 
for all the terms used here. You can find on the internet a table of primitive binary polynomials, which can 
be used as a source of polynomials as in the file primitive-polynomials.csv.

This repository contains a program for encoding and decoding an input string with a BCH code. It:
1. Builds a BCH code according to the specified parameters p - error probability in a noisy channel - and
n - a target length of a code word. The code is capable of correcting arbitrary number of errors.
2. Changes n if it doesn't fit the specified probability of error in a channel.
3. Splits a bit representation of the string onto small pieces which satisfy the built BCH code, 
therefore gets a so called message.
4. Encodes a message and gets a so called code word.
5. Emulates a noisy channel adding an arbitrary error to each code word.
6. Corrects a distorted code word and decodes it into a message.
7. Recovers a string putting all recovered messages back together.

To run the program use Python 3:
1. Clone the repo using "git clone 'repo url'".
2. Open terminal and "cd 'directory where you just have cloned the repo'".
3. Type "python3" in the prompt.
4. Type "from noisychannel import simulate".
5. Type "simulate()", it starts the program.
6. Follow instructions on the screen.