import unittest

from bch import BCH
from noisychannel import transmit_envelope_through_noisy_channel


class NoisyChannelTest(unittest.TestCase):
    power = 4
    p = 0.2

    def test_stage_encoding(self):
        code = BCH(p=self.p, n=2 ** self.power - 1)
        envelope = "Both of these issues are fixed by postponing the evaluation of annotations. Instead of compiling " \
                   "code which executes expressions in annotations at their definition time, the compiler stores the " \
                   "annotation in a string form equivalent to the AST of the expression in question. If needed, " \
                   "annotations can be resolved at runtime using typing.get_type_hints(). In the common case where " \
                   "this is not required, the annotations are cheaper to store (since short strings are interned by " \
                   "the interpreter) and make startup time faster."
        assert transmit_envelope_through_noisy_channel(code=code, envelope=envelope) == envelope


if __name__ == '__main__':
    unittest.main()
