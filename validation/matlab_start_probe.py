"""Independent default-start experiment; no MATLAB observations are inputs.

Scalar dSFMT recurrence/initialization adapted from MersenneTwister-Lab/dSFMT
(dSFMT.c, dSFMT-common.h, dSFMT-params19937.h). See DSFMT_LICENSE.txt.
Ziggurat construction follows the published equal-area normal algorithm.
Only the MATLAB R2025b seed-zero default stream has been cross-validated.
"""
import math

import numpy as np


class DefaultStream:
    def __init__(self):
        words = [0]
        for i in range(1, 768):
            words.append((1812433253 * (words[-1] ^ (words[-1] >> 30)) + i) & 0xffffffff)
        state = [words[i] | (words[i + 1] << 32) for i in range(0, 768, 2)]
        for i in range(382):
            state[i] = (state[i] & 0xfffffffffffff) | 0x3ff0000000000000
        parity = ((state[382] ^ 0x90014964b32f4329) & 0x3d84e1ac0dc82880)
        parity ^= (state[383] ^ 0x3b8d12ac548a7c7a) & 1
        if parity.bit_count() % 2 == 0:
            state[383] ^= 1
        self.state = state
        self.index = 382

    def bits(self):
        if self.index == 382:
            state = self.state
            left, right = state[382:]
            for i in range(191):
                offset = 2 * ((i + 117) % 191)
                a, b = state[2 * i:2 * i + 2]
                new_left = ((a << 19) ^ (right >> 32) ^ (right << 32) ^ state[offset]) & 0xffffffffffffffff
                new_right = ((b << 19) ^ (left >> 32) ^ (left << 32) ^ state[offset + 1]) & 0xffffffffffffffff
                state[2 * i] = (new_left >> 12) ^ (new_left & 0x000ffafffffffb3f) ^ a
                state[2 * i + 1] = (new_right >> 12) ^ (new_right & 0x000ffdfffc90fffd) ^ b
                left, right = new_left, new_right
            state[382:] = [left, right]
            self.index = 0
        value = self.state[self.index] & 0xfffffffffffff
        self.index += 1
        return value

    def uniform(self):
        return (self.bits() | 1) / 2**52

    def normal(self, count):
        edges = np.zeros(257)
        edges[255] = 3.654152885361009
        area = 0.004928673233974655
        edges[256] = area / math.exp(-edges[255]**2 / 2)
        for j in range(254, 0, -1):
            edges[j] = math.sqrt(-2 * math.log(
                area / edges[j + 1] + math.exp(-edges[j + 1]**2 / 2)))
        result = []
        while len(result) < count:
            j = ((self.bits() >> 24) & 255) + 1
            # Core uses the raw mantissa; open-interval uniforms are used
            # for rejection/tail draws, including logarithms below.
            value = (2 * (self.bits() / 2**52) - 1) * edges[j]
            if abs(value) < edges[j - 1]:
                result.append(value)
            elif j < 256:
                low = math.exp(-edges[j]**2 / 2)
                high = math.exp(-edges[j - 1]**2 / 2)
                if math.exp(-value**2 / 2) > low + self.uniform() * (high - low):
                    result.append(value)
            else:
                while True:
                    tail = -math.log(self.uniform()) / edges[255]
                    test = -math.log(self.uniform())
                    if 2 * test > tail * tail:
                        break
                result.append(math.copysign(edges[255] + tail, value))
        return np.asarray(result)


def default_start(n):
    vector = DefaultStream().normal(n)
    return vector / np.linalg.norm(vector)
