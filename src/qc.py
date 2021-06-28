import unittest
from math import pi, sin, cos, sqrt, atan2

import matplotlib.pyplot as plt
import numpy as np
from numpy import complex, arcsin, log2, floor, vdot


def padded_bin(n, k):
    return bin(k)[2:].zfill(n)


def show_graph(amplitudes, title, title_x, title_y):
    y_position = np.arange(len(amplitudes))
    plt.bar(y_position, [abs(v) for v in amplitudes.values()], align='center',
            color=[[x/255 for x in complex_to_rgb(v)] for v in amplitudes.values()])
    plt.xticks(y_position, amplitudes.keys())
    plt.xticks(rotation=90)
    plt.xlabel(title_x)
    plt.ylabel(title_y)
    plt.title(title)
    plt.tight_layout()
    plt.show()


def plot_state(state, name = ""):
    n = int(log2(len(state)))
    s = dict([(padded_bin(n, k)[::-1], state[k]) for k in range(len(state))])
    d = {k + '=' + str(int(k, 2)): # [::-1]
             complex(round(v.real, 5),round(v.imag, 5)) for k, v in s.items()}
    d = dict(sorted(d.items(), key=lambda kv: (kv[0], kv[1])))
    print(name, "amplitudes:", d)
    show_graph(d, "Quantum Simulation" + " " + name, "Outcomes", "Amplitudes")
    p = {k: round(abs(v) ** 2, 5) for k, v in d.items()}
    print(name, "probabilities:", p)
    show_graph(p, "Quantum Simulation" + " " + name, "Outcomes", "Probabilities")


def pixel_state(s, cols = 1):
    n = int(log2(len(s)))
    rows = int(len(s) / cols)
    s = {int(padded_bin(n, k)[::-1], 2):
             s[k] for k in range(len(s))}

    from PIL import Image, ImageOps
    # PIL accesses images in Cartesian co-ordinates, so it is Image[columns, rows]
    img = Image.new('RGB', (cols, rows), "black")  # create a new black image
    pixels = img.load()  # create the pixel map
    for i in range(img.size[0]):  # for every col:
        for j in range(img.size[1]):  # For every row
            rgb = [int(x) for x in complex_to_rgb(s[i*rows + j], True)]
            pixels[i, j] = (rgb[0], rgb[1], rgb[2])  # set the colour accordingly
    img = ImageOps.flip(img)
    # img.show()
    img.save("state.png")


def complex_to_rgb(c, scaled_saturation = False):
    a = c.real
    b = c.imag

    val = 100;

    hue = atan2(b, a) * 180 / pi;
    if hue < 0: hue += 360

    sat = 100
    if scaled_saturation:
        sat = sqrt(a ** 2 + b ** 2) * 100

    return hsv_to_rgb(hue, sat, val)


# https://gist.github.com/eyecatchup/9536706 Colors
def hsv_to_rgb(h, s, v):
    # Make sure our arguments stay in-range
    h = max(0, min(360, h))
    s = max(0, min(100, s))
    v = max(0, min(100, v))
    
    # We accept saturation and value arguments from 0 to 100 because that's
    # how Photoshop represents those values. Internally, however, the
    # saturation and value are calculated from a range of 0 to 1.
    # We make that conversion here.
    s /= 100
    v /= 100
    
    if s == 0:
        # Achromatic (grey)
        r = g = b = v
        return [
            round(r * 255),
            round(g * 255),
            round(b * 255)
        ]

    h /= 60 #sector 0 to 5
    i = floor(h)
    f = h - i #factorial part of h
    p = v * (1 - s)
    q = v * (1 - s * f)
    t = v * (1 - s * (1 - f))

    if i == 0:
        r = v
        g = t
        b = p
    elif i == 1:
        r = q
        g = v
        b = p
    elif i == 2:
        r = p
        g = v
        b = t
    elif i == 3:
        r = p
        g = q
        b = v
    elif i  == 4:
        r = t
        g = p
        b = v
    else: #case 5:
        r = v
        g = p
        b = q

    return [
        round(r * 255),
        round(g * 255),
        round(b * 255)
    ]


def init_state(n):
    state = [0 for k in range(2 ** n)]
    state[0] = 1
    return state


def transform(state, t, gate, cond = lambda m: True):
    n = int(log2(len(state)))
    distance = int(2 ** t)
    suffix_count = int(2 ** t)
    prefix_count = int(2 ** (n - t - 1))

    for p in range(prefix_count):
        for s in range(suffix_count):
            m0 = p * suffix_count*2 + s
            if cond(m0):
                # the 1 side of the pair
                m1 = m0 + distance
                x = state[m0]
                # print(b, "<-->", format(m1, '#0' + str(n + 2) + 'b')[2:][::-1])
                y = state[m1]
                # new amplitudes
                state[m0] = x * gate[0][0] + y * gate[0][1]
                state[m1] = x * gate[1][0] + y * gate[1][1]


def c_transform(state, c, t, gate):
    # def cond(n, m):
    #     # b = format(m, '#0' + str(n + 2) + 'b')[2:][::-1]
    #     # return b[c] == '1'
    #     return m >> c & 1
    #
    # transform(state, t, gate, cond)
    mc_transform(state, {c:1}, t, gate)

# cs: dictionary with map of qubits to be 1 or 0
def mc_transform(state, cs, t, gate):
    def cond(m):
        ret = True
        for k, v in cs.items():
            bit_set = m >> k & 1
            ret = ret and (not bit_set ^ v)
        return ret

    transform(state, t, gate, cond)


def m_transform(state, targets, gate):
    for j in targets:
        transform(state, j, gate)


def qft(state, targets):
    for j in targets:
        transform(state, j, h)
        for k in targets[:j]:
            c_transform(state, j, k, phase(pi / 2 ** (k - j)))


def iqft1(state, targets):
    l = len(targets)
    for j in targets:
        j = l - 1 - j
        transform(state, j, h) # ry(-pi/2)
        for k in targets[:j]:
            c_transform(state, j, j - 1 - k, phase(-pi / 2 ** (k + 1)))


def iqft(state, targets):
    for j in targets[::-1]:
        transform(state, j, h)
        for k in targets[:j][::-1]:
            c_transform(state, j, k, phase(-pi / 2 ** (j - k)))


# gates
h = [[1/sqrt(2), 1/sqrt(2)], [1/sqrt(2), -1/sqrt(2)]]
x = [[0, 1], [1, 0]]
z = [[1, 0], [0, -1]]


def rx(theta):
    return [[cos(theta/2), complex(0, -sin(theta/2))], [complex(0, -sin(theta/2)), cos(theta/2)]]


def ry(theta):
    return [[cos(theta/2), -sin(theta/2)], [sin(theta/2), cos(theta/2)]]


def phase(theta):
    return [[1, 0], [0, complex(cos(theta), sin(theta))]]


def param_encoding(state, targets, v):
    # theta = v * 2 * pi / 2**len(targets)
    #
    # for j in targets:
    #     transform(state, j, h)
    #
    # for i in targets:
    #     transform(state, i, phase(2 ** i * theta))
    geom_sim(state, v)

    # iqft(state, targets)
    ift_sim(state)


def fibonacci(state, targets):
    theta = 2*arcsin((sqrt(5) - 1)/2)

    n = len(targets)
    for i in range(n):
        transform(state, targets[i], ry(theta))

    for i in range(n-1):
        c_transform(state, targets[i], targets[i+1], ry(-theta))


def geom_quantum(state, v):
    targets = range(int(log2(len(state))))
    theta = v * 2 * pi / 2**len(targets)

    m_transform(state, targets, h)
    for i in targets:
        transform(state, i, phase(2 ** i * theta))


def geom_sim(state, v):
    N = len(state)
    param = cos(v*2*pi/N) + 1j*sin(v*2*pi/N)
    for k in range(N):
        state[k] = param**(k)


def ift_sim(state):
    N = len(state)
    n = int(log2(N))
    s = [state[k] for k in range(N)]
    for i in range(N):
        fourier = np.empty(N, dtype=np.complex_)
        root = cos(i*2*pi/N) + 1j*sin(i*2*pi/N)
        for k in range(N):
            fourier[k] = root**k

        state[int(padded_bin(n, i)[::-1], 2)] = vdot(fourier, s)/N


def grover_sim(state, items, iterations):
    n = int(log2(len(state)))
    for i in range(len(items)):
        items[i] = int(padded_bin(n, items[i])[::-1], 2)
    s = state.copy()

    # Grover iterate
    for _ in range(iterations):
        # oracle (reflection in bad state vector)
        for item in items:
            state[item] = -1 * state[item]

        # inversion (reflection in original state)
        inner = vdot(s, state)
        for k in range(len(state)):
            state[k] = 2*inner*s[k] - state[k]


def encode_value(state, targets, k, v):
    theta = v * 2 * pi / 2 ** len(targets)
    for i in range(len(targets)):
        d = dict([(i, int(k[i])) for i in range(len(k))])
        print(d)
        mc_transform(state, d, targets[i], phase(2 ** i * theta))


def encode_value1(state, targets, c, v):
    theta = v * 2 * pi / 2 ** len(targets)
    for i in range(len(targets)):
        c_transform(state, c, targets[i], phase(2 ** i * theta))


class TestSim(unittest.TestCase):

    def test_grover(self):
        n = 3
        s = init_state(n)
        m_transform(s, range(n), h)
        grover_sim(s, items = [3], iterations = 2)

        plot_state(s)
        pixel_state(s, 1)

    def test_dictionary(self):
        n = 5
        state = init_state(n)

        m_transform(state, [0, 1], h)

        targets = [2, 3, 4]
        m_transform(state, targets, h)

        encode_value(state, targets, '00', 1)
        encode_value(state, targets, '01', 5)
        encode_value(state, targets, '10', 2)
        encode_value(state, targets, '11', 7)

        iqft(state, targets)

        plot_state(state)
        pixel_state(state, 4)

    def test_fibonacci(self):

        n = 4
        s = init_state(n)

        fibonacci(s, range(n))

        plot_state(s)
        pixel_state(s)

    def test_param(self):
        n = 4
        s = init_state(n)
        param_encoding(s, range(n), 5.7)
        plot_state(s)
        pixel_state(s, 1)


if __name__ == "__main__":
    # unittest.main()
    TestSim.test_grover()
