import sys
import numpy as np

headers = {}
coords = {}
key = None
lines = ''
with open(sys.argv[1]) as f:
    for line in f:
        if '.' in line:
            if key:
                coords[key].append(lines)
            words = line.split()
            key = words[0] + '_' + words[6]
            v1 = (int(words[1]), float(words[5]))
            v2 = tuple(map(float, words[2:5]))
            if key not in headers:
                headers[key] = []
            headers[key].append(v1 + v2)
            if key not in coords:
                coords[key] = []
            lines = ''
        else:
            lines += line
    coords[key].append(lines)


def key_split(k):
    return map(int, k.split('_'))


def key_cmp(k1, k2):
    i1, j1 = key_split(k1)
    i2, j2 = key_split(k2)

    if i1 == i2:
        return cmp(j1, j2)
    else:
        return cmp(i1, i2)

for key in sorted(headers, cmp=key_cmp):
    f, r = key_split(key)
    h = headers[key]
    c = coords[key]
    l = len(h)

    e = np.diag([h[i][2] for i in range(l)])
    e[0][-1] = e[-1][0] = h[-1][3]
    
    fmt = '%6i' * 2 + '%10.2f' * l + '%8.3f%6i'
    for i in range(l):
        print fmt % ((f, h[i][0]) + tuple(e[i]) + (h[i][1], r))
        print c[i],
