from string import ascii_uppercase
import numpy as np
import re
from sys import stderr
from time import time, strftime, gmtime, sleep
from pkg_resources import resource_filename


# Dictionary for conversion of secondary structure from DSSP to CABS
CABS_SS = {'C': 1, 'H': 2, 'T': 3, 'E': 4, 'c': 1, 'h': 2, 't': 3, 'e': 4}

#sidechains relative coords
SIDECNT = { 'CYS': (-0.139, -1.265, 1.619, 0.019, -0.813, 1.897),
            'GLN': (-0.095, -1.674, 2.612, 0.047, -0.886, 2.991),
            'ILE': (0.094, -1.416, 1.836, -0.105, -0.659, 2.219),
            'SER': (0.121, -1.476, 1.186, 0.223, -1.042, 1.571),
            'VAL': (0.264, -1.194, 1.531, 0.077, -0.631, 1.854),
            'MET': (-0.04, -1.446, 2.587, 0.072, -0.81, 2.883),
            'PRO': (-0.751, -1.643, 0.467, -1.016, -1.228, 0.977),
            'LYS': (0.032, -1.835, 2.989, 0.002, -0.882, 3.405),
            'THR': (0.075, -1.341, 1.398, 0.051, -0.909, 1.712),
            'PHE': (0.151, -1.256, 3.161, -0.448, -0.791, 3.286),
            'ALA': (0.253, -1.133, 0.985, 0.119, -0.763, 1.312),
            'HIS': (-0.301, -1.405, 2.801, -0.207, -0.879, 3.019),
            'GLY': (0.0, -0.111, -0.111, 0.0, -0.111, -0.111),
            'ASP': (-0.287, -1.451, 1.989, 0.396, -0.798, 2.313),
            'GLU': (-0.028, -1.774, 2.546, 0.096, -0.923, 3.016),
            'LEU': (-0.069, -1.247, 2.292, 0.002, -0.462, 2.579),
            'ARG': (-0.057, -2.522, 3.639, -0.057, -1.21, 3.986),
            'TRP': (0.558, -1.694, 3.433, -0.06, -0.574, 3.834),
            'ASN': (-0.402, -1.237, 2.111, 0.132, -0.863, 2.328),
            'TYR': (0.308, -1.387, 3.492, -0.618, -0.799, 3.634)}

# Library of 1000 random peptides with up to 50 amino acids each
RANDOM_LIGAND_LIBRARY = np.reshape(
    np.fromfile(resource_filename('cabsDock', 'data/data2.dat'), sep=' '), (1000, 50, 3)
)

# Dictionary for amino acid name conversion
AA_NAMES = {
    'A': 'ALA', 'B': 'ASX', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'J': 'XLE',
    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'O': 'HOH',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T': 'THR',
    'U': 'UNK', 'V': 'VAL', 'W': 'TRP', 'X': 'XAA', 'Y': 'TYR',
    'Z': 'GLX'
}


def aa_to_long(seq):
    """Converts short amino acid name to long."""
    s = seq.upper()
    if s in AA_NAMES:
        return AA_NAMES[s]
    else:
        raise InvalidAAName(seq, 1)


def aa_to_short(seq):
    """Converts long amino acid name to short."""
    s = seq.upper()
    for short, full in AA_NAMES.items():
        if full == s:
            return short
    else:
        raise InvalidAAName(seq, 3)


def next_letter(taken_letters):
    """Returns next available letter for new protein chain."""
    return re.sub('[' + taken_letters + ']', '', ascii_uppercase)[0]


class InvalidAAName(Exception):
    """Exception raised when invalid amino acid name is used"""
    def __init__(self, name, l):
        self.name = (name, l)

    def __str__(self):
        return '%s is not a valid %d-letter amino acid code' % self.name


class ProgressBar:
    WIDTH = 60
    FORMAT = '[%s] %.1f%%\r'
    BAR0 = ' '
    BAR1 = '='

    def __init__(self, total=100, msg='', stream=stderr, delay=0):
        self.total = total
        self.current = 0
        self.stream = stream
        if msg:
            self.stream.write(msg + '\n')
        self.start_time = time()
        sleep(delay)

    def write(self):
        percent = 1.0 * self.current / self.total
        num = int(self.WIDTH * percent)
        percent = round(100. * percent, 1)
        bar = self.BAR1 * num + self.BAR0 * (self.WIDTH - num)
        self.stream.write(self.FORMAT % (bar, percent))
        self.stream.flush()

    def update(self, state=None):
        if not state:
            self.current += 1
        elif state < 0:
            self.current = self.total
        else:
            self.current = state
        if self.current > self.total:
            self.current = self.total
        self.write()

    def done(self, show_time=True):
        self.update(-1)
        self.write()
        self.stream.write('\n')
        if show_time:
            t = gmtime(time() - self.start_time)
            self.stream.write('Done in %s\n' % strftime('%H:%M:%S', t))
        self.stream.flush()


def line_count(filename):
    i = 0
    with open(filename) as f:
        for i, l in enumerate(f, 1):
            pass
    return i


def ranges(data):
    result = []
    if not data:
        return result
    idata = iter(data)
    first = prev = next(idata)
    for following in idata:
        if following - prev == 1:
            prev = following
        else:
            result.append((first, prev + 1))
            first = prev = following
    result.append((first, prev + 1))
    return result


def kabsch(t, q, concentric=False):
    if not concentric:
        t = np.subtract(t, np.average(t, 0))
        q = np.subtract(q, np.average(q, 0))
    v, s, w = np.linalg.svd(np.dot(t.T, q))
    d = np.identity(3)
    if np.linalg.det(np.dot(w.T, v.T)) < 0:
        d[2, 2] = -1
    return np.dot(np.dot(w.T, d), v.T)


def smart_flatten(l):
    """
    Function which expands and flattens a list of integers.
    m-n -> m, m+1, ..., n
    """
    fl = []
    for i in l:
        if '-' in i:
            j = i.split('-')
            if len(j) is not 2:
                raise Exception('Invalid range syntax: ' + l)
            beg = int(j[0])
            end = int(j[1])
            if beg > end:
                raise Exception('The left index(%i) is greater than the right(%i)' % (beg, end))
            for k in range(beg, end + 1):
                fl.append(k)
        else:
            fl.append(int(i))
    return fl


def kmedoids(D, k, tmax=100):
    # determine dimensions of distance matrix D
    m, n = D.shape

    if k > n:
        raise Exception('too many medoids')
    # randomly initialize an array of k medoid indices
    M = np.arange(n)
    np.random.shuffle(M)
    M = np.sort(M[:k])

    # create a copy of the array of medoid indices
    Mnew = np.copy(M)

    # initialize a dictionary to represent clusters
    C = {}
    for t in xrange(tmax):
        # determine clusters, i. e. arrays of data indices
        J = np.argmin(D[:,M], axis=1)
        for kappa in range(k):
            C[kappa] = np.where(J==kappa)[0]
        # update cluster medoids
        for kappa in range(k):
            J = np.mean(D[np.ix_(C[kappa],C[kappa])],axis=1)
            j = np.argmin(J)
            Mnew[kappa] = C[kappa][j]
        np.sort(Mnew)
        # check for convergence
        if np.array_equal(M, Mnew):
            break
        M = np.copy(Mnew)
    else:
        # final update of cluster memberships
        J = np.argmin(D[:,M], axis=1)
        for kappa in range(k):
            C[kappa] = np.where(J==kappa)[0]

    # return results
    return M, C

def rebuildCb(vec, nms):
    """Takes vector of C alpha coords and residue names and returns vector of C beta coords."""
    vec = np.insert(vec, 0, vec[0] - (vec[2] - vec[1]), axis=0)
    vec = np.append(vec, np.array([vec[-1] + (vec[-2] - vec[-3])]), axis=0)

    nvec = np.zeros((len(vec) - 2, 3))
    for i in xrange(len(vec) - 2):
        c1, c2, c3 = vec[i: i + 3]

        rdif = c3 - c1
        r2_1 = c1 - c2
        r2_3 = c3 - c2
        rsum = r2_3 + r2_1
        z = -1 * rsum / np.linalg.norm(rsum)
        x = rdif / np.linalg.norm(rdif)
        y = np.cross(x, z)

        w = np.cross(z, np.array([0, 0, 1]))

        ph_t = np.array([1, 0, 0]), w
        cph = np.dot(*ph_t)
        sph = np.linalg.norm(np.cross(*ph_t))

        cps = np.dot(w, x)
        sps = np.linalg.norm(np.cross(w, x))

        th_t = (np.array([0, 0, 1]), z)
        cth = np.dot(*th_t)
        sth = np.linalg.norm(np.cross(*th_t))

        rot = np.matrix( [  [cps * cph - sps * cth * sph, cps * sph + sph * cth * cph, sps * sth],
                            [-sps * cph - cps * cth * sph, -sps * sph + cps * cth * cph, cps * sth],
                            [sth * sph, -sth * cph, cth]])

        cb = np.dot(np.array(SIDECNT['ALA'][:3]), rot)

        nvec[i] = cb + c2

    import imp
    pdbx = imp.load_source('pdbx', '/usr/lib/python2.7/pdb.py')
    pdbx.set_trace()