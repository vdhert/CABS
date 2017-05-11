"""Module for handling CABS simualtion"""

import os
import re
import numpy as np
from os.path import join, exists, isdir, basename
from operator import attrgetter
from subprocess import Popen, PIPE, check_output
from glob import glob
from random import randint
from threading import Thread
from copy import deepcopy

from vector3d import Vector3d
from atom import Atom, Atoms
from utils import CABS_HOME, ProgressBar, line_count


class CabsLattice:
    """
    This class represent a CABS-like lattice. It is initialized with:
    grid_spacing: distance between grid nodes, default 0.61 A
    r12: tuple with min and max allowed values for CA-CA pseudo-bond length
    r13: tuple with min and max allowed values for CA-CA-CA end distance
    """
    def __init__(self, grid_spacing=0.61, r12=(3.28, 4.27), r13=(4.1, 7.35)):
        self.grid = grid_spacing
        r12min = round((r12[0] / self.grid)**2)
        r12max = round((r12[1] / self.grid)**2)
        r13min = round((r13[0] / self.grid)**2)
        r13max = round((r13[1] / self.grid)**2)
        dim = int(r12max**0.5)

        self.vectors = []
        for i in range(-dim, dim + 1):
            for j in range(-dim, dim + 1):
                for k in range(-dim, dim + 1):
                    l = i * i + j * j + k * k
                    if r12min <= l <= r12max:
                        self.vectors.append(Vector3d(i, j, k))

        n = len(self.vectors)
        self.good = np.zeros((n, n))
        for i in range(n):
            vi = self.vectors[i]
            for j in range(n):
                vj = self.vectors[j]
                if r13min < (vi + vj).mod2() < r13max and vi.cross(vj).mod2():
                    self.good[i, j] = 1

    def cast(self, ch):
        """
        Function that casts a single protein chain onto the lattice.
        Returns a list of tuples with (x, y, z) coordinates of CA atoms.
        """

        if len(ch.atoms) < 3:
            raise Exception('Protein chain too short!')

        prev = None
        coord = [Vector3d(
            round(ch.atoms[0].coord.x / self.grid),
            round(ch.atoms[0].coord.y / self.grid),
            round(ch.atoms[0].coord.z / self.grid)
        )]

        for atom in ch.atoms[1:]:
            #  iterate over atoms
            min_dr = 1e12
            min_i = -1

            for i, v in enumerate(self.vectors):
                #  iterate over all possible vectors

                if len(coord) > 2 and self.good[prev, i] == 0:
                    continue

                new = coord[-1] + v
                dr = (self.grid * new - atom.coord).mod2()
                if dr < min_dr:
                    min_dr = dr
                    min_i = i

            if min_i < 0:
                raise Exception('Unsolvable geometric problem!')
            else:
                coord.append(coord[-1] + self.vectors[min_i])
                prev = min_i

        coord.insert(0, coord[0] + coord[1] - coord[2])
        coord.append(coord[-1] + coord[-2] - coord[-3])

        return coord


class CabsRun(Thread):
    """
    Class representing single cabs run.
    TODO: add state() = initializing/running/done/error/interrupted
    Generalnie to cala klase warto przepisac.
    """
    LATTICE = CabsLattice()  # static object CabsLattice used to convert structures to CABS representation
    FORCE_FIELD = (4.0, 1.0, 1.0, 2.0, 0.125, -2.0, 0.375)  # parameters of the CABS force field

    def __init__(self, protein_complex, restraints, config):
        """
        Initialize CabsRun object.
        :param protein_complex: ProteinComplex object with initial conformation of the complex (many replicas)
        :param restraints: Restraints object with complete list of CA-CA and SG-SG restraints
        :param config: Dictionary from Job object running CabsRun
        """
        Thread.__init__(self)

        fchains, seq, ids = CabsRun.load_structure(protein_complex)
        restr, maxres = CabsRun.load_restraints(
            restraints.update_id(ids), config['ca_restraints_strength'], config['sg_restraints_strength']
        )
        ndim = max(protein_complex.chain_list.values()) + 2
        nmols = len(protein_complex.chain_list)
        nreps = config['replicas']
        inp = CabsRun.make_inp(config, nmols, CabsRun.FORCE_FIELD)
        total_lines = int(sum(1 + np.ceil((ch + 2) / 4.) for ch in protein_complex.chain_list.values())) \
            * nreps * config['mc_cycles'] * 20

        cabs_dir = join(config['work_dir'], 'CABS')
        if exists(cabs_dir):
            if not isdir(cabs_dir):
                raise Exception(cabs_dir + ' exists and is not a directory!!!')
            else:
                tra = join(cabs_dir, 'TRAF')
                if exists(tra):
                    os.remove(tra)

        else:
            os.mkdir(cabs_dir, 0755)

        with open(join(cabs_dir, 'FCHAINS'), 'w') as f:
            f.write(fchains)
        with open(join(cabs_dir, 'SEQ'), 'w') as f:
            f.write(seq)
        with open(join(cabs_dir, 'INP'), 'w') as f:
            f.write(inp + restr)

        run_cmd = CabsRun.build_exe(
            params=(ndim, nreps, nmols, maxres),
            src=join(CABS_HOME, 'src/cabs.f'),
            exe='cabs',
            build_command=config['fortran_compiler'][0],
            build_flags=config['fortran_compiler'][1],
            destination=cabs_dir
        )

        for f in glob(join(CABS_HOME, 'data/params/*')):
            l = join(cabs_dir, basename(f))
            if not exists(l):
                os.symlink(f, l)

        self.cfg = {
            'cwd': cabs_dir,
            'exe': run_cmd,
            'tra': total_lines
        }

    @staticmethod
    def load_structure(protein_complex):
        fchains = None
        seq = ''
        cabs_ids = {}
        for model in protein_complex.models():
            chains = model.chains()
            if not fchains:
                fchains = [''] * len(chains)
                ch = 1
                res = 1
                chid = model[0].chid
                for atom in model:
                    if atom.chid != chid:
                        res = 1
                        ch += 1
                        chid = atom.chid
                    cabs_ids[atom.resid_id()] = (ch, res)
                    res += 1
                    seq += '%5i%1s %1s%3s %1s%3i%6.2f\n' % (
                        atom.resnum,
                        atom.icode,
                        atom.alt,
                        atom.resname,
                        atom.chid,
                        int(atom.occ),
                        atom.bfac
                    )

            for i, chain in enumerate(chains):
                vectors = CabsRun.LATTICE.cast(chain)
                fchains[i] += str(len(vectors)) + '\n' + '\n'.join(
                    ['%i %i %i' % (int(v.x), int(v.y), int(v.z)) for v in vectors]
                ) + '\n'

        return ''.join(fchains), seq, cabs_ids

    @staticmethod
    def load_restraints(restraints, ca_weight=1.0, sg_weight=0.0):
        max_r = 0

        rest = [r for r in restraints.data if not r.sg]
        restr = '%i %.2f\n' % (len(rest), ca_weight)
        if ca_weight:
            rest.sort(key=attrgetter('id2'))
            rest.sort(key=attrgetter('id1'))
            all_ids = [r.id1 for r in rest] + [r.id2 for r in rest]
            rest_count = {i: all_ids.count(i) for i in all_ids}
            max_r = max(rest_count.values())
            rest = ['%2i %3i %2i %3i %6.2f %6.2f\n' % (r.id1 + r.id2 + (r.distance, r.weight)) for r in rest]
            restr += ''.join(rest)
    #  DO POPRAWY
        rest = [r for r in restraints.data if r.sg]
        restr += '%i %.2f\n' % (len(rest), sg_weight)
        if sg_weight:
            rest.sort(key=attrgetter('id2'))
            rest.sort(key=attrgetter('id1'))
            all_ids = [r.id1 for r in rest] + [r.id2 for r in rest]
            rest_count = {i: all_ids.count(i) for i in all_ids}
            max_r = max(rest_count.values() + [max_r,])
            rest = ['%2i %3i %2i %3i %6.2f %6.2f\n' % (r.id1 + r.id2 + (r.distance, r.weight)) for r in rest]
            restr += ''.join(rest)

        return restr, max_r

    @staticmethod
    def build_exe(params, src='cabs.f', exe='cabs', build_command='gfortran', build_flags='-O2', destination='.'):

        with open(src) as f:
            lines = f.read()

        names = ['NDIM', 'NREPS', 'NMOLS', 'MAXRES']
        for name, value in zip(names, params):
            lines = re.sub(name + '=\d+', name + '=%i' % value, lines)

        run_cmd = join(destination, exe)
        cmd = [build_command, '-o', run_cmd, build_flags, '-x', 'f77', '-']
        out, err = Popen(cmd, stdin=PIPE, stderr=PIPE).communicate(lines)
        if err:
            raise Exception(err)
        return run_cmd

    @staticmethod
    def make_inp(config, nmols, force_field):
        return '%i %i %i %i %i\n%.2f %.2f %.2f %.2f %.2f\n%.3f %.3f %.3f %.3f %.3f\n' % (
            randint(999, 10000),
            config['mc_cycles'],
            config['mc_steps'],
            config['replicas'],
            nmols,
            config['t_init'],
            config['t_final'],
            force_field[0],
            force_field[1],
            config['replicas_dtemp'],
            force_field[2],
            force_field[3],
            force_field[4],
            force_field[5],
            force_field[6]
        )

    def run(self):
        Popen(self.cfg['exe'], cwd=self.cfg['cwd']).wait()

    def status(self):
        traj = join(self.cfg['cwd'], 'TRAF')
        if not exists(traj):
            progress = 0.
        else:
            progress = 100. * int(check_output(['wc', '-l', traj]).split()[0]) / self.cfg['tra']

        return progress

    def get_trajectory(self):
        traf = join(self.cfg['cwd'], 'TRAF')
        seq = join(self.cfg['cwd'], 'SEQ')
        return CabsTrajectory(traf, seq)


class CabsTrajectory:
    """
    CABS trajectory compressed object. Holds simulation trajectory in large int array.
    Attributes:
        template: Atoms with empty coordinates. Use template.from_matrix(np.matrix(with coordinates)) to generate models
        headers: List of Headers objects extracted from trajectory with energies and temperatures.
        coordinates: large int array with coordinates in lattice units, arranged as:
        
        [x1, y1, z1][x2, y2, z2] ... [xn, yn, zn] - n is number of atoms in one model (all chains)
        [    v1    ][    v2    ] ... [    vn    ]
        
        [v1, v2, ... , vn][v1, v2, ... , vn] ... [v1, v2, ... , vn]
        [     frame1     ][     frame2     ] ... [     frameN     ] - N is number of frames in the trajectory
        
        [                       REPLICA 1                         ]
        [REPLICA 1][REPLICA 2] ... [REPLICA R] - r number of replicas
        
        This allows for purely integer indexing of arbitrary (sub)structures. Also consumes little memory.
        
        replicas: number of replicas
        frames: number of trajectory frames
        length: tuple with chain lengths
    """
    CABS_GRID = 0.61  # temporary set to constant

    class Header:
        """Trajectory header read from CABS output: energies and temperatures"""

        class CannotMerge(Exception):
            """Raised when trying to merge headers for different frames"""
            def __init__(self, h1, h2):
                self.msg = 'Cannot merge headers: %s and %s' % (h1, h2)

            def __str__(self):
                return self.msg

        def __init__(self, line):
            """
            Init header from trajectory header.
            Attributes:
                replica - replica number
                frame - number in the trajectory (time)
                length - tuple with lengths of all chains
                energy - np.matrix(N, 3) where N is number of chains. Each row contains: Etotal, Einteraction, Einternal
                    Einteraction - Energy between n-th chain and rest of the system.
                e_complex - tuple (Etotal, Ereceptor, Eligand, Ebinding)
                temperature - no surprise here - it's the temperature
            :param line: String - header for one chain in trajectory. 
            """
            header = line.split()
            self.frame = int(header[0])
            self.length = (int(header[1]) - 2,)
            self.energy = np.matrix(header[4:1:-1], float)
            self.temperature = float(header[5])
            self.replica = int(header[6])
            self.e_complex = None

        def __repr__(self):
            return 'Replica: %d Frame: %d Length: %s T: %.2f E: %s' % (
                self.replica,
                self.frame,
                self.length,
                self.temperature,
                str(self.energy)
            )

        def __add__(self, other):
            """Merges two headers from two chains of the same frame"""
            if self.replica != other.replica or self.frame != other.frame:
                raise CabsTrajectory.Header.CannotMerge(self, other)
            else:
                dt = self.temperature - other.temperature
                de = self.energy[0, 0] - other.energy[0, 0]
                if dt ** 2 > 1e-6 or de ** 2 > 1e-6:
                    raise Exception("Cannot merge headers with different T or E!!!")
                else:
                    h = deepcopy(self)
                    h.length += other.length
                    h.energy = np.concatenate([self.energy, other.energy])
            return h

    def __init__(self, traf, seq):
        """
        Load TRAF and SEQ from CABS and check if SEQ fits TRAF.
        :param traf: String with TRAF file location 
        :param seq: String with SEQ file location
        """
        self.template = Atoms(self.read_seq(seq))
        self.headers, self.coordinates = self.read_traf(traf)

        self.replicas = len(set(h.replica for h in self.headers))
        if len(self.headers) % self.replicas:
            raise Exception('Replicas have different sizes!!!')
        self.frames = len(self.headers) / self.replicas
        self.length = self.headers[0].length

        if any(self.length != h.length for h in self.headers):  # check if all frames have the same shape
            raise Exception('Invalid headers in %s!!!' % traf)

        self.flength = sum(self.length)

        if self.flength != len(self.template):
            # check if number of atoms in SEQ matches that in trajectory headers
            raise Exception('Different number of atoms in %s and %s!!!' % (traf, seq))

        # final test if information from headers agrees with number of coordinates
        size_test = self.replicas * self.frames * self.flength * 3
        if size_test != len(self.coordinates):
            raise Exception('Invalid number of atoms in %s!!!' % traf)

        self.f_size = self.flength * 3              # frame size
        self.r_size = self.f_size * self.frames     # replica size

    def r_index(self, r):
        """Returns index for the beginning of the r-th replica"""
        return r * self.r_size

    def f_index(self, r, f):
        """Returns index for the beginning of the f-th frame in the r-th replica"""
        return self.r_index(r) + f * self.f_size

    def a_index(self, r, f, i):
        """Returns index for the beginning of the i-th atom in the f-th frame in the r-th replica"""
        return self.f_index(r, f) + i

    @staticmethod
    def read_seq(filename):
        atoms = []
        with open(filename) as f:
            for i, line in enumerate(f):
                atoms.append(
                    Atom(
                        hetatm=False,
                        serial=i + 1,
                        name='CA',
                        alt=line[7],
                        resname=line[8:11],
                        chid=line[12],
                        resnum=int(line[1:5]),
                        icode=line[5],
                        occ=float(line[15]),
                        bfac=float(line[16:22])
                    )
                )
        return atoms

    @staticmethod
    def read_traf(filename):
        headers = []
        replicas = {}

        def save_header(h):
            headers.append(h)

        def save_coord(c, r):
            if r not in replicas:
                replicas[r] = []
            replicas[r].extend(c[3:-3])

        bar = ProgressBar(line_count(filename), msg='Loading TRAF')
        with open(filename) as f:
            current_header = None
            current_coord = []
            for index, line in enumerate(f):
                bar.update(index)
                if '.' in line:
                    header = CabsTrajectory.Header(line)
                    if not current_header:
                        current_header = header
                    else:
                        save_coord(current_coord, current_header.replica)
                        current_coord = []
                        try:
                            current_header = current_header + header
                        except CabsTrajectory.Header.CannotMerge:
                            save_header(current_header)
                            current_header = header
                else:
                    current_coord.extend(map(int, line.split()))
            save_header(current_header)
            save_coord(current_coord, current_header.replica)
        bar.done()

        headers.sort(key=lambda x: x.frame)
        headers.sort(key=lambda x: x.replica)
        replicas = [x for y in sorted(replicas) for x in replicas[y]]

        return headers, replicas

    def split_replicas(self):
        replicas = []
        for replica in range(self.replicas):
            atoms = Atoms()
            bar = ProgressBar(self.frames, msg='Processing replica: %d' % (replica + 1))
            for frame in range(self.frames):
                bar.update(frame)
                model = deepcopy(self.template)
                model.set_model_number(frame)
                index = self.f_index(replica, frame)
                coordinates = np.array(self.coordinates[index:index + self.f_size]).reshape((self.flength, 3))
                atoms.extend(model.from_matrix(np.matrix(coordinates * self.CABS_GRID)))
            replicas.append(atoms)
            bar.done()
        return replicas

if __name__ == '__main__':
    tra = CabsTrajectory('TRAF', 'SEQ')
    replicas = tra.split_replicas()
    for i, r in enumerate(replicas, 1):
        r.save_to_pdb('replica_%d.pdb' % i)
