"""Module for handling CABS simulation"""

import os
import re
import numpy as np
import tarfile
from os.path import join, exists, isdir
from operator import attrgetter
from subprocess import Popen, PIPE, check_output
from random import randint
from threading import Thread
from pkg_resources import resource_filename

from vector3d import Vector3d
from trajectory import Trajectory


class CabsLattice:
    """
    This class represent a CABS-like lattice. It is initialized with:
    grid_spacing: distance between grid nodes, default 0.61 A
    r12: tuple with min and max allowed values for CA-CA pseudo-bond length
    r13: tuple with min and max allowed values for CA-CA-CA end distance
    """

    def __init__(self, grid_spacing=0.61, r12=(3.28, 4.27), r13=(4.1, 7.35)):
        self.grid = grid_spacing
        r12min = round((r12[0] / self.grid) ** 2)
        r12max = round((r12[1] / self.grid) ** 2)
        r13min = round((r13[0] / self.grid) ** 2)
        r13max = round((r13[1] / self.grid) ** 2)
        dim = int(r12max ** 0.5)

        self.vectors = []
        for i in range(-dim, dim + 1):
            for j in range(-dim, dim + 1):
                for k in range(-dim, dim + 1):
                    l = i * i + j * j + k * k
                    if r12min <= float(l) <= r12max:
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
            restraints.update_id(ids), config['ca_rest_weight'], config['sc_rest_weight']
        )
        ndim = max(protein_complex.chain_list.values()) + 2
        nmols = len(protein_complex.chain_list)
        nreps = config['replicas']
        inp = CabsRun.make_inp(config, nmols, CabsRun.FORCE_FIELD)
        total_lines = int(sum(1 + np.ceil((ch + 2) / 4.) for ch in protein_complex.chain_list.values())) \
            * nreps * config['mc_cycles'] * config['mc_annealing']

        cabs_dir = join(config['work_dir'], '.CABS')
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
            f.write(inp + restr + '0 0')

        run_cmd = CabsRun.build_exe(
            params=(ndim, nreps, nmols, maxres),
            src=resource_filename('cabsDock', 'data/data0.dat'),
            exe='cabs',
            build_command=config['fortran_compiler'],
            build_flags='-O2',
            destination=cabs_dir
        )

        with tarfile.open(resource_filename('cabsDock', 'data/data1.dat')) as f:
            f.extractall(cabs_dir)

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
            max_r = max(rest_count.values() + [1])
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
            max_r = max(rest_count.values() + [max_r, ])
            rest = ['%2i %3i %2i %3i %6.2f %6.2f\n' % (r.id1 + r.id2 + (r.distance, r.weight)) for r in rest]
            restr += ''.join(rest)

        return restr, max_r

    @staticmethod
    def build_exe(params, src, exe='cabs', build_command='gfortran', build_flags='', destination='.'):
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
        return '%i\n%i %i %i %i %i\n%.2f %.2f %.2f %.2f %.2f\n%.3f %.3f %.3f %.3f %.3f\n' % (
            randint(999, 10000),
            config['mc_annealing'],
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
        return Popen(self.cfg['exe'], cwd=self.cfg['cwd']).wait()

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
        return Trajectory.read_trajectory(traf, seq)


if __name__ == '__main__':
    pass
