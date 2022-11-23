"""Module for handling CABS simulation"""

import re
import tarfile
import numpy as np

from os.path import join
from operator import attrgetter
from subprocess import Popen, PIPE
from random import randint
from threading import Thread
from pkg_resources import resource_filename
from collections import OrderedDict
from tempfile import mkdtemp
from time import strftime

from CABS import logger, _JUNK
from CABS.vector3d import Vector3d
from CABS.trajectory import Trajectory

_name = 'CABS'


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
                    _l = i * i + j * j + k * k
                    if r12min <= float(_l) <= r12max:
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
    """
    LATTICE = CabsLattice()  # static object CabsLattice used to convert structures to CABS representation
    FORCE_FIELD = (4.0, 1.0, 1.0, 2.0, 0.125, -2.0, 0.375)  # parameters of the CABS force field
    FORTRAN_COMMAND = 'gfortran -O2'
    CABS_DIR_FMT = '%y%m%d%H%M%S'

    def __init__(
            self, protein_complex, restraints, work_dir, replicas, replicas_dtemp, mc_annealing, mc_cycles, mc_steps,
            temperature, ca_rest_weight, sc_rest_weight, excluding_distance
    ):
        """
        Initialize CabsRun object.
        :param protein_complex: ProteinComplex object with initial conformation of the complex (many replicas)
        :param restraints: Restraints object with complete list of CA-CA and SG-SG restraints
        """

        Thread.__init__(self)
        logger.debug(module_name=_name, msg="Loading structures...")
        fchains, seq, ids = CabsRun.load_structure(protein_complex)
        logger.debug(module_name=_name, msg="Loading restraints...")
        restr, maxres = CabsRun.load_restraints(restraints.update_id(ids), ca_rest_weight, sc_rest_weight)
        exclude = CabsRun.load_excluding(protein_complex.protein.exclude, excluding_distance, ids)

        ndim = max(protein_complex.chain_list.values()) + 2
        nmols = len(protein_complex.chain_list)
        nreps = replicas
        inp = CabsRun.make_inp(
            nmols=nmols, force_field=CabsRun.FORCE_FIELD, replicas=replicas, replicas_dtemp=replicas_dtemp,
            mc_annealing=mc_annealing, mc_cycles=mc_cycles, mc_steps=mc_steps, temperature=temperature
        )

        cabs_dir = mkdtemp(
            prefix='.' + strftime(self.CABS_DIR_FMT),
            dir=work_dir
        )

        _JUNK.append(cabs_dir)

        with open(join(cabs_dir, 'FCHAINS'), 'w') as f:
            f.write(fchains)
        with open(join(cabs_dir, 'SEQ'), 'w') as f:
            f.write(seq)
        with open(join(cabs_dir, 'INP'), 'w') as f:
            f.write(inp + restr + exclude)

        logger.debug(module_name=_name, msg="Building exe...")
        run_cmd = CabsRun.build_exe(
            params=(ndim, nreps, nmols, maxres),
            src=resource_filename('CABS', 'data/data0.dat'),
            exe='cabs',
            build_command=self.FORTRAN_COMMAND,
            destination=cabs_dir
        )

        with tarfile.open(resource_filename('CABS', 'data/data1.dat')) as f:
            
            import os
            
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
                
            
            safe_extract(f, cabs_dir)

        self.cfg = {
            'cwd': cabs_dir,
            'exe': run_cmd
        }

    @staticmethod
    def load_structure(protein_complex):
        fchains = None
        seq = ''
        cabs_ids = OrderedDict()
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
    def load_restraints(restraints, ca_weight=1.0, sg_weight=1.0):
        max_r = 0

        rest = [r for r in restraints.data if not r.sg]
        restr = '%i %.2f\n' % (len(rest), ca_weight)
        if ca_weight and len(rest):
            rest.sort(key=attrgetter('id2'))
            rest.sort(key=attrgetter('id1'))
            all_ids = [r.id1 for r in rest] + [r.id2 for r in rest]
            rest_count = {i: all_ids.count(i) for i in all_ids}
            max_r = max(1, *rest_count.values())
            rest = ['%2i %3i %2i %3i %6.2f %6.2f\n' % (
                r.id1[0], r.id1[1], r.id2[0], r.id2[1], r.distance, r.weight
            ) for r in rest]
            restr += ''.join(rest)

        rest = [r for r in restraints.data if r.sg]
        restr += '%i %.2f\n' % (len(rest), sg_weight)
        if sg_weight and len(rest):
            rest.sort(key=attrgetter('id2'))
            rest.sort(key=attrgetter('id1'))
            all_ids = [r.id1 for r in rest] + [r.id2 for r in rest]
            rest_count = {i: all_ids.count(i) for i in all_ids}
            max_r = max(max_r, *rest_count.values())
            rest = ['%2i %3i %2i %3i %6.2f %6.2f\n' % (
                r.id1[0], r.id1[1], r.id2[0], r.id2[1], r.distance, r.weight
            ) for r in rest]
            restr += ''.join(rest)

        return restr, max_r

    @staticmethod
    def load_excluding(excl, dist, cabs_ids):
        return '%i %f\n' % (len(excl), dist) + '\n'.join(
            '%i %i %i %i' % (cabs_ids[i][0], cabs_ids[i][1], cabs_ids[j][0], cabs_ids[j][1]) for i, j in excl
        )

    @staticmethod
    def build_exe(params, src, exe='cabs', build_command=FORTRAN_COMMAND, destination='.'):
        with open(src) as f:
            lines = f.read()

        names = ['NDIM', 'NREPS', 'NMOLS', 'MAXRES']
        for name, value in zip(names, params):
            lines = re.sub(name + '=\d+', name + '=%i' % value, lines)

        run_cmd = join(destination, exe)
        cmd = build_command.split() + ['-o', run_cmd, '-x', 'f77', '-']
        try:
            out, err = Popen(cmd, stdin=PIPE, stderr=PIPE).communicate(lines)
            if err:
                raise Exception(err)
        except Exception as e:
            logger.critical(_name, msg=e.message)
        return run_cmd

    @staticmethod
    def make_inp(nmols, force_field, temperature, mc_annealing, mc_cycles, mc_steps, replicas, replicas_dtemp):
        return '%i\n%i %i %i %i %i\n%.2f %.2f %.2f %.2f %.2f\n%.3f %.3f %.3f %.3f %.3f\n' % (
            randint(999, 10000),
            mc_annealing,
            mc_cycles,
            mc_steps,
            replicas,
            nmols,
            temperature[0],
            temperature[1],
            force_field[0],
            force_field[1],
            replicas_dtemp,
            force_field[2],
            force_field[3],
            force_field[4],
            force_field[5],
            force_field[6]
        )

    def run(self):
        monitor = logger.CabsObserver(interval=0.2, progress_file=join(self.cfg['cwd'], 'PROGRESS'))
        cabs_proc = Popen(self.cfg['exe'], cwd=self.cfg['cwd'], stderr=PIPE, stdin=PIPE)
        (stdout, stderr) = cabs_proc.communicate()
        if stderr:
            logger.warning(module_name=_name, msg=stderr)
        monitor.exit()

    def get_trajectory(self):
        traf = join(self.cfg['cwd'], 'TRAF')
        seq = join(self.cfg['cwd'], 'SEQ')
        return Trajectory.read_trajectory(traf, seq)


if __name__ == '__main__':
    pass
