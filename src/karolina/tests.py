from unittest import TestCase
import unittest
from random import random, choice
from atom import Atoms, Atom
from residue import Residue
from chain import Chain
from model import Model
from trajectory import Trajectory

atom1 = Atom(line = 'ATOM      1  CA  SER A  69      -5.574  49.990   8.828  1.00  1.00           C', model = 1)
atom2 = Atom(line = 'ATOM      2  CA  LEU A  70      -2.566  48.755  10.715  1.00  1.00           C', model = 1)
atom3 = Atom(line = 'ATOM     65  CA  SER B  69      -4.897  45.106   3.961  1.00  1.00           C', model = 1)
atom4 = Atom(line = 'ATOM     66  CA  LEU B  70      -2.410  47.533   2.177  1.00  1.00           C', model = 1)
atom5 = Atom(line = 'ATOM    132  CA  GLU C   1       9.281  38.319 -37.867  2.00  0.00           C', model = 1)
atom6 = Atom(line = 'ATOM    133  CA  GLY C   2       7.380  37.719 -34.242  2.00  0.00           C', model = 1)
atom7 = Atom(line = 'ATOM    134  CA  TYR D   3       4.359  38.954 -35.520  2.00  0.00           C', model = 1)
atom8 = Atom(line = 'ATOM    135  CA  ALA D   4       1.275  36.530 -34.357  2.00  0.00           C', model = 1)
atom9 = Atom(line = 'ATOM      1  CA  SER A  69      -5.515  49.920   8.873  1.00  1.00           C', model = 2)
atom10 = Atom(line = 'ATOM      2  CA  LEU A  70      -2.505  48.681  10.757  1.00  1.00           C', model = 2)
atom11 = Atom(line = 'ATOM     65  CA  SER B  69      -4.837  45.049   3.993  1.00  1.00           C', model = 2)
atom12 = Atom(line = 'ATOM     66  CA  LEU B  70      -2.351  47.481   2.215  1.00  1.00           C', model = 2)
atom13 = Atom(line = 'ATOM    137  CA  CYS C   6       5.009  42.656 -33.655  2.00  0.00           C', model = 2)
atom14 = Atom(line = 'ATOM    138  CA  ALA C   7       1.353  40.845 -34.338  2.00  0.00           C', model = 2)
atom15 = Atom(line = 'ATOM    139  CA  THR D   8       1.287  39.008 -31.293  2.00  0.00           C', model = 2)
atom16 = Atom(line = 'ATOM    140  CA  TYR D   9       1.880  42.661 -29.443  2.00  0.00           C', model = 2)

res1 = Residue(Atoms([atom1]), 'SER', 69)
res2 = Residue(Atoms([atom2]), 'LEU', 70)
res3 = Residue(Atoms([atom3]), 'SER', 69)
res4 = Residue(Atoms([atom4]), 'LEU', 70)
res5 = Residue(Atoms([atom5]), 'GLU', 1)
res6 = Residue(Atoms([atom6]), 'GLY', 2)
res7 = Residue(Atoms([atom7]), 'TYR', 3)
res8 = Residue(Atoms([atom8]), 'ALA', 4)
res9 = Residue(Atoms([atom9]), 'SER', 69)
res10 = Residue(Atoms([atom10]), 'LEU', 70)
res11 = Residue(Atoms([atom11]), 'SER', 69)
res12 = Residue(Atoms([atom12]), 'LEU', 70)
res13 = Residue(Atoms([atom13]), 'CYS', 6)
res14 = Residue(Atoms([atom14]), 'ALA', 7)
res15 = Residue(Atoms([atom15]), 'THR', 8)
res16 = Residue(Atoms([atom16]), 'TYR', 9)

chain1 = Chain(Atoms([atom1, atom2]), 'A')
chain2 = Chain(Atoms([atom3, atom4]), 'B')
chain3 = Chain(Atoms([atom5, atom6]), 'C')
chain4 = Chain(Atoms([atom7, atom8]), 'D')
chain5 = Chain(Atoms([atom9, atom10]), 'A')
chain6 = Chain(Atoms([atom11, atom12]), 'B')
chain7 = Chain(Atoms([atom13, atom14]), 'C')
chain8 = Chain(Atoms([atom15, atom16]), 'D')

headers = [[1.95, -721.16, -63.13, 0.00, -784.29], [1.95, -748.15, -89.78, 0.00, -837.93]]
receptor_id = 'AB'
ligand_id = ['D', 'C']

model1 = Model(Atoms([atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8]), headers[0], 1, 1)
model2 = Model(Atoms([atom9, atom10, atom11, atom12, atom13, atom14, atom15, atom16]), headers[1], 2, 1)

traj1 = Trajectory(Atoms([atom1, atom2, atom3, atom4, atom5, atom6, atom7, atom8, atom9,
                          atom10, atom11, atom12, atom13, atom14, atom15, atom16]),
                   receptor_id, ligand_id, 1, headers)

class testResidue(TestCase):

    def test__init__(self):
        resnames = ['LEU', 'ALA', 'GLY', 'VAL', 'MET', 'SER', 'PHE', 'LYS', 'ARG']
        for i in range(100):
            res = Residue(Atoms(), choice(resnames), i + 1)
            assert res

    def test_get_atoms(self):
        self.assertEqual(res1.atoms, res1.get_atoms())

class testChain(TestCase):

    def test__init__(self):
        chain_ids = ['A', 'B', 'C', 'D', 'E', 'F', 'G']
        atoms = Atoms([Atom(), Atom()])
        for i in range(100):
            chain = Chain(atoms, choice(chain_ids))
            assert chain

    def test_initialize_residues(self):
        residues = chain2.initialize_residues()

        self.assertEqual(res3.atoms, residues[0].atoms)
        self.assertEqual(res3.resname, residues[0].resname)
        self.assertEqual(res3.resnum, residues[0].resnum)

        self.assertEqual(res4.atoms, residues[1].atoms)
        self.assertEqual(res4.resname, residues[1].resname)
        self.assertEqual(res4.resnum, residues[1].resnum)

    def test_get_atoms(self):
        self.assertEqual(chain2.atoms, chain2.get_atoms())

    def test_get_residues(self):
        self.assertEqual(res3.atoms, ((chain2.get_residues())[0]).atoms)
        self.assertEqual(res3.resname, ((chain2.get_residues())[0]).resname)
        self.assertEqual(res3.resnum, ((chain2.get_residues())[0]).resnum)

        self.assertEqual(res4.atoms, ((chain2.get_residues())[1]).atoms)
        self.assertEqual(res4.resname, ((chain2.get_residues())[1]).resname)
        self.assertEqual(res4.resnum, ((chain2.get_residues())[1]).resnum)

class testModel(TestCase):

    def test__init__(self):
        atoms = Atoms([Atom(), Atom(), Atom(), Atom()])
        for i in range(100):
            model = Model(atoms, [random(), random(), random(), random(), random()], i + 1, i)
            assert model

    def test__iter__(self):
        chains = [chain1, chain2, chain3, chain4]
        for i, j in zip(chains, model1):
            self.assertEqual(i.atoms, j.atoms)

    def test__getitem__(self):
        chains = [chain1, chain2, chain3, chain4]
        for i in range(len(chains)):
            self.assertEqual((chains[i]).atoms, (model1[i]).atoms)

    def test_initialize_chains(self):
        chains = model1.initialize_chains()

        self.assertEqual(chain1.atoms, chains[0].atoms)
        self.assertEqual(chain1.chid, chains[0].chid)

        self.assertEqual(chain2.atoms, chains[1].atoms)
        self.assertEqual(chain2.chid, chains[1].chid)

        self.assertEqual(chain3.atoms, chains[2].atoms)
        self.assertEqual(chain3.chid, chains[2].chid)

        self.assertEqual(chain4.atoms, chains[3].atoms)
        self.assertEqual(chain4.chid, chains[3].chid)

    def test_initialize_residues(self):
        residues = model1.initialize_residues()

        self.assertEqual(res1.atoms, residues[0].atoms)
        self.assertEqual(res1.resname, residues[0].resname)
        self.assertEqual(res1.resnum, residues[0].resnum)

        self.assertEqual(res2.atoms, residues[1].atoms)
        self.assertEqual(res2.resname, residues[1].resname)
        self.assertEqual(res2.resnum, residues[1].resnum)

        self.assertEqual(res3.atoms, residues[2].atoms)
        self.assertEqual(res3.resname, residues[2].resname)
        self.assertEqual(res3.resnum, residues[2].resnum)

        self.assertEqual(res4.atoms, residues[3].atoms)
        self.assertEqual(res4.resname, residues[3].resname)
        self.assertEqual(res4.resnum, residues[3].resnum)

        self.assertEqual(res5.atoms, residues[4].atoms)
        self.assertEqual(res5.resname, residues[4].resname)
        self.assertEqual(res5.resnum, residues[4].resnum)

        self.assertEqual(res6.atoms, residues[5].atoms)
        self.assertEqual(res6.resname, residues[5].resname)
        self.assertEqual(res6.resnum, residues[5].resnum)

        self.assertEqual(res7.atoms, residues[6].atoms)
        self.assertEqual(res7.resname, residues[6].resname)
        self.assertEqual(res7.resnum, residues[6].resnum)

        self.assertEqual(res8.atoms, residues[7].atoms)
        self.assertEqual(res8.resname, residues[7].resname)
        self.assertEqual(res8.resnum, residues[7].resnum)

    def test_get_atoms(self):
        self.assertEqual(model1.atoms, model1.get_atoms())

    def test_get_residues(self):
        residues = model1.get_residues()

        self.assertEqual(res1.atoms, residues[0].atoms)
        self.assertEqual(res1.resname, residues[0].resname)
        self.assertEqual(res1.resnum, residues[0].resnum)

        self.assertEqual(res2.atoms, residues[1].atoms)
        self.assertEqual(res2.resname, residues[1].resname)
        self.assertEqual(res2.resnum, residues[1].resnum)

        self.assertEqual(res3.atoms, residues[2].atoms)
        self.assertEqual(res3.resname, residues[2].resname)
        self.assertEqual(res3.resnum, residues[2].resnum)

        self.assertEqual(res4.atoms, residues[3].atoms)
        self.assertEqual(res4.resname, residues[3].resname)
        self.assertEqual(res4.resnum, residues[3].resnum)

        self.assertEqual(res5.atoms, residues[4].atoms)
        self.assertEqual(res5.resname, residues[4].resname)
        self.assertEqual(res5.resnum, residues[4].resnum)

        self.assertEqual(res6.atoms, residues[5].atoms)
        self.assertEqual(res6.resname, residues[5].resname)
        self.assertEqual(res6.resnum, residues[5].resnum)

        self.assertEqual(res7.atoms, residues[6].atoms)
        self.assertEqual(res7.resname, residues[6].resname)
        self.assertEqual(res7.resnum, residues[6].resnum)

        self.assertEqual(res8.atoms, residues[7].atoms)
        self.assertEqual(res8.resname, residues[7].resname)
        self.assertEqual(res8.resnum, residues[7].resnum)

    def test_get_chains(self):
        chains = model1.get_chains()

        self.assertEqual(chain1.atoms, chains[0].atoms)
        self.assertEqual(chain1.chid, chains[0].chid)

        self.assertEqual(chain2.atoms, chains[1].atoms)
        self.assertEqual(chain2.chid, chains[1].chid)

        self.assertEqual(chain3.atoms, chains[2].atoms)
        self.assertEqual(chain3.chid, chains[2].chid)

        self.assertEqual(chain4.atoms, chains[3].atoms)
        self.assertEqual(chain4.chid, chains[3].chid)

        chain = model1.get_chains(chain = 'C')

        self.assertEqual(chain3.atoms, chain[0].atoms)
        self.assertEqual(chain3.chid, chain[0].chid)

        chains = model1.get_chains(chain = 'CA')

        self.assertEqual(chain3.atoms, chains[0].atoms)
        self.assertEqual(chain3.chid, chains[0].chid)

        self.assertEqual(chain1.atoms, chains[1].atoms)
        self.assertEqual(chain1.chid, chains[1].chid)

class testTrajectory(TestCase):

    def test__init__(self):
        receptor_ids = ['ABC', 'DE', 'GQTY', 'AD', 'YRDHB', 'AK']
        ligand_ids = ['R', 'U', 'M', 'X', 'J', 'P']
        atoms = Atoms([Atom(), Atom(), Atom(), Atom(), Atom(), Atom(), Atom(), Atom()])
        for i in range(100):
            traj = Trajectory(atoms,
                              choice(receptor_ids),
                              [choice(ligand_ids), choice(ligand_ids)],
                              i + 1,
                              [[random(), random(), random(), random(), random()],
                               [random(), random(), random(), random(), random()]])
            assert traj

        atoms1 = Atoms([Atom(model = 4), Atom(model = 4)])
        atoms2 = Atoms([Atom(model = 7), Atom(model = 7)])
        models = [Model(atoms1, [random(), random(), random(), random(), random()], 4, 1),
                  Model(atoms2, [random(), random(), random(), random(), random()], 7, 1)]
        for i in range(100):
            traj = Trajectory(models,
                              choice(receptor_ids),
                              [choice(ligand_ids), choice(ligand_ids)],
                              i + 1,
                              [[random(), random(), random(), random(), random()],
                               [random(), random(), random(), random(), random()]])
            assert traj
            self.assertEqual(1, (traj.models)[0].num)
            self.assertEqual(2, (traj.models)[1].num)
            self.assertEqual(1, ((traj.models)[0].atoms)[0].model)
            self.assertEqual(1, ((traj.models)[0].atoms)[1].model)
            self.assertEqual(2, ((traj.models)[1].atoms)[0].model)
            self.assertEqual(2, ((traj.models)[1].atoms)[1].model)
            

    def test__iter__(self):
        models = [model1, model2]
        for i, j in zip(models, traj1):
            self.assertEqual(i.atoms, j.atoms)

    def test__getitem__(self):
        models = [model1, model2]
        for i in range(len(models)):
            self.assertEqual((models[i]).atoms, (traj1[i]).atoms)

    def test__len__(self):
        models = [model1, model2]
        self.assertEqual(len(models), len(traj1))

    def test_initialize_models(self):
        models = traj1.initialize_models(headers, traj1.num)

        self.assertEqual(model1.atoms, models[0].atoms)
        self.assertEqual(model2.atoms, models[1].atoms)

        self.assertEqual(model1.num, models[0].num)
        self.assertEqual(model2.num, models[1].num)

    def test_initialize_residues(self):
        residues = traj1.initialize_residues()

        self.assertEqual(res1.atoms, residues[0].atoms)
        self.assertEqual(res1.resname, residues[0].resname)
        self.assertEqual(res1.resnum, residues[0].resnum)

        self.assertEqual(res2.atoms, residues[1].atoms)
        self.assertEqual(res2.resname, residues[1].resname)
        self.assertEqual(res2.resnum, residues[1].resnum)

        self.assertEqual(res3.atoms, residues[2].atoms)
        self.assertEqual(res3.resname, residues[2].resname)
        self.assertEqual(res3.resnum, residues[2].resnum)

        self.assertEqual(res4.atoms, residues[3].atoms)
        self.assertEqual(res4.resname, residues[3].resname)
        self.assertEqual(res4.resnum, residues[3].resnum)

        self.assertEqual(res5.atoms, residues[4].atoms)
        self.assertEqual(res5.resname, residues[4].resname)
        self.assertEqual(res5.resnum, residues[4].resnum)

        self.assertEqual(res6.atoms, residues[5].atoms)
        self.assertEqual(res6.resname, residues[5].resname)
        self.assertEqual(res6.resnum, residues[5].resnum)

        self.assertEqual(res7.atoms, residues[6].atoms)
        self.assertEqual(res7.resname, residues[6].resname)
        self.assertEqual(res7.resnum, residues[6].resnum)

        self.assertEqual(res8.atoms, residues[7].atoms)
        self.assertEqual(res8.resname, residues[7].resname)
        self.assertEqual(res8.resnum, residues[7].resnum)

        self.assertEqual(res9.atoms, residues[8].atoms)
        self.assertEqual(res9.resname, residues[8].resname)
        self.assertEqual(res9.resnum, residues[8].resnum)

        self.assertEqual(res10.atoms, residues[9].atoms)
        self.assertEqual(res10.resname, residues[9].resname)
        self.assertEqual(res10.resnum, residues[9].resnum)

        self.assertEqual(res11.atoms, residues[10].atoms)
        self.assertEqual(res11.resname, residues[10].resname)
        self.assertEqual(res11.resnum, residues[10].resnum)

        self.assertEqual(res12.atoms, residues[11].atoms)
        self.assertEqual(res12.resname, residues[11].resname)
        self.assertEqual(res12.resnum, residues[11].resnum)

        self.assertEqual(res13.atoms, residues[12].atoms)
        self.assertEqual(res13.resname, residues[12].resname)
        self.assertEqual(res13.resnum, residues[12].resnum)

        self.assertEqual(res14.atoms, residues[13].atoms)
        self.assertEqual(res14.resname, residues[13].resname)
        self.assertEqual(res14.resnum, residues[13].resnum)

        self.assertEqual(res15.atoms, residues[14].atoms)
        self.assertEqual(res15.resname, residues[14].resname)
        self.assertEqual(res15.resnum, residues[14].resnum)

        self.assertEqual(res16.atoms, residues[15].atoms)
        self.assertEqual(res16.resname, residues[15].resname)
        self.assertEqual(res16.resnum, residues[15].resnum)

    def test_initialize_chains(self):
        chains = traj1.initialize_chains()

        self.assertEqual(chain1.atoms, chains[0].atoms)
        self.assertEqual(chain1.chid, chains[0].chid)

        self.assertEqual(chain2.atoms, chains[1].atoms)
        self.assertEqual(chain2.chid, chains[1].chid)

        self.assertEqual(chain3.atoms, chains[2].atoms)
        self.assertEqual(chain3.chid, chains[2].chid)

        self.assertEqual(chain4.atoms, chains[3].atoms)
        self.assertEqual(chain4.chid, chains[3].chid)

        self.assertEqual(chain5.atoms, chains[4].atoms)
        self.assertEqual(chain5.chid, chains[4].chid)

        self.assertEqual(chain6.atoms, chains[5].atoms)
        self.assertEqual(chain6.chid, chains[5].chid)

        self.assertEqual(chain7.atoms, chains[6].atoms)
        self.assertEqual(chain7.chid, chains[6].chid)

        self.assertEqual(chain8.atoms, chains[7].atoms)
        self.assertEqual(chain8.chid, chains[7].chid)

    def test_get_atoms(self):
        self.assertEqual(traj1.atoms, traj1.get_atoms())

    def test_get_residues(self):
        residues = traj1.get_residues()

        self.assertEqual(res1.atoms, residues[0].atoms)
        self.assertEqual(res1.resname, residues[0].resname)
        self.assertEqual(res1.resnum, residues[0].resnum)

        self.assertEqual(res2.atoms, residues[1].atoms)
        self.assertEqual(res2.resname, residues[1].resname)
        self.assertEqual(res2.resnum, residues[1].resnum)

        self.assertEqual(res3.atoms, residues[2].atoms)
        self.assertEqual(res3.resname, residues[2].resname)
        self.assertEqual(res3.resnum, residues[2].resnum)

        self.assertEqual(res4.atoms, residues[3].atoms)
        self.assertEqual(res4.resname, residues[3].resname)
        self.assertEqual(res4.resnum, residues[3].resnum)

        self.assertEqual(res5.atoms, residues[4].atoms)
        self.assertEqual(res5.resname, residues[4].resname)
        self.assertEqual(res5.resnum, residues[4].resnum)

        self.assertEqual(res6.atoms, residues[5].atoms)
        self.assertEqual(res6.resname, residues[5].resname)
        self.assertEqual(res6.resnum, residues[5].resnum)

        self.assertEqual(res7.atoms, residues[6].atoms)
        self.assertEqual(res7.resname, residues[6].resname)
        self.assertEqual(res7.resnum, residues[6].resnum)

        self.assertEqual(res8.atoms, residues[7].atoms)
        self.assertEqual(res8.resname, residues[7].resname)
        self.assertEqual(res8.resnum, residues[7].resnum)

        self.assertEqual(res9.atoms, residues[8].atoms)
        self.assertEqual(res9.resname, residues[8].resname)
        self.assertEqual(res9.resnum, residues[8].resnum)

        self.assertEqual(res10.atoms, residues[9].atoms)
        self.assertEqual(res10.resname, residues[9].resname)
        self.assertEqual(res10.resnum, residues[9].resnum)

        self.assertEqual(res11.atoms, residues[10].atoms)
        self.assertEqual(res11.resname, residues[10].resname)
        self.assertEqual(res11.resnum, residues[10].resnum)

        self.assertEqual(res12.atoms, residues[11].atoms)
        self.assertEqual(res12.resname, residues[11].resname)
        self.assertEqual(res12.resnum, residues[11].resnum)

        self.assertEqual(res13.atoms, residues[12].atoms)
        self.assertEqual(res13.resname, residues[12].resname)
        self.assertEqual(res13.resnum, residues[12].resnum)

        self.assertEqual(res14.atoms, residues[13].atoms)
        self.assertEqual(res14.resname, residues[13].resname)
        self.assertEqual(res14.resnum, residues[13].resnum)

        self.assertEqual(res15.atoms, residues[14].atoms)
        self.assertEqual(res15.resname, residues[14].resname)
        self.assertEqual(res15.resnum, residues[14].resnum)

        self.assertEqual(res16.atoms, residues[15].atoms)
        self.assertEqual(res16.resname, residues[15].resname)
        self.assertEqual(res16.resnum, residues[15].resnum)

    def test_get_chains(self):
        chains = traj1.get_chains()

        self.assertEqual(chain1.atoms, chains[0].atoms)
        self.assertEqual(chain1.chid, chains[0].chid)

        self.assertEqual(chain2.atoms, chains[1].atoms)
        self.assertEqual(chain2.chid, chains[1].chid)

        self.assertEqual(chain3.atoms, chains[2].atoms)
        self.assertEqual(chain3.chid, chains[2].chid)

        self.assertEqual(chain4.atoms, chains[3].atoms)
        self.assertEqual(chain4.chid, chains[3].chid)

        self.assertEqual(chain5.atoms, chains[4].atoms)
        self.assertEqual(chain5.chid, chains[4].chid)

        self.assertEqual(chain6.atoms, chains[5].atoms)
        self.assertEqual(chain6.chid, chains[5].chid)

        self.assertEqual(chain7.atoms, chains[6].atoms)
        self.assertEqual(chain7.chid, chains[6].chid)

        self.assertEqual(chain8.atoms, chains[7].atoms)
        self.assertEqual(chain8.chid, chains[7].chid)

        chains = traj1.get_chains(chain = 'B')

        self.assertEqual(chain2.atoms, chains[0].atoms)
        self.assertEqual(chain2.chid, chains[0].chid)

        self.assertEqual(chain6.atoms, chains[1].atoms)
        self.assertEqual(chain6.chid, chains[1].chid)

        chains = traj1.get_chains(chain = 'BD')

        self.assertEqual(chain2.atoms, chains[0].atoms)
        self.assertEqual(chain2.chid, chains[0].chid)

        self.assertEqual(chain6.atoms, chains[1].atoms)
        self.assertEqual(chain6.chid, chains[1].chid)

        self.assertEqual(chain4.atoms, chains[2].atoms)
        self.assertEqual(chain4.chid, chains[2].chid)

        self.assertEqual(chain8.atoms, chains[3].atoms)
        self.assertEqual(chain8.chid, chains[3].chid)

    def test_get_models(self):
        models = traj1.get_models()

        self.assertEqual(model1.atoms, models[0].atoms)
        self.assertEqual(model1.num, models[0].num)

        self.assertEqual(model2.atoms, models[1].atoms)
        self.assertEqual(model2.num, models[1].num)

if __name__ == '__main__':
    unittest.main()
