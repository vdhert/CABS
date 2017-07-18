from abc import ABCMeta
from abc import abstractmethod
from tempfile import mkstemp
from os import remove
from subprocess import check_output

from utils import aa_to_short

class AlnMeta(ABCMeta):
    def __init__(self, *args, **kwargs):
        super(AlnMeta, self).__init__(*args, **kwargs)
        try:
            self.methodname
        except AttributeError:
            self.methodname = self.__name__.lower()

class AbstractAlignMethod(object):

    __metaclass__ = AlnMeta

    @abstractmethod
    def execute(self, mers1, mers2, **kwargs):
        pass

    @classmethod
    def get_subclass_dict(cls):
        return {i.methodname: i() for i in cls.__subclasses__()}

class TrivialAlign(AbstractAlignMethod):

    methodname = 'trivial'

    def execute(self, atoms1, atoms2, **kwargs):
        return zip(atoms1.atoms, atoms2.atoms)

class BLASTpAlign(AbstractAlignMethod):

    methodname = 'blastp'

    def execute(self, atoms1, atoms2, short=False):
        get_seq = lambda stc: [aa_to_short(r.resname) for r in stc.atoms]
        tn1 = mkstemp()[1]
        tn2 = mkstemp()[1]
        with open(tn1, 'w') as f1:
            seq1 = get_seq(atoms1)
            f1.write(">tmp1\n%s" % ''.join(seq1))
        with open(tn2, 'w') as f2:
            seq2 = get_seq(atoms2)
            f2.write(">tmp2\n%s" % ''.join(seq2))
        with open(tn1) as f1:
            with open(tn2) as f2:
                task = ['-task', 'blastp-short'] if short else []
                res = check_output(['blastp', '-subject', f1.name, '-query', f2.name] + task)
        remove(tn1)
        remove(tn2)
        bhit = [i for i in map(str.strip, res.split('\n\n\n')) if i.startswith('Score')][0].split('\n\n')[1:]
        mk_ind = lambda x: int(x) - 1
        aln_res = []
        for prt in bhit:
            query, match, subject = map(str.split, prt.split('\n'))
            indqs, indqe, indss, indse = map(mk_ind, (query[1], query[-1], subject[1], subject[-1]))
            it1 = iter(atoms1.atoms[indss: indse + 1])
            it2 = iter(atoms2.atoms[indqs: indse + 1])
            for i, j in zip(query[2], subject[2]):
                m1 = it1.next() if i.isalpha() else None
                m2 = it2.next() if i.isalpha() else None
                if None in (m1, m2): continue
                aln_res.append((m1, m2))
        return aln_res