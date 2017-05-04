"""Module to handle pdb files."""

import re
from os.path import exists
from gzip import GzipFile
from urllib2 import urlopen
from StringIO import StringIO
from atom import *
from subprocess import Popen, PIPE

__all__ = ['Pdb', 'download_pdb', 'PDBFileEmpty']


class Pdb:
    """
    Pdb parser. Initialized by:
    1. pdb filename
    2. gzipped pdb filename
    3. 4-letter pdb code
    """

    def __init__(self, *args, **kwargs):
        self.file_name = None
        self.pdb_code = None
        self.lines = None
        self.atoms = []
        self.header = []
        self.missed = []

        if args and len(args) == 1:
            if exists(args[0]):
                self.file_name = args[0]
            else:
                self.pdb_code = args[0]
        elif kwargs:
            if 'pdb_file' in kwargs:
                self.file_name = kwargs['pdb_file']
            elif 'pdb_code' in kwargs:
                self.pdb_code = kwargs['pdb_code']

        if self.file_name:
            try:
                self.lines = GzipFile(self.file_name).readlines()
            except IOError:
                self.lines = open(self.file_name).readlines()
        elif self.pdb_code:
            self.lines = download_pdb(self.pdb_code).readlines()

        if self.lines:
            current_model = 0
            for line in self.lines:
                match = re.match(r'^(ATOM|HETATM)', line)
                if match:
                    self.atoms.append(Atom(line, current_model))
                elif line[:5] == 'MODEL':
                    current_model = int(line.split()[1])
                elif line[:3] == 'END' or line[:3] == "TER":
                    pass
                elif len(self.atoms) == 0:
                    self.header.append(line)
                else:
                    self.missed.append(line)

        if not len(self.atoms):
            raise PDBFileEmpty

    def dssp(self, dssp_command='dssp'):
        """Runs dssp on the read pdb file and returns a dictionary with secondary structure"""
        proc = Popen([dssp_command, '/dev/stdin'], stdin=PIPE, stdout=PIPE)
        out = proc.communicate(input=''.join(self.lines))[0].split('\n')
        sec = {}
        p = '^([0-9 ]{5}) ([0-9 ]{4}.)([A-Z ]) ([A-Z])  ([HBEGITS ])(.*)$'
        for line in out:
            m = re.match(p, line)
            if m:
                key = m.group(2).strip() + ':' + m.group(3)
                if m.group(5) in 'HGI':
                    val = 'H'
                elif m.group(5) in 'BE':
                    val = 'E'
                elif m.group(5) in 'T':
                    val = 'T'
                else:
                    val = 'C'
                sec[key] = val
        return sec


def download_pdb(pdb_code):
    url = 'http://www.rcsb.org/pdb/files/' + pdb_code.lower() + '.pdb.gz'
    s = StringIO()
    s.write(urlopen(url).read())
    s.seek(0)
    return GzipFile(fileobj=s)


class PDBFileEmpty(Exception):
    """Exception raised when Pdb contains no atoms - usually when file other than pdb was parsed."""
    pass


if __name__ == '__main__':
    pass
