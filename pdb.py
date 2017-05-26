"""Module to handle pdb files."""

import re
from os.path import exists
from gzip import GzipFile
from urllib2 import urlopen, HTTPError
from StringIO import StringIO
from subprocess import Popen, PIPE

from atom import Atom, Atoms


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
        self.atoms = Atoms()
        self.header = []
        self.missed = {}

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
        else:
            raise Exception('Cannot create a Pdb object with no arguments!!!')

        if self.file_name:
            try:
                self.lines = GzipFile(self.file_name).readlines()
            except IOError:
                self.lines = open(self.file_name).readlines()
        elif self.pdb_code:
            self.lines = download_pdb(self.pdb_code).readlines()

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
                if current_model not in self.missed:
                    self.missed[current_model] = []
                self.missed[current_model].append(line)

        if not len(self.atoms):
            if self.file_name:
                raise PdbFileEmpty(self.file_name)
            elif self.pdb_code:
                raise PdbFileEmpty(self.pdb_code + '.pdb')

    def __repr__(self):
        return ''.join(self.lines)

    def dssp(self, dssp_command):
        """Runs dssp on the read pdb file and returns a dictionary with secondary structure"""
        try:
            proc = Popen([dssp_command, '/dev/stdin'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        except OSError:
            raise Exception('Dssp not found!!!')
        out, err = proc.communicate(input=''.join(self.lines))
        if err:
            return None
        sec = {}
        p = '^([0-9 ]{5}) ([0-9 ]{4}.)([A-Z ]) ([A-Z])  ([HBEGITS ])(.*)$'
        for line in out.split('\n'):
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
    try:
        s.write(urlopen(url).read())
    except HTTPError:
        raise InvalidPdbCode(pdb_code)
    s.seek(0)
    return GzipFile(fileobj=s)


class PdbFileEmpty(Exception):
    """Exception raised when Pdb contains no atoms - usually when file other than pdb was parsed."""
    def __init__(self, filename):
        self.filename = filename

    def __str__(self):
        return 'File: \'' + self.filename + '\' empty!!!'


class InvalidPdbCode(Exception):
    """Exception raised when invalid pdb_code is used with download_pdb"""
    def __init__(self, pdb_code):
        self.pdbCode = pdb_code

    def __str__(self):
        return self.pdbCode + ' is not a valid pdb code!!!'


if __name__ == '__main__':
    pass
