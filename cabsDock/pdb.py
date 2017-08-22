"""Module to handle pdb files."""

import re
import os
from os.path import exists, expanduser
from gzip import GzipFile
from urllib2 import urlopen, HTTPError, URLError
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
        self.selection = ""
        self.remove_alternative_locations = True
        self.atoms = Atoms()
        self.header = []
        self.missed = {}


        if args and len(args) == 1:
            if os.path.isfile(re.split(":",args[0])[0]):
                self.file_name = args[0]
            else:
                self.pdb_code = args[0]
        if kwargs:
            if 'pdb_file' in kwargs:
                self.file_name = kwargs['pdb_file']
            elif 'pdb_code' in kwargs:
                self.pdb_code = kwargs['pdb_code']
            if 'selection' in kwargs:
                self.selection += kwargs['selection']
            if 'remove_alternative_locations' in kwargs:
                self.remove_alternative_locations = kwargs['remove_alternative_locations']
        if not self.file_name and not self.pdb_code:
            raise Exception('Cannot create a Pdb object with no arguments!!!')

        if self.file_name:
            m = re.match(r'[^:]*:([A-Z]*)', self.file_name)
            self.file_name = re.split(":",self.file_name)[0]
            try:
                self.lines = GzipFile(self.file_name).readlines()
            except IOError:
                self.lines = open(self.file_name).readlines()
        elif self.pdb_code:
            m = re.match(r'.{4}:([A-Z]*)', self.pdb_code)
            self.lines = download_pdb(re.split(":",self.pdb_code)[0]).readlines()

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

        if m:
            self.selection += ' and chain ' + ','.join(m.group(1))

        if self.remove_alternative_locations:
            self.atoms.remove_alternative_locations()
        self.atoms = self.atoms.select(self.selection)

    def __repr__(self):
        return ''.join(self.lines)

    def dssp(self, dssp_command='mkdssp'):
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


def download_pdb(pdb_code, work_dir=expanduser('~'), force_download=False):
    path = work_dir + '/cabsPDBcache/%s' % pdb_code[1:3]
    try: os.makedirs(path)
    except OSError: pass
    fname = path + '/%s.pdb' % pdb_code
    try:
        if force_download:
            raise IOError
        file_ = open(fname)
    except IOError:
        try:
            gz_string = urlopen('http://www.rcsb.org/pdb/files/' + pdb_code.lower() + '.pdb.gz').read()
        except HTTPError:
            raise InvalidPdbCode(pdb_code)
        except URLError:
            raise CannotConnectToPdb()
        with open(fname, 'w') as fobj:
            fobj.write(gz_string)
    file_ = open(fname)
    return GzipFile(fileobj=file_)


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
        return self.pdbCode + ' is not a valid pdb code! (perhaps you specified a file that does not exist)'

class CannotConnectToPdb(Exception):
    """Exception raised when the PDB database is not accessible"""
    def __str__(self):
        return 'Cannot connect to the PDB database!!!'

if __name__ == '__main__':
    pass
