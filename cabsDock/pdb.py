"""Module to handle pdb files."""

import re
import os
import logger

from os.path import exists, expanduser
from os.path import expanduser
from gzip import GzipFile
from urllib2 import urlopen, HTTPError, URLError
from subprocess import Popen, PIPE

from atom import Atom, Atoms

__all__ = ["PDB"]
_name = "PDB"

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
        self.name = ""


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
            logger.exit_program(module_name=_name, msg="No PDB file/code provided. Quitting.",traceback=False)

        if self.file_name:
            m = re.match(r'[^:]*:([A-Z]*)', self.file_name)
            self.file_name = re.split(":",self.file_name)[0]
            try:
                self.lines = GzipFile(self.file_name).readlines()
            except IOError:
                self.lines = open(self.file_name).readlines()
            self.name = self.file_name.split(".")[0]
        elif self.pdb_code:
            m = re.match(r'.{4}:([A-Z]*)', self.pdb_code)
            self.name = re.split(":", self.pdb_code)[0]
            self.lines = download_pdb(self.name).readlines()

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
                logger.exit_program(module_name=_name, msg="Failed to read %s (perhaps its not a valid pdp file?). Quitting." % self.file_name)
            elif self.pdb_code:
                logger.exit_program(module_name=_name, msg="Failed to read %s. Quitting." % self.pdb_code)

        if m:
            self.selection += ' and chain ' + ','.join(m.group(1))

        if self.remove_alternative_locations:
            self.atoms.remove_alternative_locations()
        self.atoms = self.atoms.select(self.selection)

    def __repr__(self):
        return ''.join(self.lines)

    def dssp(self, dssp_command='mkdssp', output = ''):
        """Runs dssp on the read pdb file and returns a dictionary with secondary structure"""
        try:
            proc = Popen([dssp_command, '/dev/stdin'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
        except OSError:
            logger.exit_program(module_name=_name,
                                msg="DSSP was not found. Quitting.",
                                traceback=False)

        logger.info(module_name=_name,
                    msg = "DSSP running on %s. Selected chains: %s" % (self.name, "".join(self.atoms.list_chains().keys())))
        out, err = proc.communicate(input=''.join(self.lines))
        if err:
            logger.critical(module_name=_name, msg="DSSP returned an error: %s" % err)
            return None
        else:
            if logger.log_level >=2 and output:
                output += "/output_data/DSSP_output_%s.txt" % self.name
                logger.to_file(filename=output, content=out, msg="Saving DSSP output to %s" % output)


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
        except URLError as e:
            logger.exit_program(module_name=_name,
                                msg="Could not download the pdb file. Can't connect to the PDB database - quitting",
                                traceback=True,exc=e)
        with open(fname, 'w') as fobj:
            fobj.write(gz_string)
    file_ = open(fname)
    return GzipFile(fileobj=file_)


class InvalidPdbCode(Exception):
    """Exception raised when invalid pdb_code is used with download_pdb"""
    def __init__(self, pdb_code):
        self.pdbCode = pdb_code

    def __str__(self):
        return self.pdbCode + ' is not a valid pdb code! (perhaps you specified a file that does not exist)'




if __name__ == '__main__':
    pass
