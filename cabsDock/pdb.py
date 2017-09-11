"""Module to handle pdb files."""

import re
import os
import logger
import gzip
import json
import requests as req
from tempfile import mkstemp
from time import sleep
from requests.exceptions import HTTPError, ConnectionError
from os.path import expanduser, basename, isfile, join, isdir
from subprocess import Popen, PIPE
from atom import Atom, Atoms
from collections import OrderedDict

_name = 'PDB'  # module name for logger

PDB_CACHE = join(expanduser('~'), 'cabsPDBcache')
try:
    os.makedirs(PDB_CACHE)
except OSError:
    pass


class Pdb:
    """
    Pdb parser.

    Usage:
    pdb_file=<PATH TO FILE>:<chains>
    pdb_code=<PDB_CODE>:<chains>
    """

    def __init__(self, remove_alternative_locations=True, **kwargs):
        self.file_name = kwargs.get('pdb_file')
        self.pdb_code = kwargs.get('pdb_code')
        self.selection = kwargs.get('selection')
        self.remove_alternative_locations = remove_alternative_locations
        self.atoms = Atoms()
        self.header = []
        self.missed = OrderedDict()
        chains = None

        if self.file_name:
            if ':' in self.file_name:
                filename, chains = self.file_name.split(':')
            else:
                filename = self.file_name
            try:
                logger.info(
                    module_name=_name,
                    msg='Loading file %s' % filename
                )
                self.lines = gzip.open(filename).readlines()
            except IOError:
                try:
                    self.lines = open(filename).readlines()
                except IOError:
                    logger.exit_program(
                        module_name=_name,
                        msg='ERROR: Cannot read %s.' % filename,
                        traceback=False
                    )
            self.name = basename(filename).rsplit('.', 1)[0]

        elif self.pdb_code:
            if ':' in self.pdb_code:
                pdbcode, chains = self.pdb_code.split(':')
            else:
                pdbcode = self.pdb_code
            try:
                self.lines = download_pdb(pdbcode).readlines()
            except IOError:
                logger.exit_program(
                    module_name=_name,
                    msg='ERROR: Could not download %s.' % pdbcode,
                    traceback=False
                )
            self.name = pdbcode[:4]
        else:
            logger.exit_program(
                module_name=_name,
                msg='ERROR: No PDB file/code provided. Quitting.',
                traceback=False
            )

        current_model = 0
        for line in self.lines:
            match = re.match(r'^(ATOM|HETATM)', line)
            if match:
                self.atoms.append(Atom(line, current_model))
            elif line[:5] == 'MODEL':
                current_model = int(line.split()[1])
            elif line[:3] == 'END' or line[:3] == 'TER':
                pass
            elif len(self.atoms) == 0:
                self.header.append(line)
            else:
                if current_model not in self.missed:
                    self.missed[current_model] = []
                self.missed[current_model].append(line)

        if len(self.atoms):
            logger.debug(
                module_name=_name,
                msg='Loaded %i atoms from %s' % (len(self.atoms), self.name)
            )
        else:
            logger.exit_program(
                module_name=_name,
                msg='ERROR: Structure %s contains no atoms.' % self.name,
                traceback=False
            )

        if chains:
            if not self.selection:
                self.selection = 'chain %s' % chains
            else:
                self.selection = '(%s) and chain %s' % (self.selection, chains)

        if self.remove_alternative_locations:
            self.atoms.remove_alternative_locations()
            logger.debug(
                module_name=_name,
                msg='Removing atoms at alternative locations'
            )
        if self.selection:
            self.atoms = self.atoms.select(self.selection)
            logger.debug(
                module_name=_name,
                msg='Selection: \'%s\'' % self.selection
            )

    def __repr__(self):
        return ''.join(self.lines)

    def dssp(self, dssp_command='mkdssp', output=''):
        """Runs dssp on the read pdb file and returns a dictionary with secondary structure"""

        out = err = ''

        try:
            proc = Popen([dssp_command, '/dev/stdin'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate(input=str(self))
            logger.info(
                module_name=_name,
                msg='DSSP running on %s.' % self.name
            )
        except OSError:
            logger.warning(
                module_name=_name,
                msg='DSSP not found.'
            )

            tempfile = mkstemp(suffix='.pdb', prefix='.tmp.dssp.', dir=PDB_CACHE)[1]
            with open(tempfile, 'wb') as f:
                f.write(str(self))

            try:
                logger.debug(
                    module_name=_name,
                    msg='Submitting structure to the DSSP server'
                )
                out, err = dssp_server(tempfile)

            except (HTTPError, ConnectionError):
                logger.warning(
                    module_name=_name,
                    msg='Cannot connect to the DSSP server. II structure set to \'coil\' for %s' % self.name,
                )
            finally:
                try:
                    os.remove(tempfile)
                except OSError:
                    pass

        if err:
            logger.critical(
                module_name=_name,
                msg='DSSP returned an error: %s' % err
            )
            return None
        else:
            if logger.log_level >= 2 and output:
                output += '/output_data/DSSP_output_%s.txt' % self.name
                logger.to_file(
                    filename=output,
                    content=out,
                    msg='Saving DSSP output to %s' % output
                )

        sec = OrderedDict()
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


def download_pdb(pdb_code, force_download=False):
    path = join(PDB_CACHE, pdb_code[1:3])
    try:
        os.makedirs(path)
    except OSError:
        pass
    filename = join(path, '%s.pdb.gz' % pdb_code)

    if not isfile(filename) or force_download:
        logger.info(_name, 'Downloading %s' % pdb_code)
        url = 'http://files.rcsb.org/download/%s.pdb.gz' % pdb_code
        try:
            r = req.get(url)
            r.raise_for_status()
            with open(filename, 'wb') as f:
                f.write(r.content)
        except HTTPError:
            logger.exit_program(
                module_name=_name,
                msg='ERROR: Invalid PDB code %s.' % pdb_code,
                traceback=False
            )
        except ConnectionError:
            logger.exit_program(
                module_name=_name,
                msg='ERROR: Cannot connect to the PDB database.' % pdb_code,
                traceback=False
            )
    return gzip.open(filename, 'rb')


def dssp_server(filename, server='http://www.cmbi.umcn.nl/xssp'):
    url_api = server + '/api/%s/pdb_file/dssp/'

    files = {'file_': open(filename, 'rb')}

    r = req.post(url=url_api % 'create', files=files)
    r.raise_for_status()
    job_id = json.loads(r.content)['id']

    while True:
        r = req.get(url_api % 'status' + job_id)
        r.raise_for_status()
        status = json.loads(r.content)['status']

        if status == 'SUCCESS':
            r = req.get(url_api % 'result' + job_id)
            r.raise_for_status()
            out = json.loads(r.content)['result']
            err = ''
            break
        elif status in ['FAILURE', 'REVOKED']:
            err = json.loads(r.content)['message']
            out = ''
            break
        else:
            sleep(1)

    return out, err


if __name__ == '__main__':
    pass