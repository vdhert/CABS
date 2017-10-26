"""Module to handle pdb files."""

import re
import os
import gzip
import json
import requests as req

from tempfile import mkstemp
from time import sleep
from requests.exceptions import HTTPError, ConnectionError
from os.path import expanduser, isfile, join, isdir
from subprocess import Popen, PIPE
from collections import OrderedDict

from CABS import logger
from CABS.atom import Atom, Atoms
from CABS.utils import AA_NAMES, AA_SUB_NAMES

_name = 'PDB'  # module name for logger
_DSSP_COMMAND = 'dssp'


class Pdb(object):
    """
    Pdb parser.
    """
    PDB_CACHE = join(expanduser('~'), 'cabsPDBcache')

    class InvalidPdbInput(Exception):
        pass

    def __init__(
            self,
            source,
            selection='',
            remove_alternative_locations=True,
            fix_non_standard_aa=True,
            remove_water=True,
            remove_hetero=True,
            verify=False,
            no_exit=False,  # does not exit on error, raises InvalidPdbInput instead
    ):

        logger.debug(_name, 'Creating Pdb object from {}'.format(source))
        self.atoms = Atoms()

        words = source.split(':')
        try:
            name, rec, pep = words
            chains = rec + pep
        except ValueError:
            try:
                name, chains = words
            except ValueError:
                name = words[0]
                chains = ''

        try:
            self.body = self.read(name)
            self.name = os.path.basename(name).split('.')[0]
        except IOError:
            try:
                self.body = self.read(self.fetch(name))
                self.name = name
            except ConnectionError as e:
                if no_exit:
                    raise Pdb.InvalidPdbInput(e.message)
                else:
                    logger.exit_program(
                        module_name=_name,
                        msg='Cannot connect to the PDB database',
                        exc=e
                    )
            except HTTPError as e:
                if no_exit:
                    raise Pdb.InvalidPdbInput(e.message)
                else:
                    logger.exit_program(
                        module_name=_name,
                        msg='Invalid PDB code: {}'.format(name),
                        exc=e
                    )
            except IOError as e:
                if no_exit:
                    raise Pdb.InvalidPdbInput(e.message)
                else:
                    logger.exit_program(
                        module_name=_name,
                        msg='File {} not found'.format(name),
                        exc=e
                    )

        try:
            logger.debug(_name, 'Processing {}'.format(name))
            current_model = 0
            for line in self.body.split('\n'):
                match = re.match(r'(ATOM|HETATM)', line)
                if match:
                    self.atoms.append(Atom(line, current_model))
                else:
                    match = re.match(r'MODEL\s+(\d+)', line)
                    if match:
                        current_model = int(match.group(1))

            if chains:
                logger.debug(_name, 'Selected chains {}'.format(chains))
                if selection:
                    selection = '({}) and chain {}'.format(selection, ','.join(chains))
                else:
                    selection = 'chain {}'.format(','.join(chains))

            if remove_alternative_locations:
                logger.debug(_name, 'Removing alternative locations from {}'.format(name))
                self.atoms.remove_alternative_locations()

            if remove_water:
                logger.debug(_name, 'Removing water molecules from {}'.format(name))
                self.atoms = self.atoms.drop('resname HOH')

            if fix_non_standard_aa:
                logger.debug(_name, 'Scanning {} for non-standard amino acids'.format(name))
                aa_names = [AA_NAMES[k] for k in AA_NAMES]
                for residue in self.atoms.residues():
                    resname = residue[0].resname
                    if resname not in aa_names:
                        if resname not in AA_SUB_NAMES:
                            logger.warning(
                                module_name=_name,
                                msg='Unknown residue {} at {} in {}'.format(
                                    resname, residue[0].resid_id(), name
                                )
                            )
                        else:
                            sub_name = AA_SUB_NAMES[resname]
                            for atom in residue:
                                atom.resname = sub_name
                                atom.hetatm = False
                            logger.warning(
                                module_name=_name,
                                msg='Replacing {} -> {} for {} in {}'.format(
                                    resname, sub_name, residue[0].resid_id(), name
                                )
                            )

            if remove_hetero:
                logger.debug(_name, 'Removing heteroatoms from {}'.format(name))
                self.atoms = self.atoms.drop('hetero')

            if selection:
                logger.debug(_name, 'Selecting [{}] from {}'.format(selection, name))
                self.atoms = self.atoms.select(selection)

            if not len(self.atoms):
                raise Exception('{} contains no atoms'.format(source))

            if chains and verify:
                actual_chains = ''.join(self.atoms.list_chains().keys())
                logger.debug(
                    module_name=_name,
                    msg='Matching declared [{}] with actual [{}] chain IDs in {}.'.format(chains, actual_chains, name)
                )
                if set(chains) != set(actual_chains):
                    msg = 'Mismatch in chain IDs in {}: {} differs from {}'.format(name, chains, actual_chains)
                    logger.warning(_name, msg)
                    raise Exception(msg)

        except Exception as e:
            if no_exit:
                raise Pdb.InvalidPdbInput(e.message)
            else:
                logger.exit_program(
                    module_name=_name,
                    msg=e.message,
                    exc=e
                )

    @staticmethod
    def fetch(pdb_code, force_download=False):

        if not re.match(r'[1-9][0-9A-Za-z]{3}', pdb_code):
            raise IOError

        pdb_low = pdb_code.lower()
        path = join(Pdb.PDB_CACHE, pdb_low[1:3])
        try:
            os.makedirs(path)
        except OSError:
            pass

        filename = join(path, '%s.pdb.gz' % pdb_low)

        if not isfile(filename) or force_download:
            logger.debug(_name, 'Downloading {}'.format(pdb_low))
            url = 'http://files.rcsb.org/download/%s.pdb.gz' % pdb_low
            r = req.get(url)
            r.raise_for_status()
            with open(filename, 'wb') as f:
                f.write(r.content)

        return filename

    @staticmethod
    def read(filename):
        try:
            with gzip.open(filename, 'rb') as f:
                content = f.read()
        except IOError:
            with open(filename, 'rb') as f:
                content = f.read()
        return content

    def dssp(self, output=''):
        """Runs dssp on the read pdb file and returns a dictionary with secondary structure"""

        out = err = None

        try:
            proc = Popen([_DSSP_COMMAND, '/dev/stdin'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
            out, err = proc.communicate(input=self.body)
            logger.debug(
                module_name=_name,
                msg='Running DSSP'
            )
        except OSError:
            logger.warning(
                module_name=_name,
                msg='DSSP not found.'
            )

            tempfile = mkstemp(suffix='.pdb', prefix='.tmp.dssp.', dir=Pdb.PDB_CACHE)[1]
            with open(tempfile, 'wb') as f:
                f.write(self.body)

            try:
                logger.debug(
                    module_name=_name,
                    msg='Submitting structure to the DSSP server'
                )
                out, err = self.xssp(tempfile)

            except (HTTPError, ConnectionError):
                logger.warning(
                    module_name=_name,
                    msg='Cannot connect to the DSSP server'
                )
            finally:
                try:
                    os.remove(tempfile)
                except OSError:
                    pass

        if err:
            logger.critical(
                module_name=_name,
                msg='DSSP ERROR: %s' % err.replace('\n', ' ')
            )
            return None
        else:
            logger.debug(_name, 'DSSP successful')
            if logger.log_level >= 2 and output:
                output = os.path.join(output, 'output_data', 'DSSP_output_%s.txt' % self.name)
                d = os.path.dirname(output)
                if not isdir(d):
                    os.makedirs(d)
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
                key = str(m.group(2).strip() + ':' + m.group(3))
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

    @staticmethod
    def xssp(filename, server='http://www.cmbi.umcn.nl/xssp'):
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

    def __str__(self):
        return self.body

    def __repr__(self):
        return "<PDB from %s, %i atoms>" % (self.name, len(self.atoms))
