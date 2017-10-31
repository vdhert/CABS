import argparse
import textwrap
import re
from copy import deepcopy as dc

_HELPW = 100
_wrapper = textwrap.TextWrapper(width=_HELPW, break_long_words=True, expand_tabs=False)


def _wrap(text):
    lines = []
    for line in text.split('\n'):
        _list = _wrapper.wrap(line)
        if not _list:
            _list = ['\n']
        lines.extend(_list)
    return '\n'.join(lines)


class CABSFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, prog, indent_increment=2, max_help_position=4, width=_HELPW):
        super(CABSFormatter, self).__init__(prog, indent_increment, max_help_position, width)

    def _split_lines(self, text, width):
        return super(CABSFormatter, self)._split_lines(text, width) + [' ']


class ConfigFileParser:
    OPTIONRE = re.compile(
        r'(?P<option>[^:=]*)'
        r'[:=]'
        r'(?P<value>.*)$'
    )

    def __init__(self, filename):
        self.args = []
        with open(filename) as f:
            for line in f:
                if line == '' or line[0] in ';#\n':
                    continue
                match = self.OPTIONRE.match(line)
                if match:
                    option, value = match.groups()
                    self.args.append('--' + option.strip())
                    self.args.extend(value.split('#')[0].split(';')[0].split())
                else:
                    self.args.extend(line.split('#')[0].split(';')[0].split())


def mk_parser(parser_dict, group_dict, option_dict):

    parser_dict['description'] = _wrap(parser_dict['description'])
    parser_dict['epilog'] = _wrap(parser_dict['epilog'])

    groups = parser_dict.pop('groups')
    defaults = parser_dict.pop('defaults', {})
    parser = argparse.ArgumentParser(
        formatter_class=CABSFormatter,
        add_help=False,
        **parser_dict
    )
    for group_name, options in groups:
        group = group_dict[group_name]
        group['description'] = _wrap(group['description'])
        group = parser.add_argument_group(title=group_name, **group)
        for opt_name in options:
            option = next(v for k, v in option_dict.items() if k == opt_name)
            name = ['--' + opt_name]
            flag = option.pop('flag', None)
            if flag:
                name.insert(0, flag)
            if opt_name in defaults:
                option['default'] = defaults[opt_name]
            option['help'] = _wrap(option['help'])
            group.add_argument(*name, **option)
    return parser


dock_dict = {
    'prog': 'CABSdock',
    'description':
        'CABSdock application is a versatile tool for molecular docking of peptides to proteins. It allows for flexible'
        ' docking (also with large-scale conformational changes) and without the knowledge about the binding site. '
        'CABSdock enables peptide docking using only information about the peptide sequence and the protein protein '
        'structure. Additionally many advanced options are available that allow for manipulation of the simulation '
        'setup, the degree of protein flexibility or guiding the peptide binding etc.\n\n'
        'CABSdock method has been first made available as a web server [Nucleic Acids Research, 43(W1): W419-W424, 2015'
        '; web server website: http://biocomp.chem.uw.edu.pl/CABSdock]. The standalone application [submitted to '
        'publication] provides the same modeling methodology equipped with many additional features and customizable '
        'options.',
    'epilog': 'CABSdock repository: https://bitbucket.org/lcbio/cabsdock',
    'defaults': {'temperature': (2.0, 1.0), 'replicas': 10, 'protein-restraints': ('all', 5, 5.0, 15.0)},
    'groups': [
        ('BASIC OPTIONS', ['input-protein', 'peptide', 'config']),
        ('PROTEIN OPTIONS',['exclude', 'excluding-distance', 'protein-flexibility',
                            'protein-restraints', 'weighted-fit', 'weighted-fit-file']),
        ('PEPTIDE OPTIONS', ['add-peptide', 'separation', 'insertion-clash', 'insertion-attempts']),
        ('RESTRAINTS OPTIONS', ['ca-rest-add', 'sc-rest-add', 'ca-rest-weight',
                                'sc-rest-weight', 'ca-rest-file','sc-rest-file']),
        ('SIMULATION OPTIONS', ['mc-annealing', 'mc-cycles', 'mc-steps', 'replicas',
                                'replicas-dtemp', 'temperature', 'random-seed']),
        ('ALL_ATOM RECONSTRUCTION OPTIONS', ['modeller-iterations', 'aa-rebuild']),
        ('ANALYSIS OPTIONS', ['reference-pdb', 'filtering-count', 'filtering-mode', 'clustering-medoids',
                              'clustering-iterations', 'contact-maps', 'contact-threshold', 'align',
                              'align-options', 'align-peptide-options']),
        ('OUTPUT OPTIONS', ['save-cabs-files', 'load-cabs-files', 'save-config', 'pdb-output']),
        ('MISCELLANEOUS OPTIONS', ['work-dir',  'dssp-command', 'fortran-command', 'image-file-format',
                                   'contact-map-colors', 'verbose', 'version', 'help'])
    ]
}


flex_dict = {
    'prog': 'CABSflex',
    'description':
        'CABSflex: versatile tool for the simulation of structure flexibility of folded globular proteins.',
    'epilog': 'CABSdock repository: https://bitbucket.org/lcbio/cabsdock',
    'defaults': {'temperature': (1.4, 1.4), 'replicas': 1, 'protein-restraints': ('ss2', 3, 3.8, 8.0)},
    'groups': [
        ('BASIC OPTIONS', ['input-protein', 'config']),
        ('PROTEIN OPTIONS',['exclude', 'excluding-distance', 'protein-flexibility',
                            'protein-restraints', 'weighted-fit', 'weighted-fit-file']),
        ('RESTRAINTS OPTIONS', ['ca-rest-add', 'sc-rest-add', 'ca-rest-weight',
                                'sc-rest-weight', 'ca-rest-file','sc-rest-file']),
        ('SIMULATION OPTIONS', ['mc-annealing', 'mc-cycles', 'mc-steps', 'replicas',
                                'replicas-dtemp', 'temperature', 'random-seed']),
        ('ALL_ATOM RECONSTRUCTION OPTIONS', ['modeller-iterations', 'aa-rebuild']),
        ('ANALYSIS OPTIONS', ['reference-pdb', 'filtering-count', 'filtering-mode', 'clustering-medoids',
                              'clustering-iterations', 'contact-maps', 'contact-threshold', 'align',
                              'align-options', 'align-peptide-options']),
        ('OUTPUT OPTIONS', ['save-cabs-files', 'load-cabs-files', 'save-config', 'pdb-output']),
        ('MISCELLANEOUS OPTIONS', ['work-dir',  'dssp-command', 'fortran-command', 'image-file-format',
                                   'contact-map-colors', 'verbose', 'version', 'help'])
    ]
}


groups = {
    'BASIC OPTIONS': {
        'description':
            'Basic options allow tu run protein-peptide docking using only information about the peptide sequence and '
            'the protein structure. In addition to peptide sequence, a peptide secondary structure information may be '
            'provided.'
    },
    'PROTEIN OPTIONS': {
        'description':
            'Protein options allow to exclude some of the protein residues from docking search and to manipulate the '
            'degree of protein flexibility (using either --protein-flexibility or --protein-restraints).'
    },
    'PEPTIDE OPTIONS': {
        'description':
            'Peptide options, apart from providing peptide sequence, allow to introduce starting peptide '
            'conformations and locations (when these additional options are not used, starting peptide '
            'conformations and locations are random).'
    },
    'RESTRAINTS OPTIONS': {
        'description':
            'Restraints options allow to set up distance restraints, either between C-alpha atoms (CA) or Side Chains '
            '(SC), where SCs are geometric centers of side chains atoms (as defined in the CABS coarse-grained model).'
    },
    'SIMULATION OPTIONS': {
        'description':
            'Simulation options allow to modify different parameters of Replica Exchange Monte Carlo simulation '
            'procedure.'
    },
    'ALL_ATOM RECONSTRUCTION OPTIONS': {
        'description':
            'All-atom reconstruction options allow to set up details of all-atom reconstruction and refinement '
            'procedure, which is performed using MODELLER package.'
        },
    'ANALYSIS OPTIONS': {
        'description':
            'Analysis options allow to perform comparison analyses (to provided reference complex structure) and for '
            'repeated scoring and analysis of CABSdock trajectories.'
        },
    'OUTPUT OPTIONS': {'description': 'Output options.'},
    'MISCELLANEOUS OPTIONS': {'description': 'Miscellaneous options.'}
}


options = {
    'aa-rebuild': {
        'flag': '-A',
        'action': 'store_true',
        'help':
            'Rebuild final models to all-atom representation. (default: %(default)s)'
    },
    'add-peptide': {
        'flag': '-P',
        'nargs': 3,
        'action': 'append',
        'metavar': ('PEPTIDE', 'CONFORMATION', 'LOCATION'),
        'help':
            'Add peptide to the complex. This option can be used multiple times to add multiple peptides.\n\n'
            'PEPTIDE must be either:\n\n'
            '[1] amino acid sequence in one-letter code (optionally annotated with secondary structure: '
            'H - helix, E - sheet, C - coil)\n'
            'i.e. \'-p HKILHRLLQD:CHHHHHHHHC\' loads HKILHRLLQD peptide sequence with the secondary structure '
            'assignment: CHHHHHHHHC\n\n'
            'HINT: If possible, it is always recommended to use secondary structure information/prediction. For '
            'residues with ambiguous secondary structure prediction assignment it is better to assign coil (C) than '
            'the regular (H - helix or E - extended) type of structure.\n\n'
            '[2] pdb file (may be gzipped)\n'
            '[3] pdb code (optionally with chain_id i.e. 1abc:D)\n\n'
            'CONFORMATION sets initial conformation of the peptide. Must be either:\n\n'
            '[1] \'random\' - random conformation is generated\n'
            '[2] \'keep\' - preserve conformation from file. This has no effect if PEPTIDE=SEQUENCE.\n\n'
            'LOCATION sets initial location for the peptide. Must be either:\n\n'
            '[1] \'random\' - peptide is placed in a random location on the surface of a sphere centered at the '
            'protein\'s geometrical center at distance defined by the \'--separation\' option from the surface of '
            'the protein.\n'
            '[2] keep - preserve location from file. This has no effect if PEPTIDE=SEQUENCE\n'
            '[3] PATCH - list of protein\'s residues (i.e 123:A+125:A+17:B). Peptide will be placed above the '
            'geometrical center of listed residues at distance defined by the \'--separation\' option from the '
            'surface of the protein.\n'
            'WARNING: residues listed in patch should be on the surface of the protein and close to each other.'
    },
    'align': {
        'default': 'SW',
        'metavar': 'METHOD',
        'help':
            'Select the method to be used to align target with reference pdb.\n\n'
            'Available options are: (default %(default)s)\n\n'
            '[1] SW - Smith-Waterman\n'
            '[2] blastp - protein BLAST (requires NCBI+ package installed)\n'
            '[3] trivial - simple sequential alignment, used only in case of one-chain input and reference of the same '
            'length'
    },
    'align-options': {
        'metavar': 'KEY=VAL',
        'nargs': '+',
        'default': [],
        'type': lambda x: x.split('='),
        'help': 'Path to alignment with reference structure. If set, the \'--align\' option is ignored'
    },
    'align-peptide-options': {
        'metavar': 'KEY=VAL',
        'nargs': '+',
        'default': [],
        'type': lambda x: x.split('='),
        'help': 'Path to alignment with reference structure. If set, the \'--align\' option is ignored'
    },
    'ca-rest-add': {
        'action': 'append',
        'metavar': ('RESI', 'RESJ', 'DIST', 'WEIGHT'),
        'nargs': 4,
        'help':
            'Add distance restraint between C-alpha (CA) atom in residue RESI and C-alpha atom in residue RESJ.\n'
            'DIST is a distance between these atoms and WEIGHT is restraint\'s weight (number from [0, 1]).\n'
            'In order to add restraints between the peptide and the protein, or between two peptides, use PEP1, '
            'PEP2, ... as chain identifiers of the peptides (even when peptide is read from a pdb file its chain '
            'identifier is ignored).\n'
            'i.e. \'--ca-rest-add 123:A 5:PEP1 8.7 1.0\' adds restraint between the C-alpha atom of the residue number '
            '123 in the chain A of the protein and the C-alpha atom of the 5th residue of the peptide. If you add only '
            'one peptide both \'PEP\' and \'PEP1\' is a valid chain identifier. If you add multiple peptides they will '
            'be ordered as follows:\n\n'
            '[1] from config file added by the \'peptide\' option\n'
            '[2] from config file added by the \'add-peptide\' option\n'
            '[3] from command line added by the \'-p, --peptide\' option\n'
            '[4] from command line added by the \'--add-peptide\' option.\n\n'
            'Peptides added by the same method preserve the order by which they appear in the config file, or on the '
            'command line. Option can be used multiple times to add multiple restraints.'
    },
    'ca-rest-file': {
        'action': 'append',
        'metavar': 'FILE',
        'help': 'Read C-alpha restraints from file (use multiple times to add multiple files)'
    },
    'ca-rest-weight': {
        'default': 1.0,
        'type': float,
        'metavar': 'WEIGHT',
        'help':
            'Set global weight for all C-alpha restraints (including automatically generated restraints for the '
            'protein) (default: %(default)s)',
    },
    'clustering-iterations': {
        'default': 100,
        'type': int,
        'metavar': 'NUM',
        'help':
            'Set number of iterations of the k-medoids clustering algorithm. (default: %(default)s)'
    },
    'clustering-medoids': {
        'flag': '-k',
        'default': 10,
        'type': int,
        'metavar': 'NUM',
        'help':
            'Set number of medoids in k-medoids clustering algorithm. This option also sets number of final models '
            'to be generated. (default: %(default)s)'
    },
    'config': {
        'flag': '-c',
        'metavar': 'FILE',
        'help': 'read options from FILE'
    },
    'contact-maps': {
        'flag': '-M',
        'action': 'store_true',
        'help': 'Store contact maps matrix plots and histograms of contact frequencies. (default: %(default)s)'
    },
    'contact-map-colors': {
        'type': str,
        'nargs': 6,
        'default': ['#ffffff', '#f2d600', '#4b8f24', '#666666', '#e80915', '#000000'],
        'metavar': 'HEX',
        'help':
            'Sets 6 colors (hex code, e.g. #00FF00 for green etc.) to be used in contact map color bars.\n'
            '(default: %(default)s)'
    },
    'contact-threshold': {
        'flag': '-T',
        'default': 6.5,
        'type': float,
        'metavar': 'DIST',
        'help':
            'Set contact distance between side chains pseudoatoms (SC) for contact map plotting. (default: %(default)s)'
    },
    'dssp-command': {
        'default': 'dssp',
        'metavar': 'PATH',
        'help': 'Provide path to dssp binary (default is \'%(default)s\')'
    },
    'exclude': {
        'flag': '-e',
        'action': 'append',
        'metavar': 'RESIDUES',
        'help':
            'Exclude proteins residues listed in RESIDUES from the docking search, therefore enforces more effective '
            'search in other areas of the protein surface, for example, it may be known that some parts of the protein '
            'are not accessible to peptide (due to binding to other proteins) and therefore it could be useful to '
            'exclude these regions from the search procedure.\n\n'
            'RESIDUES must be a single string of characters (no whitespaces) consisting of residue identifiers (i. e. '
            '123:A) or chain identifiers (i. e. A) joined with the \'+\' sign. \'-\' is also allowed to specify a '
            'continuous range of residues, or chains.\n\n'
            'Examples:\n'
            '\'-e 123:A\'            excludes residue 123 from chain A\n'
            '\'-e 123:A+125:A\'      residues 123 and 125 from chain A\n'
            '\'-e 123:A-125:A\'      residues 123, 124 and 125 from chain A\n'
            '\'-e A\'                whole chain A\n'
            '\'-e A+C\'              chains A and C\n'
            '\'-e A-C\'              chains A, B and C\n\n'
            'Adding @PEP<N> at the end of the string limits the excluding to only N-th peptide.\n'
            'i.e. \'-e 123:A@PEP1\' will exclude residue 123 in chain A for binding with the first peptide only.\n'
            'If @PEP<N> is omitted the exclusion list affects all peptides.\n\n'
            'This option can be used multiple times to add multiple sets of excluded residues.'
    },
    'excluding-distance': {
        'default': 5.0,
        'type': float,
        'metavar': 'DISTANCE',
        'help':
            'Set minimum distance between side chain atoms of peptide(s) and proteins residues marked as \'excluded\''
            '(default: %(default)s)'
    },
    'filtering-count': {
        'flag': '-n',
        'default': 1000,
        'type': int,
        'metavar': 'NUM',
        'help': 'Set number of low-energy models from trajectories to be clustered (default %(default)s)'
    },
    'filtering-mode': {
        'choices': ['each', 'all'],
        'default': 'each',
        'metavar': 'MODE',
        'help':
            'Choose the filtering mode to select NUM (set by \'--filtering-count\') models for clustering.\n\n'
            'MODE can be either: (default: %(default)s)\n\n'
            '[1] \'each\' - models are ordered by protein-peptide(s) binding energy and top n = [NUM / R] (R is the '
            'number of replicas) is selected from EACH replica\n'
            '[2] \'all\' - models are ordered by protein-peptide(s) binding energy and top NUM is selected from ALL '
            'replicas combined'
    },
    'fortran-command': {
        'default': 'gfortran',
        'metavar': 'PATH',
        'help': 'Provide path to fortran compiler binary (default is \'%(default)s\')'
    },
    'help': {
        'flag': '-h',
        'action': 'store_true',
        'help': 'Print help and exit program'
    },
    'image-file-format': {
        'default': 'svg',
        'metavar': 'FMT',
        'help': 'Save all of the image files in given format. (default: %(default)s)'
    },
    'input-protein': {
        'flag': '-i',
        'metavar': 'INPUT',
        'required': True,
        'help':
            'Load input protein structure.\n\n'
            'INPUT can be either:\n\n'
            '[1] PDB code (optionally with chain IDs)\n'
            'i.e. \'-P 1CE1:HL\' loads chains H and L of 1CE1 protein structure downloaded from the PDB database\n'
            '[2] PDB file (optionally gzipped)',
    },
    'insertion-attempts': {
        'default': 1000,
        'type': int,
        'metavar': 'NUM',
        'help':
            'This option enables advanced settings of building starting conformations of modelled complexes. The '
            'option sets number of attempts to insert peptide while building inital complex (default: %(default)s)',
    },
    'insertion-clash': {
        'default': 1.0,
        'type': float,
        'metavar': 'DIST',
        'help':
            'This option enables advanced settings of building starting conformations of modelled complexes. The '
            'option sets distance in Angstroms between any two atoms (of different modeled chains) at which a clash '
            'occurs while building initial complex (default: %(default)s Angstrom)'
    },
    'load-cabs-files': {
        'flag': '-L',
        'metavar': 'FILE',
        'nargs': 2,
        'help':
            'Load CABSdock simulation data file for repeated scoring and analysis of CABSdock trajectories (with new '
            'settings, for example using a reference complex structure and --reference option).'
    },
    'mc-annealing': {
        'flag': '-a',
        'default': 20,
        'type': int,
        'metavar': 'NUM',
        'help': 'Set number of Monte Carlo temperature annealing cycles to NUM (NUM > 0, default: %(default)s)'
    },
    'mc-cycles': {
        'flag': '-y',
        'default': 50,
        'type': int,
        'metavar': 'NUM',
        'help':
            'Set number of Monte Carlo cycles to NUM (NUM>0, default: %(default)s).\n'
            'Total number of snapshots generated for each replica/trajectory = [mc-annealing] x [mc-cycles], '
            '(default: 20 x 50 = 1000)'
    },
    'mc-steps': {
        'flag': '-s',
        'default': 50,
        'type': int,
        'metavar': 'NUM',
        'help':
            'Set number of Monte Carlo cycles between trajectory frames to NUM (NUM > 0, default: %(default)s).\n'
            'NUM = 1 means that every generated conformation will occur in trajectory. This option enables to increase '
            'the simulation length (between printed snapshots) and doesn\'t impact the number of snapshots in '
            'trajectories (see also --mc-cycles description).'
    },
    'modeller-iterations': {
        'flag': '-m',
        'default': 3,
        'type': int,
        'metavar': 'NUM',
        'help':
            'Set number of iterations for reconstruction procedure in MODELLER package (default: %(default)s). Bigger '
            'numbers may result in more accurate models, but reconstruction takes longer.'
    },
    'pdb-output': {
        'flag': '-o',
        'default': 'A',
        'metavar': 'SELECTION',
        'help':
            'Select structures to be saved in the pdb format.\n\n'
            'Available options are:\n'
            '[1] \'A\' - all (default)\n'
            '[2] \'R\' - replicas\n'
            '[3] \'F\' - filtered\n'
            '[4] \'C\' - clusters\n'
            '[5] \'M\' - models\n'
            '[6] \'N\' - none\n\n'
            'i. e. \'-o RM\' - saves replicas and models'
    },
    'peptide': {
        'flag': '-p',
        'action': 'append',
        'metavar': 'PEPTIDE',
        'help':
            'Load peptide sequence and optionally peptide secondary structure in one-letter code (can be used multiple'
            ' times to add multiple peptides).\n\n'
            'PEPTIDE can be either:\n\n'
            '[1] amino acid sequence in one-letter code (optionally annotated with secondary structure: '
            'H - helix, E - sheet, C - coil)\n'
            'i.e. \'-p HKILHRLLQD:CHHHHHHHHC\' loads HKILHRLLQD peptide sequence with the secondary structure '
            'assignment: CHHHHHHHHC\n\n'
            'HINT: If possible, it is always recommended to use secondary structure information/prediction. For '
            'residues with ambiguous secondary structure prediction assignment it is better to assign coil (C) than the'
            'regular (H - helix or E - extended) type of structure.\n\n'
            '[2] PDB code (optionally with chain ID)\n'
            'i.e. \'-p 1CE1:P\' loads the sequence of the chain P from 1CE1 protein\n'
            '[3] PDB file with peptide\'s coordinates, loads only a peptide sequence from a PDB file\n\n'
            '\'--peptide PEPTIDE\' is an alias for \'--add-peptide PEPTIDE random random\''
    },
    'protein-flexibility': {
        'flag': '-f',
        'default': 1.0,
        'metavar': 'FLEXIBILITY',
        'help':
            'Modify flexibility of selected protein\'s residues:\n\n'
            'f = 0 - fully flexible backbone\n'
            'f = 1 - almost stiff backbone (default value)\n\n'
            'FLEXIBILITY can be either:\n\n'
            '[1] positive real number - all protein residues are assigned the same flexiblity equal that number\n'
            '[2] \'bf\' - flexibility for each residue is read from the beta factor column of the CA atom in the '
            'pdb input file\n'
            '(Note that standard beta factor in pdb file has opposite meaning to CABSdock flexibillity, '
            'edit pdb accordingly or use FLEXIBILITY = \'bfi\')\n'
            '[3] \'bfi\' - flexibility is assigned from the inverted beta factors in the input pdb file so that:\n'
            'bf <= 0.0 -> f = 1.0; bf >= 1.0 -> f = 0.0; f = 1 - bf otherwise.\n'
            '[4] \'bfg\' - flexibility is assigned from the beta factors in the input pdb file, so that:\n'
            'f = exp(-bf * bf / 2.) and f = 1.0 if bf < 0.0\n'
            '[5] <filename> - flexibility is read from file <filename> in the following format:\n\n'
            'default <default flexibility value> (if omitted default f = %(default)s)\n'
            'resid_ID <flexibility> i.e. 12:A 0.75 OR resid_ID - resid_ID <flexibility> i.e. 12:A - 15:A  0.75\n\n'
            'Multiple entries can be used.',
    },
    'protein-restraints': {
        'flag': '-g',
        'default': ('all', 5, 5.0, 15.0),
        'nargs': 4,
        'metavar': ('MODE', 'GAP', 'MIN', 'MAX'),
        'help':
            'This options allows to generate a set of binary distance restraints for C-alpha atoms, that keep the '
            'protein in predefined conformation. (default: %(default)s)\n\n'
            'MODE can be either:\n'
            '[1] \'all\' - generate restraints for all protein residues\n'
            '[2] \'ss1\' - generate restraints only when at least one restrained residue is assigned regular secondary '
            'structure (helix or sheet)\n'
            '[3] \'ss2\' - generate restraints only when both restrained residues are assigned regular secondary '
            'structure (helix, sheet)\n\n'
            'GAP specifies minimal gap along the main chain for two resiudes to be restrained.\n'
            'MIN and MAX are min and max values in Angstroms for two residues to be restrained.'
    },
    'random-seed': {
        'flag': '-z',
        'type': int,
        'metavar': 'SEED',
        'help': 'Set seed for random number generator.'
    },
    'reference-pdb': {
        'flag': '-R',
        'metavar': 'REF',
        'help':
            'Load reference complex structure. This option allows for comparison with the reference complex structure '
            'and triggers additional analysis features.\n\n'
            'REF must be either:\n\n'
            '[1] [pdb code]:[protein chains]:[peptide1 chain][peptide2 chain]...\n'
            '[2] [pdb file]:[protein chains]:[peptide1 chain][peptide2 chain]...\n\n'
            'i.e 1abc:AB:C, 1abc:AB:CD, myfile.pdb:AB:C, myfile.pdb.gz:AB:CDE'
    },
    'replicas': {
        'flag': '-r',
        'default': 10,
        'type': int,
        'metavar': 'NUM',
        'help':
            'Set number of replicas to be used in Replica Exchange Monte Carlo. (NUM > 0, default: %(default)s)'
    },
    'replicas-dtemp': {
        'flag': '-D',
        'default': 0.5,
        'type': float,
        'metavar': 'DELTA',
        'help': 'Set temperature increment between replicas. (DELTA > 0, default: %(default)s)'
    },
    'save-cabs-files': {
        'flag': '-S',
        'action': 'store_true',
        'help': 'Save CABSdock simulation files. (default: %(default)s)'
    },
    'save-config': {
        'flag': '-C',
        'action': 'store_true',
        'help': 'Save simulation parameters in config file. (default: %(default)s)'
    },
    'sc-rest-add':{
        'action': 'append',
        'nargs': 4,
        'metavar': ('RESI', 'RESJ', 'DIST', 'WEIGHT'),
        'help':
            'Add distance restraint between SC pseudoatom in residue RESI and SC pseudoatom in residue RESJ.\n'
            'For more details see help for \'--ca-rest-add\''
    },
    'sc-rest-file': {
        'action': 'append',
        'metavar': 'FILE',
        'help': 'Read SC restraints from file (use multiple times to add multiple files)'
    },
    'sc-rest-weight': {
        'default': 1.0,
        'type': float,
        'metavar': 'WEIGHT',
        'help': 'Set global weight for all SC restraints (default: %(default)s)'
    },
    'separation': {
        'flag': '-d',
        'default': 20.0,
        'type': float,
        'metavar': 'SEP',
        'help':
            'This option enables advanced settings of building starting conformations of modelled complexes (to be '
            'used only in specific protocols). The option sets separation distance in Angstroms between the peptide '
            'and the surface of the protein (default: %(default)s Angstroms)'
    },
    'temperature': {
        'flag': '-t',
        'default': (2.0, 1.0),
        'nargs': 2,
        'type': float,
        'metavar': ('TINIT', 'TFINAL'),
        'help':
            'Set temperature range for simulated annealing.\n'
            'TINIT - initial temperature, TFINAL - final temperature (default: (TINIT, TFINAL) = %(default)s.\n'
            'CABSdock uses temperature units, which do not correspond straightforwardly to real temperatures. '
            'Temperature around 1.0 roughly corresponds to nearly frozen conformation, folding temperature of a small '
            'proteins in the CABS model is usually around 2.0'
    },
    'verbose': {
        'flag': '-v',
        'choices': [-1, 0, 1, 2, 3],
        'default': 1,
        'type': int,
        'metavar': 'LEVEL',
        'help': 'Set verbosity LEVEL: (default: %(default)s)\n'
                '-1 - silent mode (only critical messages)\n'
                ' 3 - maximum verbosity'
    },
    'version': {
        'action': 'store_true',
        'help': 'print version and exit program'
    },
    'weighted-fit': {
        'action': 'store_true',
        'help':
            'This option fixes a problem occurring when part of the protein is assigned increased flexibility, '
            'which in default mode results in sometimes incorrect structural fit between generated models, which in '
            'consequence affects clustering and scoring. Structural fit can be corrected by using weights assigned '
            'to the protein residues, which modify the significance of a given residue in the structural fit process. '
            'By default all residues are assigned uniform weight w = 1.0.\n\n'
            'Use \'--weighted-fit\' in combination with \'--receptor-flexibility\' to use flexibility settings as '
            'weights for structural fit\n\nOR\n\n'
            'use \'--weighted-fit-file <path to file>\' to provide custom weights.\n'
            'The file should follow this format:\n\n'
            'default 1.0 (default value, if omitted w = 1.0 is assumed)\n1:A 0.5\n2:A 0.1\n...\n1:B 0.99\netc.\n\n'
            'Weights should be numbers from [0:1] with \'0\' meaning irrelevant in fitting process.'
    },
    'weighted-fit-file': {
        'metavar': 'FILE',
        'help': 'See help entry for \'--weighted-fit\''
    },
    'work-dir': {
        'flag': '-w',
        'default': '.',
        'metavar': 'DIR',
        'help': 'Set working directory to DIR. (default: \'%(default)s\')'
    },
}


DockParser = mk_parser(dock_dict, dc(groups), dc(options))
FlexParser = mk_parser(flex_dict, dc(groups), dc(options))
