from cabsDock.optparser import ParserFactory
from pkg_resources import resource_filename
from ConfigParser import SafeConfigParser
from StringIO import StringIO
from sys import argv


class Config(dict):
    def __init__(self, config):
        """
        Smart dictionary that reads argparse.Namespace returned by argparse.parse_args().
        also checks / updates some of the parsed options
        :param config: argparse.Namespace
        """

        # config.peptide + config.add_peptide -> config.ligand
        config.ligand = []
        if config.peptide:
            config.ligand.extend([[p, 'random', 'random'] for p in config.peptide])
        if config.add_peptide:
            config.ligand.extend([p for p in config.add_peptide if p])

        dict.__init__(self, vars(config))

    def __repr__(self):
        return '\n'.join([k + ': ' + str(v) for k, v in sorted(self.items())])


def run_job():

    parser = ParserFactory(
        filecsv=resource_filename('cabsDock', 'data/data3.dat')
    ).parser
    args = parser.parse_args()

    cfg_args = []
    if args.config:
        with open(args.config) as f:
            f = StringIO('[config]\n' + f.read())
            config_parser = SafeConfigParser(allow_no_value=True)
            # config_parser.optionxform = lambda x: x
            config_parser.readfp(f)
            for k, v in config_parser.items('config'):
                cfg_args.append('--' + k)
                cfg_args.extend(v.split())

    parser = ParserFactory(
        filecsv=resource_filename('cabsDock', 'data/data3.dat'), required=['receptor']
    ).parser

    args = parser.parse_args(cfg_args + argv[1:])
    config = Config(args)

    from cabsDock.job import Job
    job = Job(**config)

    # tutaj dopisac warunek na cabsflexa

    # wypisz config
    # print job

    # start docking
    job.cabsdock()

    
if __name__ == '__main__':
    run_job()
