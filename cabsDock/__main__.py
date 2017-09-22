import logger
import os
from cabsDock.optparser import ParserFactory, ConfigFileParser
from pkg_resources import resource_filename
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


def run_dock():

    junk = []  # put here filepaths to whatever should be deleted if cabs crashes

    parser = ParserFactory(
        filecsv=resource_filename('cabsDock', 'data/data3.dat')
    ).parser
    args = parser.parse_args()

    cfg_args = []
    if args.config:
        cfg_args = ConfigFileParser(args.config).args

    parser = ParserFactory(
        filecsv=resource_filename('cabsDock', 'data/data3.dat'), required=['receptor']
    ).parser

    args = parser.parse_args(cfg_args + argv[1:])
    config = Config(args)

    from cabsDock.job import DockTask
    job = DockTask(**config)

    # start docking
    try:
        job.run()
    except KeyboardInterrupt:
        logger.info(
            module_name='CABSdock',
            msg='Interrupted by user.'
        )
    except Exception as e:
        logger.exit_program(
            module_name='CABSdock',
            msg=e.message,
            exc=e,
            traceback=(logger.log_level > 2)
        )
    finally:
        map(os.removedirs,junk)


def run_flex():

    junk = []  # put here filepaths to whatever should be deleted if cabs crashes

    parser = ParserFactory(
        filecsv=resource_filename('cabsDock', 'data/data4.dat')
    ).parser
    args = parser.parse_args()

    cfg_args = []
    if args.config:
        cfg_args = ConfigFileParser(args.config).args

    parser = ParserFactory(
        filecsv=resource_filename('cabsDock', 'data/data4.dat'), required=['structure']
    ).parser

    args = parser.parse_args(cfg_args + argv[1:])
    config = dict(vars(args))

    from cabsDock.job import FlexTask
    job = FlexTask(**config)

    # start flexing
    try:
        job.run()
    except KeyboardInterrupt:
        logger.info(
            module_name='CABSflex',
            msg='Interrupted by user.'
        )
    except Exception as e:
        logger.exit_program(
            module_name='CABSflex',
            msg=e.message,
            exc=e,
            traceback=(logger.log_level > 2)
        )
    finally:
        map(os.removedirs, junk)