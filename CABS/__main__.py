import os
import re
import argparse
from sys import argv, exit
from CABS import logger, __version__


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
    from CABS.optparser import DockParser as parser, ConfigFileParser

    preparser = argparse.ArgumentParser(add_help=False)
    preparser.add_argument('-c', '--config')
    preparser.add_argument('--version', action='store_true')
    preparser.add_argument('-h', '--help', action='store_true')

    preargs, remains = preparser.parse_known_args()
    if preargs.help:
        _help = parser.format_help()
        print re.sub("\n( *)\n( *)\n", "\n\n", _help)
        exit(0)
    elif preargs.version:
        print __version__
        exit(0)
    elif preargs.config:
        remains = ConfigFileParser(preargs.config).args + remains

    config = Config(parser.parse_args(remains))

    from CABS.job import DockTask
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
        map(os.removedirs, junk)


def run_flex():
    junk = []  # put here filepaths to whatever should be deleted if cabs crashes
    from CABS.optparser import FlexParser as parser, ConfigFileParser

    preparser = argparse.ArgumentParser(add_help=False)
    preparser.add_argument('-c', '--config')
    preparser.add_argument('--version', action='store_true')
    preparser.add_argument('-h', '--help', action='store_true')

    preargs, remains = preparser.parse_known_args()
    if preargs.help:
        _help = parser.format_help()
        print re.sub("\n( *)\n( *)\n", "\n\n", _help)
        exit(0)
    elif preargs.version:
        print __version__
        exit(0)
    elif preargs.config:
        remains = ConfigFileParser(preargs.config).args + remains

    config = Config(parser.parse_args(remains))

    from CABS.job import FlexTask
    job = FlexTask(**config)

    # start flexing
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
        map(os.removedirs, junk)
