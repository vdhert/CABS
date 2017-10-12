import os
import re
from sys import argv
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
    from CABS.optparser import DockParser as parser
    from CABS.optparser import preparse

    config = Config(
        parser.parse_args(
            preparse(argv[1:])
        )
    )

    if config['help']:
        _help = parser.format_help()
        print re.sub("\n( *)\n( *)\n", "\n\n", _help)
        exit(0)
    elif config['version']:
        print __version__
        exit(0)
    else:
        print config
        exit(0)

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

    # parser = mk_flex_parser()
    #
    # config = Config(
    #     parser.parse_args(
    #         collect_args(argv)
    #     )
    # )

    from CABS.job import FlexTask
    job = FlexTask(**config)

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


if __name__ == '__main__':
    run_dock()