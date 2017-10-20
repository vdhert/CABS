import os
import re
import sys
import argparse
from CABS import logger, __version__


def run_dock(cmd_line=sys.argv[1:]):

    junk = []  # put here filepaths to whatever should be deleted if cabs crashes
    from CABS.optparser import DockParser as parser, ConfigFileParser

    preparser = argparse.ArgumentParser(add_help=False)
    preparser.add_argument('-c', '--config')
    preparser.add_argument('--version', action='store_true')
    preparser.add_argument('-h', '--help', action='store_true')

    preargs, remains = preparser.parse_known_args(cmd_line)
    if preargs.help:
        _help = parser.format_help()
        print re.sub("\n( *)\n( *)\n", "\n\n", _help)
        sys.exit(0)
    elif preargs.version:
        print __version__
        sys.exit(0)
    elif preargs.config:
        remains = ConfigFileParser(preargs.config).args + remains

    config = vars(parser.parse_args(remains))

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


def run_flex(cmd_line=sys.argv[1:]):

    junk = []  # put here filepaths to whatever should be deleted if cabs crashes
    from CABS.optparser import FlexParser as parser, ConfigFileParser

    preparser = argparse.ArgumentParser(add_help=False)
    preparser.add_argument('-c', '--config')
    preparser.add_argument('--version', action='store_true')
    preparser.add_argument('-h', '--help', action='store_true')

    preargs, remains = preparser.parse_known_args(cmd_line)
    if preargs.help:
        _help = parser.format_help()
        print re.sub("\n( *)\n( *)\n", "\n\n", _help)
        sys.exit(0)
    elif preargs.version:
        print __version__
        sys.exit(0)
    elif preargs.config:
        remains = ConfigFileParser(preargs.config).args + remains

    config = vars(parser.parse_args(remains))

    from CABS.job import FlexTask
    job = FlexTask(**config)

    # start docking
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


if __name__ == '__main__':
    try:
        cmd = sys.argv[1]
        options = sys.argv[2:]

        if cmd == 'dock':
            run_dock(options)
        elif cmd == 'flex':
            run_flex(options)
        else:
            raise IndexError

    except IndexError:
        print 'usage: python CABS <cmd> <options>\n\tcmd: dock or flex.\n\t' \
              'For the list of <options> run \'python CABS <cmd> -h\''
        sys.exit(0)
