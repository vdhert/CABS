import re
import sys
import argparse
import imp
import traceback as _tr

try:
    from CABS import logger, __version__, _JUNK
except ImportError:
    cabs_module = imp.find_module("CABS", ["."])
    imp.load_module("CABS", *cabs_module)
    from CABS import logger, __version__, _JUNK

from shutil import rmtree

# TODO:     run_dock and run_flex became long enough, to write run_task(command, argv) and change run_dock
# TODO:     and run_flex to run_dock(argv): run_task('dock', argv) etc.


def run_dock(cmd_line=sys.argv[1:]):

    from CABS.optparser import DockParser as Parser, ConfigFileParser

    preparser = argparse.ArgumentParser(add_help=False)
    preparser.add_argument('-c', '--config')
    preparser.add_argument('--version', action='store_true')
    preparser.add_argument('-h', '--help', action='store_true')

    preargs, remains = preparser.parse_known_args(cmd_line)
    if preargs.help:
        _help = Parser.format_help()
        print re.sub("\n( *)\n( *)\n", "\n\n", _help)
        sys.exit(0)
    elif preargs.version:
        print __version__
        sys.exit(0)
    elif preargs.config:
        remains = ConfigFileParser(preargs.config).args + remains

    config = vars(Parser.parse_args(remains))

    from CABS.job import DockTask
    job = DockTask(**config)

    # start docking
    try:
        job.run()
    except KeyboardInterrupt:
        logger.info('CABSdock', 'Interrupted by user.')
    except Exception as e:
        logger.exit_program(module_name='CABSdock', msg='Error occured', traceback=_tr.format_exc(), exc=e)
    finally:
        for _file in _JUNK:
            rmtree(_file, ignore_errors=True)


def run_flex(cmd_line=sys.argv[1:]):
    from CABS.optparser import FlexParser as Parser, ConfigFileParser

    preparser = argparse.ArgumentParser(add_help=False)
    preparser.add_argument('-c', '--config')
    preparser.add_argument('--version', action='store_true')
    preparser.add_argument('-h', '--help', action='store_true')

    preargs, remains = preparser.parse_known_args(cmd_line)
    if preargs.help:
        _help = Parser.format_help()
        print re.sub("\n( *)\n( *)\n", "\n\n", _help)
        sys.exit(0)
    elif preargs.version:
        print __version__
        sys.exit(0)
    elif preargs.config:
        remains = ConfigFileParser(preargs.config).args + remains

    config = vars(Parser.parse_args(remains))

    from CABS.job import FlexTask
    job = FlexTask(**config)

    # start flexing
    try:
        job.run()
    except KeyboardInterrupt:
        logger.info('CABSflex', 'Interrupted by user.')
    except Exception as e:
        logger.exit_program(module_name='CABSflex', msg='Error occured', traceback=_tr.format_exc(), exc=e)
    finally:
        for _file in _JUNK:
            rmtree(_file, ignore_errors=True)


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
