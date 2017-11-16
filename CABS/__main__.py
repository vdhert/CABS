import re
import sys
import argparse
import imp
import traceback as _tr
from shutil import rmtree

try:
    from CABS import optparser, logger, __version__, _JUNK
except ImportError:
    cabs_module = imp.find_module('CABS', ['.'])
    imp.load_module('CABS', *cabs_module)
    from CABS import optparser, logger, __version__, _JUNK


def run(cabs_cmd, cmd_line=sys.argv[1:]):
    if cabs_cmd not in ['dock', 'flex']:
        raise IndexError

    module_name = 'CABS' + cabs_cmd
    parser = getattr(optparser, cabs_cmd + '_parser')

    pre_parser = argparse.ArgumentParser(add_help=False)
    pre_parser.add_argument('-c', '--config')
    pre_parser.add_argument('--version', action='store_true')
    pre_parser.add_argument('-h', '--help', action='store_true')

    pre_args, remains = pre_parser.parse_known_args(cmd_line)
    if pre_args.help:
        _help = parser.format_help()
        print re.sub('\n( *)\n( *)\n', '\n\n', _help)
        sys.exit(0)
    elif pre_args.version:
        print __version__
        sys.exit(0)
    elif pre_args.config:
        remains = optparser.ConfigFileParser(pre_args.config).args + remains

    config = vars(parser.parse_args(remains))
    import CABS.job
    task = getattr(CABS.job, cabs_cmd.title() + 'Task')
    job = task(**config)

    try:
        job.run()
    except KeyboardInterrupt:
        logger.critical(module_name, 'Interrupted by user.')
    except Exception as e:
        logger.exit_program(module_name, 'Error occurred', traceback=_tr.format_exc(), exc=e)
    finally:
        logger.close_log()
        for _file in _JUNK:
            rmtree(_file, ignore_errors=True)


def run_dock(cmd_line=sys.argv[1:]):
    run('dock', cmd_line)


def run_flex(cmd_line=sys.argv[1:]):
    run('flex', cmd_line)


if __name__ == '__main__':
    try:
        run(sys.argv[1], sys.argv[2:])

    except IndexError:
        print 'usage: python CABS <cmd> <options>\n\tcmd: dock or flex.\n\t' \
              'For the list of <options> run \'python CABS <cmd> -h\''
        sys.exit(0)
