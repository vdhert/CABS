import csv
import argparse
import re
from pkg_resources import require

__version__ = require('CABS')[0].version


class CustomFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, *args, **kwargs):
        argparse.HelpFormatter.__init__(self, indent_increment=1, max_help_position=4, *args, **kwargs)


def string_cast(string):
    try:
        var = int(string)
        return var
    except ValueError:
        try:
            var = float(string)
            return var
        except ValueError:
            s = string.lower()
            if s == 'true':
                return True
            elif s == 'false':
                return False
            else:
                return string


class ParserFactory:

    def __init__(self, filecsv, required=()):

        self.required = required

        with open(filecsv) as f:
            lines = [l for l in csv.reader(f, delimiter=',')]
        p, g = ParserFactory._parse_csv(lines)
        self.parser = ParserFactory._build_parser(p)
        self._populate_parser(g, self.parser)

    @staticmethod
    def _parse_csv(lines):
        objects = []
        current = None
        for line in lines:
            if not any(line):
                pass
            elif line[0].startswith('#'):
                if current:
                    objects.append(current)
                    current = [line]
                else:
                    current = [line]
            else:
                if current:
                    current.append(line)
        objects.append(current)

        parser = [o for o in objects if o[0][0] == '#parser']
        groups = [o for o in objects if o[0][0] == '#group']
        if len(parser) != 1:
            raise Exception('Multiple \'#parser\' directives in csv file.')
        elif len(groups) == 0:
            raise Exception('No groups found in csv file.')
        return parser[0], groups

    @staticmethod
    def _build_parser(lines):
        d = {l[0]: l[1] for l in lines[1:]}
        d['version'] = __version__
        return argparse.ArgumentParser(formatter_class=CustomFormatter, **d)

    def _populate_parser(self, _groups, _parser):
        for g in _groups:
            self._add_group(g, _parser)

    def _add_group(self, lines, _parser):
        group = _parser.add_argument_group(
            title=lines[0][1],
            description=lines[0][2]
        )
        for line in lines[1:]:
            self._add_arg(line, group)

    def _add_arg(self, line, _group):

        name = [n.strip() for n in line[0].split(',')]
        usage = line[1].split()
        help = line[2]
        default = line[3].split()
        action = line[4]
        argtype = line[5]

        kwargs = {}
        if name[-1][2:] in self.required:
            kwargs['required'] = True

        ile = len(usage)
        if ile > 1:
            kwargs['nargs'] = ile
            kwargs['metavar'] = tuple(usage)
            if len(default) == ile:
                kwargs['default'] = tuple(string_cast(s) for s in default)
        elif ile == 1:
            kwargs['metavar'] = usage[0]
            if len(default) == 1:
                kwargs['default'] = string_cast(default[0])
        else:
            pass
        if any('=' in i for i in usage):
            kwargs['nargs'] = '+'

        kwargs['help'] = help
        if action:
            kwargs['action'] = action

        if action == 'count':
            kwargs.pop('metavar')

        if argtype == 'int':
            kwargs['type'] = int
        elif argtype == 'float':
            kwargs['type'] = float
        elif 'lambda' in argtype:
            kwargs["type"] = lambda x: x.split('=')

        _group.add_argument(*name, **kwargs)


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


if __name__ == '__main__':
    pass
