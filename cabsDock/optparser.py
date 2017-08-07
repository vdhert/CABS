import csv
import argparse


class CustomFormatter(argparse.RawTextHelpFormatter):
    def __init__(self, *args, **kwargs):
        argparse.HelpFormatter.__init__(self, indent_increment=1, max_help_position=4, *args, **kwargs)
    # TODO zawijanie dlugich lini, indent dla list


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

    FIELDS = {}

    def __init__(self, filecsv, fields, sep):

        ParserFactory.FIELDS['name'] = fields[0] - 1
        ParserFactory.FIELDS['usage'] = fields[1] - 1
        ParserFactory.FIELDS['help'] = fields[2] - 1
        ParserFactory.FIELDS['default'] = fields[3] - 1
        ParserFactory.FIELDS['action'] = fields[4] - 1

        self.defaults = {}

        with open(filecsv) as f:
            lines = [l for l in csv.reader(f, delimiter=sep)]
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

        name = [n.strip() for n in line[ParserFactory.FIELDS['name']].split(',')]
        usage = line[ParserFactory.FIELDS['usage']].split()
        help = line[ParserFactory.FIELDS['help']]
        default = line[ParserFactory.FIELDS['default']].split()
        action = line[ParserFactory.FIELDS['action']]

        kwargs = {}
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
                kwargs['type'] = type(kwargs['default'])
        else:
            pass

        kwargs['help'] = help
        if action:
            kwargs['action'] = action

        _group.add_argument(*name, **kwargs)

if __name__ == '__main__':
   pass
