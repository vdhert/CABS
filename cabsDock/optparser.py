import argparse


class ParserFactory:
    def __init__(self, parser_config):
        self.parser = argparse.ArgumentParser(
            **parser_config
        )
        self.groups = {}

    def add_group(self, name, title, desc):
        group = self.parser.add_argument_group(
            title=title,
            description=desc
        )
        self.groups[name] = group

    def add_option(self, option, groupname):
        self.groups[groupname].add_argument(**option)

