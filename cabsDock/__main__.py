from cabsDock.optparser import ParserFactory
from pkg_resources import resource_filename
from ConfigParser import ConfigParser
from StringIO import StringIO
from cabsDock.job import Job


class Config(dict):
    def __init__(self, cl_dict, cfg_dict={}):
        config = cfg_dict

        dict.__init__(self, **config)


def run_job():
    args = ParserFactory(
        filecsv=resource_filename('cabsDock', 'data/data3.dat'),
        fields=(1, 2, 3, 4, 5),
        sep=','
    ).parser.parse_args()

    if args.config:
        with open(args.config) as f:
            f = StringIO('[config]\n' + f.read())
            config_parser = ConfigParser()
            config_parser.readfp(f)
            cfg_args = dict(config_parser.items('config'))
    else:
        cfg_args = {}

    config = Config(
        cl_dict=vars(args),
        cfg_dict=cfg_args
    )

    for k in config:
        print k + ':', config[k]
    # job = Job(**config)
    # job.run_job()


if __name__ == '__main__':
    run_job()
