from optparser import ParserFactory
from pkg_resources import resource_filename


def run_job():
    args = ParserFactory(
        filecsv=resource_filename('cabsDock', 'data/data3.dat'),
        fields=(1, 2, 3, 4, 5),
        sep=','
    ).parser.parse_args()

if __name__ == '__main__':
    run_job()
