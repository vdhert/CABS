"""
CABSdock is a ...
"""


import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='CABSdock 2.0', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-r', '--receptor', type=str, required=True, metavar='RECEPTOR_INPUT')
    parser.add_argument('-c', '--config', type=str, metavar='CONFIG_FILE')
    parser.add_argument('-d', '--dir', type=str, default='.', metavar='PATH')
    parser.add_argument('-l', '--add_ligand', action='append', type=str, nargs='+', metavar='LIGAND_INPUT')

    args = parser.parse_args()
    print args
