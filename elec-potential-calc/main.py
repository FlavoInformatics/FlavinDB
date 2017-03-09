import argparse
import potential_calc
import os
import sys

from itertools import chain


def main(arguments):
    codes = arguments.pdb_codes

    if arguments.codes_file:
        with open(arguments.codes_file) as f:
            lines = [line.strip().split() for line in list(f)]
            codes += list(chain(*lines))

    potential_calc.get_potentials(codes, arguments.output_path)


if __name__ == "__main__":
    # set up command line arguments
    parser = argparse.ArgumentParser(description="Compute electrostatic potential for a given set of proteins.")
    parser.add_argument("--codes", "-c", dest="pdb_codes", type=str, help="Four letter PDB protein identifiers.", nargs="+", default=[])
    parser.add_argument("--output", "-o", dest="output_path", type=str, help="Target path for output files", default="")
    parser.add_argument("--get-codes-from-file", dest="codes_file", type=str, help="Path to whitespace delineated PDB identifiers on disk", default=None)
    args = parser.parse_args()

    # make sure that at least one PDB file is specified by the user
    if not args.pdb_codes and not args.codes_file:
        parser.print_help()
        sys.exit(os.EX_USAGE)

    main(args)
