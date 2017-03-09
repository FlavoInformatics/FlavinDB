import argparse
import potential_calc
import os
import sys

from itertools import chain


def main(arguments):
    # make sure that at least one protein is specified by the user
    if not arguments.pdb_codes and not arguments.codes_file:
        print "Please specify the four letter PDB identifier for the protein(s), either on the command line with the" \
            "-c' arg, or as a file with the --get-codes-from-file arg"
        sys.exit(os.EX_USAGE)

    pdb_codes = arguments.pdb_codes

    if arguments.codes_file:
        with open(arguments.codes_file) as f:
            lines = [line.strip().split() for line in list(f)]
            pdb_codes += list(chain(*lines))

    potential_calc.get_potentials(pdb_codes, arguments.output_path)


if __name__ == "__main__":
    # set up command line arguments
    parser = argparse.ArgumentParser(description="Compute electrostatic potential for a given set of proteins.")
    parser.add_argument("--codes", "-c", dest="pdb_codes", type=str, help="Four letter PDB protein identifiers.", nargs="+", default=[])
    parser.add_argument("--output", "-o", dest="output_path", type=str, help="Target path for output files", default=".")
    parser.add_argument("--get-codes-from-file", dest="codes_file", type=str, help="Path to whitespace delineated PDB codes on disk", default=None)
    args = parser.parse_args()


    main(args)
