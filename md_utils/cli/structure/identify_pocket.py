import os
import sys
import time
import argparse
import warnings
import MDAnalysis as mda
from general_utils import utils


def initialize(args):
    parser = argparse.ArgumentParser(
        description="Identify pocket residues around a ligand in a given binding complex structure."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=False,
        help="Input PDB/GRO file containing the complex structure."
    )
    parser.add_argument(
        "-r",
        "--resname",
        type=str,
        default="LIG",
        help="Residue name of the ligand. The default is 'LIG'."
    )
    parser.add_argument(
        "-c",
        "--cutoff",
        type=float,
        default=6.0,
        help="Cutoff distance in Å to define the pocket. The default is 6.0 Å."
    )
    parser.add_argument(
        "-l",
        "--log",
        type=str,
        default="identify_pocket.log",
        help="Log file to record the output. The default is identify_pocket.log."
    )

    args = parser.parse_args()

    return args


def main():
    t1 = time.time()
    args = initialize(sys.argv[1:])
    sys.stdout = utils.Logger(args.log)
    sys.stderr = utils.Logger(args.log)

    print(f"\nCommand line: {' '.join(sys.argv)}")
    print(f"Current working directory: {os.getcwd()}\n")

    u = mda.Universe(args.input)

    ligand = u.select_atoms(f"resname {args.resname}")
    pocket = u.select_atoms(f"protein and around {args.cutoff} group ligand", ligand=ligand)
    if len(ligand.residues) == 0:
        warnings.warn(
            "No pocket identification will be performed as no residues with the specified "
            "ligand resname were found. Please check the input file and ligand resname."
        )
    else:
        if len(pocket.residues) == 0:
            print("No pocket residues found within the specified cutoff distance.")
        else:
            print(f"List of pocket residues ({len(pocket.residues)}):")
            for i in pocket.residues:
                print(f"  {i.resname} {i.resid}")

            pymol_selection = f"select pocket, resi {'+'.join([str(i.resid) for i in pocket.residues])}"
            residues = " ".join([str(i.resid) for i in pocket.residues])
            ndx_selection = f'a N CA C O & r {residues}'
            vmd_selection = f"resid {' '.join([str(i.resid) for i in pocket.residues])}"

            print(f"\n- PyMOL selection:\n{pymol_selection}")
            print(f"\n- GROMACS index file selection:\n{ndx_selection}")
            print(f"\n- VMD selection:\n{vmd_selection}")

    print(f"\nElapsed time: {utils.format_time(time.time() - t1)} seconds")
