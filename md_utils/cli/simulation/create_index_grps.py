import os
import sys
import tempfile
import argparse
from general_utils import utils
from md_utils.simulation import gmx_utils, gmx_parser


def initialize(args):
    parser = argparse.ArgumentParser(
        description="Create GROMACS index groups given a text file containing selections."
    )
    parser.add_argument(
        "-f",
        "--gro",
        type=str,
        required=True,
        help="Path to the input GROMACS .gro file."
    )
    parser.add_argument(
        "-s",
        "--selections",
        type=str,
        default="selections.txt",
        help="Path to the input text file containing selections. Each line in the file should contain a \
            GROMACS selection string commented (#) with the desired group name. For example:\
            '1 | 13 # ICL1'. Prepend '#' to comment out any lines that should be ignored. The default is\
            'selections.txt'."
    )
    parser.add_argument(
        "-n",
        "--ndx",
        type=str,
        help="Path to the input GROMACS index (.ndx) file. If not provided, a new index file will be created \
            from the .gro file."
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="index.ndx",
        help="Path to the output GROMACS index (.ndx) file. The default is 'index.ndx'."
    )
    parser.add_argument(
        "-l",
        "--log",
        type=str,
        default="create_index_grps.log",
        help="Log file to record the output of the index group creation steps. The default is create_index_grps.log."
    )
    args = parser.parse_args(args)

    return args


def main():
    args = initialize(sys.argv[1:])
    sys.stdout = utils.Logger(args.log)
    sys.stderr = utils.Logger(args.log)

    print(f"\nCommand line: {' '.join(sys.argv)}")
    print(f"Current working directory: {os.getcwd()}")

    # Test GROMACS installation
    gmx_utils.run_gmx_cmd(['gmx', '--version'], print_output=False)

    gmx_args = [
        'gmx', 'make_ndx',
        '-f', args.gro,
        '-o', args.output
    ]

    if args.ndx:
        gmx_args.extend(['-n', args.ndx])

    # Count the currently available number of groups
    if args.ndx:
        groups, group_str = gmx_parser.parse_ndx(args.ndx)
    else:
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_ndx = os.path.join(temp_dir, "temp.ndx")
            gmx_utils.run_gmx_cmd([
                'gmx', 'make_ndx',
                '-f', args.gro,
                '-o', temp_ndx
            ], print_output=False, prompt_input='q\n')
            groups, group_str = gmx_parser.parse_ndx(temp_ndx)

    n_groups = len(groups)  # Current number of groups
    print(f"Current index groups ({n_groups}):")
    print(group_str)

    print(f"\nNow creating GROMACS index groups from selections in {args.selections}...")
    with open(args.selections, 'r') as f:
        selections = f.readlines()

    prompt_input = ''
    for selection in selections:
        if not selection.strip() or selection.strip().startswith('#'):
            print(f"Skipping line: {selection.strip()}")
            continue
        selection = selection.strip()
        if not selection or selection.startswith('#'):
            continue
        sel_str, grp_name = selection.split('#')
        sel_str = sel_str.strip()
        grp_name = grp_name.strip()
        if grp_name in groups:
            print(f"Skipping group '{grp_name}': name already exists in index file.")
            continue
        n_groups += 1
        prompt_input += f"{sel_str}\nname {n_groups-1} {grp_name}\n"  # Note that the index starts from 0

    prompt_input += 'q\n'
    print(f"Prompt input for make_ndx:\n{prompt_input}")
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args, prompt_input=prompt_input, print_output=False)

    new_groups, new_group_str = gmx_parser.parse_ndx(args.output)
    print(f"Current index groups ({len(new_groups)}):")
    print(new_group_str)
