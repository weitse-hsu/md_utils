import os
import sys
import time
import argparse
from general_utils import utils
from md_utils.structure import protein
from pymol import cmd


def initialize(args):
    parser = argparse.ArgumentParser(
        description="Identify corresponding residues composing secondary structure domains in a given \
            protein structure."
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input PDB/GRO file containing the protein structure."
    )
    parser.add_argument(
        "-s",
        "--ss_type",
        type=str,
        choices=["H", "S", "T"],
        required=True,
        help="Type of secondary structure domain to identify: helix, sheet, or turn."
    )
    parser.add_argument(
        "-n",
        "--min_length",
        type=int,
        default=5,
        help="Minimum length of the secondary structure domain to consider. The default is 5 residues."
    )
    parser.add_argument(
        "-l",
        "--log",
        type=str,
        default="identify_ss_domains.log",
        help="Log file to record the output. The default is identify_ss_domains.log."
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

    cmd.reinitialize()
    cmd.load(args.input, "structure")
    model = cmd.get_model(f"structure and ss {args.ss_type} and name CA")  # Use CA atoms to get one entry per residue

    segment_dict = {}  # maps chain ID to list of (start_res, end_res) segment tuples
    for i, a in enumerate(model.atom):
        res = f"{protein.convert_res_code(a.resn)}{a.resi}"
        chain = a.chain
        if chain not in segment_dict:
            segment_dict[chain] = []
        if i == 0:
            start_res = res
            prev_resi = int(a.resi)
        else:
            curr_resi = int(a.resi)
            if curr_resi == prev_resi + 1:
                prev_resi = curr_resi
            else:
                # New segment
                end_res = f"{protein.convert_res_code(model.atom[i-1].resn)}{model.atom[i-1].resi}"
                segment_dict[model.atom[i-1].chain].append((start_res, end_res))
                start_res = res
                prev_resi = curr_resi

    end_res = f"{protein.convert_res_code(model.atom[-1].resn)}{model.atom[-1].resi}"
    segment_dict[chain].append((start_res, end_res))

    # Filter segments by minimum length
    for chain in segment_dict:
        filtered_segments = []
        for start, end in segment_dict[chain]:
            start_idx = int(''.join(filter(str.isdigit, start)))
            end_idx = int(''.join(filter(str.isdigit, end)))
            if (end_idx - start_idx + 1) >= args.min_length:
                filtered_segments.append((start, end))
        segment_dict[chain] = filtered_segments

    for chain, segments in segment_dict.items():
        print(f"Chain {chain}:")
        for start, end in segments:
            print(f"  {start} to {end}")

    print(segment_dict)

    # Generate PyMol selection strings (name objects as H_1, S_1, T_1, etc.)
    for chain, segments in segment_dict.items():
        for idx, (start, end) in enumerate(segments, start=1):
            selection_name = f"{args.ss_type}_{chain}_{idx}"
            pymol_selection = f"chain {chain} and resi {start[1:]}-{end[1:]}"
            print(f"select {selection_name}, {pymol_selection}")

    print(f"Elapsed time: {utils.format_time(time.time() - t1)}")
