import csv
import os
import sys
import time
import argparse
from general_utils import utils
from md_utils.structure import protein


def _parse_component(value: str):
    """Parse 'NAME:ID1,ID2,...' into (name, [ids])."""
    if ":" not in value:
        raise argparse.ArgumentTypeError(
            f"Component must be formatted as 'NAME:ID1,ID2,...', got: {value!r}"
        )
    name, ids_str = value.split(":", 1)
    name = name.strip()
    ids = [s.strip() for s in ids_str.split(",") if s.strip()]
    if not name:
        raise argparse.ArgumentTypeError("Component name cannot be empty.")
    if not ids:
        raise argparse.ArgumentTypeError(f"No IDs provided for component {name!r}.")
    return name, ids


def _parse_landmark(value: str):
    """Parse 'COMPONENT:RESID1,RESID2,...' into (component, [resids]).

    Each token may be a single residue number or a range (e.g. 118-132).
    """
    if ":" not in value:
        raise argparse.ArgumentTypeError(
            f"Landmark must be formatted as 'COMPONENT:RESID1,RESID2,...', got: {value!r}"
        )
    name, resids_str = value.split(":", 1)
    name = name.strip()
    tokens = []
    try:
        for token in resids_str.split(","):
            token = token.strip()
            if not token:
                continue
            if not token.startswith("-") and "-" in token:
                start_str, end_str = token.split("-", 1)
                start, end = int(start_str.strip()), int(end_str.strip())
                if start > end:
                    raise argparse.ArgumentTypeError(
                        f"Invalid range {token!r} in landmark {value!r}: start > end"
                    )
                tokens.append((start, end))
            else:
                tokens.append(int(token))
    except argparse.ArgumentTypeError:
        raise
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid residue number in landmark {value!r}: {exc}"
        ) from exc
    return name, tokens


def initialize(args):
    parser = argparse.ArgumentParser(
        description=(
            "Map residue numbers from a reference PDB to a simulation GRO file using "
            "a fitting Needleman-Wunsch sequence alignment. Each molecular component "
            "(e.g. a protein chain) is aligned independently against the full GRO "
            "protein sequence."
        )
    )
    parser.add_argument(
        "-p", "--pdb",
        type=str,
        required=True,
        help="Reference PDB file with correct residue numbers."
    )
    parser.add_argument(
        "-g", "--gro",
        type=str,
        required=True,
        help="Simulation GRO file to map residues onto."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="residue_mappings.tsv",
        help="Output TSV file for the mapped residue pairs. Default: residue_mappings.tsv."
    )
    parser.add_argument(
        "-c", "--component",
        type=_parse_component,
        action="append",
        dest="components",
        required=True,
        metavar="NAME:ID[,ID...]",
        help=(
            "Define one molecular component and the PDB segment IDs (or chain IDs, see "
            "--pdb-id-field) that belong to it. Repeat for each component. "
            "Example: -c ATP8B:PROA,PROC -c CDC50A:PROB"
        )
    )
    parser.add_argument(
        "--pdb-id-field",
        choices=["segid", "chain"],
        default="segid",
        help=(
            "PDB field used to identify components: 'segid' reads columns 73-76 "
            "(standard segment ID), 'chain' reads column 22 (chain ID). Default: segid."
        )
    )
    parser.add_argument(
        "--landmark",
        type=_parse_landmark,
        action="append",
        dest="landmarks",
        default=[],
        metavar="COMPONENT:RESID[,RESID...|START-END...]",
        help=(
            "PDB residue numbers of interest to report after mapping. Each token may "
            "be a single residue or a range (e.g. 118-132). Repeat for multiple "
            "components. Example: --landmark ATP8B:237,260,118-132"
        )
    )
    parser.add_argument(
        "--match-score",
        type=int,
        default=2,
        help="Alignment match score. Default: 2."
    )
    parser.add_argument(
        "--mismatch-score",
        type=int,
        default=-1,
        help="Alignment mismatch penalty. Default: -1."
    )
    parser.add_argument(
        "--gap-score",
        type=int,
        default=-4,
        help="Alignment gap penalty. Default: -4."
    )
    parser.add_argument(
        "--min-mapped-fraction",
        type=float,
        default=0.95,
        help=(
            "Minimum fraction of PDB residues in a component that must be aligned to "
            "a GRO residue (not a gap). Mapping is rejected below this threshold. Default: 0.95."
        )
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.90,
        help=(
            "Minimum sequence identity over the aligned (non-gap) residue pairs. "
            "Mapping is rejected below this threshold. Default: 0.90."
        )
    )
    parser.add_argument(
        "--max-problems",
        type=int,
        default=30,
        help="Maximum number of mismatches or gaps to print per component. Default: 30."
    )
    parser.add_argument(
        "-l", "--log",
        type=str,
        default="map_residues.log",
        help="Log file to record all output. Default: map_residues.log."
    )
    return parser.parse_args(args)


_TSV_FIELDNAMES = [
    "component", "pdb_segid", "pdb_resid", "pdb_icode", "pdb_resname", "pdb_aa",
    "pdb_component_ordinal", "gro_resid", "gro_segment", "gro_resname", "gro_aa",
    "gro_protein_ordinal", "gro_file_ordinal", "aa_match",
]


def _write_tsv(rows, output_tsv):
    with open(output_tsv, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=_TSV_FIELDNAMES, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    t1 = time.time()
    args = initialize(sys.argv[1:])
    sys.stdout = utils.Logger(args.log)
    sys.stderr = utils.Logger(args.log)

    print(f"\nCommand line: {' '.join(sys.argv)}")
    print(f"Current working directory: {os.getcwd()}\n")

    component_ids: dict = {}
    for name, ids in args.components:
        if name in component_ids:
            raise ValueError(
                f"Duplicate component name {name!r}. Each component must appear only once."
            )
        component_ids[name] = ids

    landmark_pdb_residues: dict = {}
    for name, resids in args.landmarks:
        landmark_pdb_residues.setdefault(name, []).extend(resids)

    print(f"Reference PDB:  {args.pdb}")
    print(f"Simulation GRO: {args.gro}")
    print(f"Output TSV:     {args.output}")
    print(f"PDB id field:   {args.pdb_id_field}")
    print(f"\nComponents ({len(component_ids)}):")
    for name, ids in component_ids.items():
        print(f"  {name}: {', '.join(ids)}")

    pdb_components = protein.read_pdb_components(
        args.pdb, component_ids, id_field=args.pdb_id_field
    )
    gro_residues = protein.read_gro_protein_residues(args.gro)

    print("\nPDB component sizes:")
    for component, res_list in pdb_components.items():
        ids = ", ".join(component_ids[component])
        print(f"  {component} ({ids}): {len(res_list)} protein residues")
        if not res_list:
            raise RuntimeError(
                f"No PDB residues found for component {component!r}. "
                f"Check the --pdb-id-field and the IDs provided with --component."
            )

    print(f"\nGRO protein-like residues detected: {len(gro_residues)}")
    if not gro_residues:
        raise RuntimeError(
            "No protein-like residues found in GRO. "
            "Check that the GRO file contains standard amino acid residue names."
        )

    gro_sequence = "".join(r.aa for r in gro_residues)

    all_rows = []
    for component, pdb_residues in pdb_components.items():
        pdb_sequence = "".join(r.aa for r in pdb_residues)
        aln = protein.fitting_needleman_wunsch(
            pdb_sequence,
            gro_sequence,
            match_score=args.match_score,
            mismatch_score=args.mismatch_score,
            gap_score=args.gap_score,
        )
        protein.summarize_alignment(
            component,
            pdb_residues,
            gro_residues,
            aln,
            min_fraction_mapped=args.min_mapped_fraction,
            min_identity=args.min_identity,
            max_problems=args.max_problems,
        )
        all_rows.extend(
            protein.build_mapping_rows(component, pdb_residues, gro_residues, aln)
        )

    protein.check_nonoverlapping_mappings(all_rows, max_problems=args.max_problems)
    protein.check_monotonicity(all_rows, list(component_ids))

    if landmark_pdb_residues:
        protein.print_landmarks(all_rows, landmark_pdb_residues)

    _write_tsv(all_rows, args.output)
    print(f"\nWrote {len(all_rows)} mapped residue pairs to {args.output}")
    print(f"\nElapsed time: {utils.format_time(time.time() - t1)}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"\nFATAL: {exc}", file=sys.stderr)
        sys.exit(1)
