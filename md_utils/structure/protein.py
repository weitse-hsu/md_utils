from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Amino acid lookup tables
# ---------------------------------------------------------------------------

AA3_TO_1: Dict[str, str] = {
    # Standard amino acids
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Common protonation / force-field variants
    "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H", "HIE": "H", "HIP": "H",
    "ASH": "D", "GLH": "E", "LYN": "K", "CYX": "C", "CYM": "C",
    # N-terminal variants
    "NALA": "A", "NARG": "R", "NASN": "N", "NASP": "D", "NCYS": "C",
    "NGLN": "Q", "NGLU": "E", "NGLY": "G", "NHIS": "H", "NILE": "I",
    "NLEU": "L", "NLYS": "K", "NMET": "M", "NPHE": "F", "NPRO": "P",
    "NSER": "S", "NTHR": "T", "NTRP": "W", "NTYR": "Y", "NVAL": "V",
    # C-terminal variants
    "CALA": "A", "CARG": "R", "CASN": "N", "CASP": "D", "CCYS": "C",
    "CGLN": "Q", "CGLU": "E", "CGLY": "G", "CHIS": "H", "CILE": "I",
    "CLEU": "L", "CLYS": "K", "CMET": "M", "CPHE": "F", "CPRO": "P",
    "CSER": "S", "CTHR": "T", "CTRP": "W", "CTYR": "Y", "CVAL": "V",
    # Modified residues commonly encountered in MD work
    "D2P": "D", "SEP": "S", "TPO": "T", "PTR": "Y",
}

_AA1_TO_3: Dict[str, str] = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}


def convert_res_code(res_code):
    """
    Converts an amino acid code between three-letter and one-letter formats.
    Non-standard amino acids or invalid codes are denoted by 'X'.

    Parameters
    ----------
    res_code : str
        The input amino acid code (either three-letter like 'ALA' or one-letter like 'A').

    Returns
    -------
    converted_code : str
        The converted amino acid code (one-letter if input is three-letter, or three-letter if input is one-letter).
        Returns 'X' for invalid or non-standard codes.
    """
    res_code = res_code.upper()
    if len(res_code) == 3:
        return AA3_TO_1.get(res_code, 'X')
    elif len(res_code) == 1:
        return _AA1_TO_3.get(res_code, 'X')
    else:
        raise ValueError(
            f"Invalid amino acid code: {res_code}. It must be either a three-letter or one-letter code."
        )


# ---------------------------------------------------------------------------
# Data classes for residue mapping
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class PDBResidue:
    component: str
    pdb_id: str          # segid or chain ID, depending on id_field used
    pdb_resid: int
    icode: str
    resname: str
    aa: str
    pdb_component_ordinal: int


@dataclass(frozen=True)
class GROResidue:
    gro_resid: int
    resname: str
    aa: str
    gro_protein_ordinal: int
    gro_file_ordinal: int
    gro_segment: int = 0


@dataclass(frozen=True)
class AlignmentResult:
    score: int
    pairs: List[Tuple[Optional[int], Optional[int]]]
    gro_start_index: int
    gro_end_index: int


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _normalize_resname(resname: str) -> str:
    return resname.strip().upper()


def _resname_to_aa(resname: str) -> Optional[str]:
    return AA3_TO_1.get(_normalize_resname(resname))


def _parse_int_field(text: str) -> Optional[int]:
    text = text.strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


# ---------------------------------------------------------------------------
# PDB / GRO parsers
# ---------------------------------------------------------------------------

def read_pdb_components(
    pdb_file: str,
    component_ids: Dict[str, List[str]],
    id_field: str = "segid",
) -> Dict[str, List[PDBResidue]]:
    """
    Parse a PDB file and return protein residues grouped by component.

    Parameters
    ----------
    pdb_file:
        Path to the reference PDB file.
    component_ids:
        Mapping from component name to a list of identifier strings (segment IDs or
        chain IDs) belonging to that component.
    id_field:
        Which PDB field identifies components: ``'segid'`` (columns 73-76) or
        ``'chain'`` (column 22).

    Returns
    -------
    Dict mapping component name to an ordered list of PDBResidue objects (one per
    unique residue, in the order they appear in the PDB).
    """
    if id_field not in ("segid", "chain"):
        raise ValueError(f"id_field must be 'segid' or 'chain', got {id_field!r}.")

    id_to_component: Dict[str, str] = {}
    for component, ids in component_ids.items():
        for id_ in ids:
            if id_ in id_to_component:
                raise ValueError(f"ID {id_!r} is assigned to multiple components.")
            id_to_component[id_] = component

    components_raw: Dict[str, List[Tuple[str, int, str, str, str]]] = {
        component: [] for component in component_ids
    }
    seen: set = set()
    all_ids_seen: set = set()

    with open(pdb_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            resname = _normalize_resname(line[17:20])
            aa = _resname_to_aa(resname)
            if aa is None:
                continue

            pdb_resid = _parse_int_field(line[22:26])
            if pdb_resid is None:
                continue
            icode = line[26].strip() if len(line) > 26 else ""

            if id_field == "segid":
                pdb_id = line[72:76].strip() if len(line) > 72 else ""
            else:
                pdb_id = line[21].strip() if len(line) > 21 else ""

            all_ids_seen.add(pdb_id)
            if pdb_id not in id_to_component:
                continue

            component = id_to_component[pdb_id]
            key = (component, pdb_id, pdb_resid, icode, resname)
            if key in seen:
                continue
            seen.add(key)
            components_raw[component].append((pdb_id, pdb_resid, icode, resname, aa))

    components: Dict[str, List[PDBResidue]] = {}
    for component, rows in components_raw.items():
        components[component] = [
            PDBResidue(
                component=component,
                pdb_id=pdb_id,
                pdb_resid=pdb_resid,
                icode=icode,
                resname=resname,
                aa=aa,
                pdb_component_ordinal=i + 1,
            )
            for i, (pdb_id, pdb_resid, icode, resname, aa) in enumerate(rows)
        ]

    expected_ids = set(id_to_component)
    missing = expected_ids - all_ids_seen
    if missing:
        print(f"WARNING: Expected PDB {id_field}(s) not seen in ATOM/HETATM records: {sorted(missing)}")

    extra = all_ids_seen - expected_ids - {""}
    if extra:
        print(f"NOTE: PDB has additional {id_field}(s) not used for mapping: {sorted(extra)}")

    return components


def read_gro_protein_residues(gro_file: str) -> List[GROResidue]:
    """
    Parse a GROMACS GRO file and return all protein-like residues in order.

    Protein-like residues are those whose residue name maps to a standard amino acid
    via :data:`AA3_TO_1`. One entry is returned per unique (resid, resname) pair.
    """
    residues: List[GROResidue] = []
    last_key = None
    last_unique_resid: Optional[int] = None
    gro_file_ordinal = 0
    gro_segment = 0

    with open(gro_file, "r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    if len(lines) < 3:
        raise ValueError(f"{gro_file!r} does not look like a valid GRO file.")

    natoms = _parse_int_field(lines[1])
    if natoms is None:
        raise ValueError(f"Could not read atom count from second line of {gro_file!r}.")

    atom_lines = lines[2:2 + natoms]
    if len(atom_lines) != natoms:
        raise ValueError(f"Expected {natoms} atom lines in {gro_file!r}, found {len(atom_lines)}.")

    for line in atom_lines:
        gro_resid = _parse_int_field(line[0:5])
        if gro_resid is None:
            continue
        resname = _normalize_resname(line[5:10])
        aa = _resname_to_aa(resname)

        key = (gro_resid, resname)
        if key != last_key:
            gro_file_ordinal += 1
            if last_unique_resid is not None and gro_resid < last_unique_resid:
                gro_segment += 1
            last_unique_resid = gro_resid
            last_key = key
            if aa is not None:
                residues.append(
                    GROResidue(
                        gro_resid=gro_resid,
                        resname=resname,
                        aa=aa,
                        gro_protein_ordinal=len(residues) + 1,
                        gro_file_ordinal=gro_file_ordinal,
                        gro_segment=gro_segment,
                    )
                )

    return residues


# ---------------------------------------------------------------------------
# Sequence alignment
# ---------------------------------------------------------------------------

def fitting_needleman_wunsch(
    query: str,
    target: str,
    match_score: int = 2,
    mismatch_score: int = -1,
    gap_score: int = -4,
) -> AlignmentResult:
    """
    Align the full *query* sequence to the best-scoring subsequence of *target*.

    Free gaps are allowed at the start and end of *target* so that each component
    sequence can be located anywhere within the full GRO protein sequence. The entire
    *query* must be accounted for (no free leading/trailing gaps in the query).
    """
    m = len(query)
    n = len(target)
    if m == 0:
        raise ValueError("Empty query sequence.")
    if n == 0:
        raise ValueError("Empty target sequence.")

    dp = [[0] * (n + 1) for _ in range(m + 1)]
    tb = [[""] * (n + 1) for _ in range(m + 1)]

    for j in range(n + 1):
        dp[0][j] = 0
        tb[0][j] = "L"

    for i in range(1, m + 1):
        dp[i][0] = i * gap_score
        tb[i][0] = "U"

    for i in range(1, m + 1):
        qi = query[i - 1]
        for j in range(1, n + 1):
            tj = target[j - 1]
            diag = dp[i - 1][j - 1] + (match_score if qi == tj else mismatch_score)
            up = dp[i - 1][j] + gap_score
            left = dp[i][j - 1] + gap_score
            best = max(diag, up, left)
            dp[i][j] = best
            if best == diag:
                tb[i][j] = "D"
            elif best == up:
                tb[i][j] = "U"
            else:
                tb[i][j] = "L"

    end_j = max(range(n + 1), key=lambda j: dp[m][j])
    score = dp[m][end_j]

    pairs_rev: List[Tuple[Optional[int], Optional[int]]] = []
    i, j = m, end_j
    while i > 0:
        move = tb[i][j]
        if move == "D":
            pairs_rev.append((i - 1, j - 1))
            i -= 1
            j -= 1
        elif move == "U":
            pairs_rev.append((i - 1, None))
            i -= 1
        elif move == "L":
            pairs_rev.append((None, j - 1))
            j -= 1
        else:
            raise RuntimeError(f"Traceback failed at i={i}, j={j}.")

    pairs = list(reversed(pairs_rev))
    mapped_target_indices = [tj for _, tj in pairs if tj is not None]
    gro_start_index = min(mapped_target_indices) if mapped_target_indices else -1
    gro_end_index = max(mapped_target_indices) if mapped_target_indices else -1

    return AlignmentResult(
        score=score,
        pairs=pairs,
        gro_start_index=gro_start_index,
        gro_end_index=gro_end_index,
    )


# ---------------------------------------------------------------------------
# Mapping helpers
# ---------------------------------------------------------------------------

def summarize_alignment(
    component: str,
    pdb_residues: List[PDBResidue],
    gro_residues: List[GROResidue],
    aln: AlignmentResult,
    min_fraction_mapped: float = 0.95,
    min_identity: float = 0.90,
    max_problems: int = 30,
) -> None:
    """
    Print alignment statistics and raise RuntimeError if quality thresholds are not met.
    """
    mapped = []
    matches = 0
    mismatches = []
    pdb_gaps = []
    gro_gaps = []

    for pdb_i, gro_i in aln.pairs:
        if pdb_i is not None and gro_i is not None:
            pdb_r = pdb_residues[pdb_i]
            gro_r = gro_residues[gro_i]
            mapped.append((pdb_i, gro_i))
            if pdb_r.aa == gro_r.aa:
                matches += 1
            else:
                mismatches.append((pdb_r, gro_r))
        elif pdb_i is not None:
            gro_gaps.append(pdb_residues[pdb_i])
        elif gro_i is not None:
            pdb_gaps.append(gro_residues[gro_i])

    fraction_mapped = len(mapped) / len(pdb_residues) if pdb_residues else 0.0
    identity = matches / len(mapped) if mapped else 0.0

    print(f"\n[{component}] alignment summary")
    print(f"  PDB residues:                  {len(pdb_residues)}")
    print(f"  Mapped PDB residues:           {len(mapped)} ({fraction_mapped:.3f})")
    print(f"  Sequence identity over mapped: {identity:.3f}")
    print(f"  Alignment score:               {aln.score}")
    if aln.gro_start_index >= 0:
        start = gro_residues[aln.gro_start_index]
        end = gro_residues[aln.gro_end_index]
        if start.gro_segment > 0 or end.gro_segment != start.gro_segment:
            if start.gro_segment == end.gro_segment:
                seg_info = f", segment {start.gro_segment}"
            else:
                seg_info = f", segments {start.gro_segment}-{end.gro_segment}"
        else:
            seg_info = ""
        print(
            f"  GRO protein-ordinal span:      "
            f"{start.gro_protein_ordinal}-{end.gro_protein_ordinal} "
            f"(GRO resid {start.gro_resid}-{end.gro_resid}{seg_info})"
        )

    if mismatches:
        print(f"  WARNING: {len(mismatches)} mapped residue AA mismatches. First examples:")
        for pdb_r, gro_r in mismatches[:max_problems]:
            print(
                f"    PDB {component} {pdb_r.pdb_id}:{pdb_r.pdb_resid}{pdb_r.icode} "
                f"{pdb_r.resname}/{pdb_r.aa}  <->  "
                f"GRO resid {gro_r.gro_resid} {gro_r.resname}/{gro_r.aa}"
            )

    if gro_gaps:
        print(f"  WARNING: {len(gro_gaps)} PDB residues aligned to gaps in GRO. First examples:")
        for pdb_r in gro_gaps[:max_problems]:
            print(
                f"    PDB {component} {pdb_r.pdb_id}:{pdb_r.pdb_resid}{pdb_r.icode} "
                f"{pdb_r.resname}/{pdb_r.aa}"
            )

    if pdb_gaps:
        print(f"  NOTE: {len(pdb_gaps)} GRO residues aligned to gaps in this PDB component. First examples:")
        for gro_r in pdb_gaps[:max_problems]:
            print(f"    GRO resid {gro_r.gro_resid} {gro_r.resname}/{gro_r.aa}")

    if fraction_mapped < min_fraction_mapped:
        raise RuntimeError(
            f"{component}: only {fraction_mapped:.3f} of PDB residues were mapped "
            f"(threshold: {min_fraction_mapped}). Mapping is likely wrong."
        )
    if identity < min_identity:
        raise RuntimeError(
            f"{component}: sequence identity {identity:.3f} is below threshold "
            f"{min_identity}. Mapping is likely wrong."
        )


def build_mapping_rows(
    component: str,
    pdb_residues: List[PDBResidue],
    gro_residues: List[GROResidue],
    aln: AlignmentResult,
) -> List[Dict[str, object]]:
    """Return one row dict per aligned (PDB residue, GRO residue) pair."""
    rows = []
    for pdb_i, gro_i in aln.pairs:
        if pdb_i is None or gro_i is None:
            continue
        pdb_r = pdb_residues[pdb_i]
        gro_r = gro_residues[gro_i]
        rows.append({
            "component": component,
            "pdb_segid": pdb_r.pdb_id,
            "pdb_resid": pdb_r.pdb_resid,
            "pdb_icode": pdb_r.icode,
            "pdb_resname": pdb_r.resname,
            "pdb_aa": pdb_r.aa,
            "pdb_component_ordinal": pdb_r.pdb_component_ordinal,
            "gro_resid": gro_r.gro_resid,
            "gro_segment": gro_r.gro_segment,
            "gro_resname": gro_r.resname,
            "gro_aa": gro_r.aa,
            "gro_protein_ordinal": gro_r.gro_protein_ordinal,
            "gro_file_ordinal": gro_r.gro_file_ordinal,
            "aa_match": pdb_r.aa == gro_r.aa,
        })
    return rows


def check_nonoverlapping_mappings(
    all_rows: List[Dict[str, object]],
    max_problems: int = 30,
) -> None:
    """Raise RuntimeError if any GRO residue is mapped from more than one component."""
    gro_to_rows: Dict[int, List[Dict[str, object]]] = {}
    for row in all_rows:
        ordinal = int(row["gro_protein_ordinal"])
        gro_to_rows.setdefault(ordinal, []).append(row)

    overlaps = {k: v for k, v in gro_to_rows.items() if len(v) > 1}
    if overlaps:
        print("\nERROR: Multiple PDB residues/components mapped to the same GRO protein residue.")
        for _, rows in list(overlaps.items())[:max_problems]:
            print("  Overlap:")
            for row in rows:
                print(
                    f"    {row['component']} {row['pdb_segid']}:{row['pdb_resid']} "
                    f"-> GRO resid {row['gro_resid']} ({row['gro_resname']})"
                )
        raise RuntimeError("Non-unique GRO mapping detected.")


def check_monotonicity(
    all_rows: List[Dict[str, object]],
    component_names: List[str],
) -> None:
    """Raise RuntimeError if GRO ordinal order is not strictly increasing within any component."""
    print("\nMonotonicity checks")
    for component in component_names:
        component_rows = [row for row in all_rows if row["component"] == component]
        gro_ordinals = [int(row["gro_protein_ordinal"]) for row in component_rows]
        is_strictly_increasing = all(b > a for a, b in zip(gro_ordinals, gro_ordinals[1:]))
        print(f"  {component}: GRO order strictly increasing = {is_strictly_increasing}")
        if not is_strictly_increasing:
            raise RuntimeError(f"{component}: mapping is not monotonic in GRO order.")


def print_landmarks(
    all_rows: List[Dict[str, object]],
    landmark_pdb_residues: Dict[str, List],
) -> None:
    """Print the GRO mapping for each requested landmark PDB residue or range."""
    print("\nLandmark residue mappings")
    for component, tokens in landmark_pdb_residues.items():
        if not tokens:
            continue
        print(f"  {component}:")
        for token in tokens:
            if isinstance(token, tuple):
                start, end = token
                hits = [
                    row for row in all_rows
                    if row["component"] == component
                    and start <= int(row["pdb_resid"]) <= end
                ]
                if hits:
                    gro_resids = [int(row["gro_resid"]) for row in hits]
                    n_total = end - start + 1
                    n_mapped = len(hits)
                    note = "" if n_mapped == n_total else f" ({n_mapped}/{n_total} PDB residues mapped)"
                    print(
                        f"    PDB {start}-{end} "
                        f"-> GRO resid {min(gro_resids)}-{max(gro_resids)}"
                        f"{note}"
                    )
                else:
                    print(f"    PDB resid {start}-{end}: NOT MAPPED")
            else:
                pdb_resid = token
                hits = [
                    row for row in all_rows
                    if row["component"] == component and int(row["pdb_resid"]) == pdb_resid
                ]
                if hits:
                    for row in hits:
                        print(
                            f"    PDB {row['pdb_segid']}:{row['pdb_resid']} {row['pdb_resname']} "
                            f"-> GRO resid {row['gro_resid']} {row['gro_resname']} "
                            f"(GRO protein ordinal {row['gro_protein_ordinal']})"
                        )
                else:
                    print(f"    PDB resid {pdb_resid}: NOT MAPPED")
