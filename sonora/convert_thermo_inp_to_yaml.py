#!/usr/bin/env python3
"""Convert a CEA/NASA thermo .inp file into the YAML format used by equilibrate."""

from __future__ import annotations

import argparse
import math
import pathlib
import re
from collections import OrderedDict
from typing import Dict, List, Sequence, Tuple
from photochem.utils._format import FormatReactions_main, yaml, MyDumper

FLOAT_RE = re.compile(r"[+-]?(?:\d+\.\d*|\d*\.\d+|\d+)(?:[DE][+-]?\d+)?")
ELECTRON_MASS = 5.48579909e-4


def clean_float(text: str) -> float:
    """Parse a float while handling Fortran's D exponent."""
    return float(text.replace("D", "E").replace("d", "E"))


def extract_floats(text: str, expected: int | None = None) -> List[float]:
    values = [clean_float(x) for x in FLOAT_RE.findall(text)]
    if expected is not None and len(values) != expected:
        raise ValueError(f"Expected {expected} floats, found {len(values)} in '{text}'")
    return values


def normalise_atom(symbol: str) -> str:
    """Normalise atom symbols (AL -> Al, CL -> Cl)."""
    if not symbol:
        return symbol
    return symbol[0].upper() + symbol[1:].lower()


def parse_header(line: str) -> Tuple[int, Dict[str, float], bool, float | None]:
    """Parse the header line that encodes composition and bookkeeping fields."""
    n_ranges = int(line[:2].strip())
    composition_fields = [line[10 + i * 8 : 10 + (i + 1) * 8] for i in range(5)]

    composition: OrderedDict[str, float] = OrderedDict()
    for field in composition_fields:
        symbol = normalise_atom(field[:2].strip())
        count_str = field[2:].strip()
        if not symbol or not count_str:
            continue
        count_val = clean_float(count_str)
        rounded = round(count_val)
        if math.isclose(count_val, rounded):
            count_val = int(rounded)
        composition[symbol] = count_val

    rest_tokens = line[50:].split()
    condensate = bool(rest_tokens and rest_tokens[0] != "0")
    mol_weight = clean_float(rest_tokens[1]) if len(rest_tokens) > 1 else None

    return n_ranges, composition, condensate, mol_weight


def clean_temperature(value: float) -> float:
    """Snap temperatures suffering from concatenated digit artifacts (e.g., 300.0007 -> 300.0)."""
    rounded_int = round(value)
    if abs(value - rounded_int) < 5e-3:
        return float(rounded_int)
    return value


def parse_species(lines: List[str], start: int) -> Tuple[Dict, int]:
    """Parse a single species block starting at `start` (name line)."""
    name_line = lines[start].rstrip("\n")
    if not name_line.strip():
        raise ValueError(f"Empty species name line at {start}")
    name = name_line[:18].strip()

    header = lines[start + 1].rstrip("\n")
    n_ranges, composition, condensate, mol_weight = parse_header(header)
    if name.strip().lower() in {"e", "e-", "electron"}:
        condensate = False

    temp_ranges: List[float] = []
    thermo_data: List[List[float]] = []

    idx = start + 2
    for i in range(n_ranges):
        temp_line = lines[idx].rstrip("\n")
        temps = extract_floats(temp_line)
        if len(temps) < 2:
            raise ValueError(f"Could not find two temperatures for {name} on line {idx+1}")
        low, high = clean_temperature(temps[0]), clean_temperature(temps[1])
        if i == 0:
            temp_ranges.extend([low, high])
        else:
            temp_ranges.append(high)

        coeffs_first = extract_floats(lines[idx + 1])
        coeffs_second_raw = extract_floats(lines[idx + 2])

        if len(coeffs_first) != 5:
            raise ValueError(
                f"Unexpected thermo coefficient count for {name} (range {i+1}) on line 1: "
                f"{len(coeffs_first)}"
            )

        if len(coeffs_second_raw) == 5:
            coeffs_second = coeffs_second_raw[:2] + coeffs_second_raw[3:]
        elif len(coeffs_second_raw) == 4:
            coeffs_second = coeffs_second_raw
        else:
            raise ValueError(
                f"Unexpected thermo coefficient count for {name} (range {i+1}) on line 2: "
                f"{len(coeffs_second_raw)}"
            )

        coeffs = coeffs_first + coeffs_second
        thermo_data.append(coeffs)
        idx += 3

    species = {
        "name": name,
        "composition": composition,
        "condensate": condensate,
        "thermo": {
            "model": "NASA9",
            "temperature-ranges": temp_ranges,
            "data": thermo_data,
        },
        "mass": mol_weight,
    }
    return species, idx


def standardize_species_name(name: str, composition: Dict[str, float]) -> str:
    """Build a canonical species name from its composition (e.g., HCL -> HCl).

    Preserves trailing phase/annotation suffixes like "(c)" or "(L)".
    """
    if not composition:
        return name

    suffix = ""
    m = re.match(r"^(.*?)(\([^()]*\))$", name.strip())
    base_name = name.strip()
    if m:
        base_name, suffix = m.group(1), m.group(2)

    charge = 0
    comp_no_e = OrderedDict()
    for atom, count in composition.items():
        if atom == "E":
            charge -= count
        else:
            comp_no_e[atom] = count

    if not comp_no_e:
        return name

    parts: List[str] = []
    for atom, count in comp_no_e.items():
        atom_norm = normalise_atom(atom)
        cnt_round = round(count)
        cnt = int(cnt_round) if math.isclose(count, cnt_round) else count
        if cnt == 1:
            parts.append(atom_norm)
        else:
            parts.append(f"{atom_norm}{cnt}")

    charge_suffix = ""
    if charge:
        if isinstance(charge, float) and not charge.is_integer():
            charge_suffix = f"{charge:+}"
        else:
            q = int(charge)
            if q > 0:
                charge_suffix = "+" * q
            else:
                charge_suffix = "-" * abs(q)

    return "".join(parts) + charge_suffix + suffix


def maybe_normalize_species_name(
    name: str,
    composition: Dict[str, float],
    existing_names: set[str],
    original_names: set[str],
) -> str:
    """
    Optionally normalize a species name based on casing/atoms while avoiding duplicates.

    Keeps names unchanged unless:
    - the name is all uppercase,
    - the species uses at least one multi-letter atom (e.g., Al, Ca),
    - the standardized name differs, and
    - the standardized name does not collide with existing or other original names.
    """
    candidate = standardize_species_name(name, composition)
    same_letters_ignore_case = candidate.upper() == name.upper()
    if (
        candidate != name
        and same_letters_ignore_case
        and candidate not in existing_names
        and (candidate == name or candidate not in original_names)
    ):
        return candidate
    return name


def format_number(value) -> str:
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float) and value == 0.0:
        return "0.0"
    return f"{value:.12g}"


def to_dict(atoms: List[Tuple[str, float]], species: List[Dict], output: pathlib.Path) -> None:
    atoms_yaml = [{"name": name, "mass": mass} for name, mass in atoms]

    species_yaml = []
    for sp in species:
        entry = {
            "name": sp["name"],
            "composition": dict(sp["composition"]),
            "thermo": sp["thermo"],
        }
        if sp["condensate"]:
            entry["condensate"] = True
        species_yaml.append(entry)

    content = {"atoms": atoms_yaml, "species": species_yaml}
    return content


def nasa9_cp_over_r(coeffs: Sequence[float], temperature: float) -> float:
    """Dimensionless Cp/R from NASA9 coefficients at a temperature."""
    a1, a2, a3, a4, a5, a6, a7, _, _ = coeffs
    t = temperature
    return a1 / (t * t) + a2 / t + a3 + a4 * t + a5 * t * t + a6 * t**3 + a7 * t**4


def sample_temperatures(t_low: float, t_high: float, n_samples: int = 32) -> List[float]:
    """Log-space samples including interval endpoints."""
    if t_high <= t_low:
        return [t_low]
    if t_low <= 0.0:
        step = (t_high - t_low) / (n_samples - 1)
        return [t_low + i * step for i in range(n_samples)]
    log_low = math.log(t_low)
    log_high = math.log(t_high)
    return [math.exp(log_low + i * (log_high - log_low) / (n_samples - 1)) for i in range(n_samples)]


def species_has_unphysical_cp(
    sp: Dict, cp_fix_tmin: float | None = None, cp_fix_tmax: float | None = None
) -> bool:
    """Return True if a species has non-positive or non-finite Cp in the requested T window."""
    temp_ranges = sp["thermo"]["temperature-ranges"]
    blocks = sp["thermo"]["data"]

    t_window_low = cp_fix_tmin if cp_fix_tmin is not None else temp_ranges[0]
    t_window_high = cp_fix_tmax if cp_fix_tmax is not None else temp_ranges[-1]
    if t_window_high <= t_window_low:
        raise ValueError("cp_fix_tmax must be greater than cp_fix_tmin")

    for i, coeffs in enumerate(blocks):
        t_low = max(temp_ranges[i], t_window_low)
        t_high = min(temp_ranges[i + 1], t_window_high)
        if t_high <= t_low:
            continue
        for t in sample_temperatures(t_low, t_high):
            cp_r = nasa9_cp_over_r(coeffs, t)
            if not math.isfinite(cp_r) or cp_r <= 0.0:
                return True
    return False


def dof_cp_over_r(composition: Dict[str, float]) -> float:
    """
    Simple ideal-gas equipartition fallback Cp/R.

    Monatomic: 2.5, Diatomic: 3.5, Polyatomic: 4.0 (rotational DoF active).
    """
    n_atoms = 0.0
    for atom, count in composition.items():
        if atom == "E":
            continue
        n_atoms += float(count)
    n_atoms = int(round(n_atoms))

    if n_atoms <= 1:
        return 2.5
    if n_atoms == 2:
        return 3.5
    return 4.0


def replace_species_cp_with_dof(sp: Dict) -> None:
    """Replace all NASA9 blocks with constant-Cp fallback based on composition DoF."""
    cp_r = dof_cp_over_r(sp["composition"])
    n_blocks = len(sp["thermo"]["data"])
    sp["thermo"]["data"] = [[0.0, 0.0, cp_r, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] for _ in range(n_blocks)]


def convert_file(
    input_path: pathlib.Path,
    output_path: pathlib.Path,
    species: Sequence[str] | None = None,
    normalize_species_case: bool = False,
    dedupe_names: bool = False,
    cp_fix_mode: str | None = None,
    cp_fix_tmin: float | None = None,
    cp_fix_tmax: float | None = None,
    thermo_tmin_override: float | None = None,
    thermo_tmax_override: float | None = None,
) -> None:
    raw_lines = input_path.read_text().splitlines()
    species_list: List[Dict] = []
    atom_masses: OrderedDict[str, float] = OrderedDict()
    atoms_seen: OrderedDict[str, None] = OrderedDict()

    # Fallback atomic masses for atoms that may not have isolated entries in the .inp file
    periodic_mass = {
        "H": 1.00794,
        "He": 4.002602,
        "C": 12.0107,
        "N": 14.0067,
        "O": 15.9994,
        "Na": 22.98977,
        "Mg": 24.305,
        "Al": 26.981538,
        "Si": 28.0855,
        "P": 30.973761,
        "S": 32.065,
        "Cl": 35.453,
        "K": 39.0983,
        "Ca": 40.078,
        "Ti": 47.867,
        "V": 50.9415,
        "Fe": 55.845,
        "Ni": 58.6934,
    }

    # First pass: gather original species names to avoid renaming collisions
    def collect_original_names(lines: List[str]) -> set[str]:
        names = set()
        i = 0
        n = len(lines)
        while i < n - 1:
            hdr = lines[i + 1] if i + 1 < n else ""
            if not hdr[:2].strip().isdigit():
                i += 1
                continue
            name_line = lines[i].rstrip("\n")
            names.add(name_line[:18].strip())
            n_ranges, *_ = parse_header(hdr.rstrip("\n"))
            i += 2 + n_ranges * 3
        return names

    original_names = collect_original_names(raw_lines)

    requested_species = set(species) if species is not None else None
    matched_species: set[str] = set()

    idx = 0
    n_lines = len(raw_lines)
    existing_names = set()

    while idx < n_lines:
        line = raw_lines[idx]
        if not line.strip():
            idx += 1
            continue
        # Detect start of a species: name line followed by a header that starts with a number
        header_candidate = raw_lines[idx + 1] if idx + 1 < n_lines else ""
        if not header_candidate[:2].strip().isdigit():
            idx += 1
            continue

        sp, idx = parse_species(raw_lines, idx)
        if normalize_species_case:
            sp["name"] = maybe_normalize_species_name(
                sp["name"], sp["composition"], existing_names, original_names
            )
        if requested_species is not None and sp["name"] not in requested_species:
            continue

        if dedupe_names and sp["name"] in existing_names:
            pass  # skip duplicate name
        else:
            species_list.append(sp)
            existing_names.add(sp["name"])
            if requested_species is not None:
                matched_species.add(sp["name"])

        # Record atomic masses for elemental species (single atom with count 1)
        if len(sp["composition"]) == 1:
            atom, count = next(iter(sp["composition"].items()))
            if count == 1 and atom not in atom_masses and sp["mass"] is not None:
                atom_masses[atom] = sp["mass"]
        # Track atom usage order
        for atom in sp["composition"].keys():
            atoms_seen.setdefault(atom, None)

    if requested_species is not None:
        missing_species = [name for name in species if name not in matched_species]
        if missing_species:
            raise ValueError(
                "Requested species were not found in input file: "
                + ", ".join(missing_species)
            )

    if thermo_tmin_override is not None:
        for sp in species_list:
            temp_ranges = sp["thermo"]["temperature-ranges"]
            if thermo_tmin_override >= temp_ranges[1]:
                raise ValueError(
                    f"thermo_tmin_override={thermo_tmin_override} must be less than the "
                    f"second temperature bound ({temp_ranges[1]}) for species '{sp['name']}'."
                )
            temp_ranges[0] = float(thermo_tmin_override)

    if thermo_tmax_override is not None:
        for sp in species_list:
            temp_ranges = sp["thermo"]["temperature-ranges"]
            if thermo_tmax_override <= temp_ranges[-2]:
                raise ValueError(
                    f"thermo_tmax_override={thermo_tmax_override} must be greater than the "
                    f"second-to-last temperature bound ({temp_ranges[-2]}) for species "
                    f"'{sp['name']}'."
                )
            temp_ranges[-1] = float(thermo_tmax_override)

    if thermo_tmin_override is not None and thermo_tmax_override is not None:
        for sp in species_list:
            temp_ranges = sp["thermo"]["temperature-ranges"]
            if temp_ranges[0] >= temp_ranges[-1]:
                raise ValueError(
                    f"thermo_tmin_override ({temp_ranges[0]}) must be less than "
                    f"thermo_tmax_override ({temp_ranges[-1]}) for species '{sp['name']}'."
                )

    if cp_fix_mode is not None:
        if cp_fix_mode != "dof":
            raise ValueError(f"Unsupported cp_fix_mode '{cp_fix_mode}'. Supported: dof")
        for sp in species_list:
            if species_has_unphysical_cp(sp, cp_fix_tmin=cp_fix_tmin, cp_fix_tmax=cp_fix_tmax):
                replace_species_cp_with_dof(sp)

    if "E" not in atom_masses:
        atom_masses["E"] = ELECTRON_MASS
    atoms_seen.setdefault("E", None)

    atoms: List[Tuple[str, float]] = []
    for atom in atoms_seen.keys():
        if atom in atom_masses:
            atoms.append((atom, atom_masses[atom]))
        elif atom in periodic_mass:
            atoms.append((atom, periodic_mass[atom]))
        else:
            raise ValueError(f"Missing atomic mass for atom '{atom}'. Add to periodic_mass map.")

    out = to_dict(atoms, species_list, output_path)
    out = FormatReactions_main(out)

    with open(output_path,'w') as f:
        yaml.dump(out, f, Dumper=MyDumper, sort_keys=False, width=80)

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a CEA/NASA thermo .inp file into the equilibrate YAML format."
    )
    parser.add_argument("input", type=pathlib.Path, help="Path to the .inp file to convert")
    parser.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        help="Output YAML path (defaults to input name with .yaml extension)",
    )
    parser.add_argument(
        "--normalize-species-case",
        action="store_true",
        help="Rename species using standardized atom casing (e.g., HCL -> HCl)",
    )
    parser.add_argument(
        "--dedupe-names",
        action="store_true",
        help="Skip duplicate species entries with the same (possibly normalized) name",
    )
    parser.add_argument(
        "--species",
        nargs="+",
        help="Optional list of species names to include; errors if any are missing",
    )
    parser.add_argument(
        "--cp-fix-mode",
        choices=["dof"],
        help="Optional Cp fix mode for unphysical species (climate-only output).",
    )
    parser.add_argument(
        "--cp-fix-tmin",
        type=float,
        help="Lower temperature bound for Cp physicality checks (default: species range minimum).",
    )
    parser.add_argument(
        "--cp-fix-tmax",
        type=float,
        help="Upper temperature bound for Cp physicality checks (default: species range maximum).",
    )
    parser.add_argument(
        "--thermo-tmin-override",
        type=float,
        help="Override the lowest temperature bound for all species thermo polynomials.",
    )
    parser.add_argument(
        "--thermo-tmax-override",
        type=float,
        help="Override the highest temperature bound for all species thermo polynomials.",
    )
    args = parser.parse_args()

    input_path = args.input
    output_path = args.output or input_path.with_suffix(".yaml")

    convert_file(
        input_path,
        output_path,
        species=args.species,
        normalize_species_case=args.normalize_species_case,
        dedupe_names=args.dedupe_names,
        cp_fix_mode=args.cp_fix_mode,
        cp_fix_tmin=args.cp_fix_tmin,
        cp_fix_tmax=args.cp_fix_tmax,
        thermo_tmin_override=args.thermo_tmin_override,
        thermo_tmax_override=args.thermo_tmax_override,
    )
    print(f"Wrote YAML to {output_path}")


if __name__ == "__main__":
    main()
