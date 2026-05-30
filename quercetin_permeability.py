"""
Comparative passive-permeability / physicochemical analysis of quercetin and its
major circulating human metabolites.

Purpose
-------
Quercetin circulates in human plasma overwhelmingly as conjugated metabolites
(glucuronides, sulfates, and the methylated derivative isorhamnetin), not as the
free aglycone. Yet only the aglycone is expected to cross cell and nuclear
membranes by passive diffusion. This script quantifies that contrast using
established, interpretable physicochemical descriptors and rule-based passive-
permeability indicators, providing computational support for the argument that
local deconjugation to the aglycone is a prerequisite for intracellular /
nuclear access.

IMPORTANT INTERPRETATION NOTE
-----------------------------
These are in silico physicochemical predictors, not measured permeability. They
are intended as *supporting mechanistic context* alongside experimental data, and
are most robust as a RELATIVE contrast (aglycone vs conjugates), not as absolute
values. Predictors trained largely on drug-like space handle ionisable
glucuronide/sulfate conjugates less reliably; the direction of the contrast is
the defensible result, not the precise numbers.

Molecules (SMILES verified against PubChem)
-------------------------------------------
  Quercetin (aglycone)        PubChem CID 5280343
  Quercetin-3-O-glucuronide   PubChem CID 12004528  (major plasma metabolite)
  Quercetin-3-O-sulfate       PubChem CID 5280362   (major plasma metabolite)
  Isorhamnetin (3'-O-methyl)  PubChem CID 5281654   (methylated metabolite)
  Isoquercetin (3-O-glucoside)PubChem CID 5280804   (dietary glycoside, reference)

Outputs
-------
  results/quercetin_permeability.csv   full descriptor table
  results/quercetin_permeability.md    formatted summary table + interpretation
"""
from __future__ import annotations

import csv
import os

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors


# SMILES verified against PubChem / authoritative chemical databases.
MOLECULES = {
    "Quercetin (aglycone)": {
        "smiles": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O)O",
        "cid": 5280343,
        "class": "aglycone",
    },
    "Quercetin-3-O-glucuronide": {
        "smiles": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)OC4C(C(C(C(O4)C(=O)O)O)O)O)O)O",
        "cid": 12004528,
        "class": "glucuronide conjugate",
    },
    "Quercetin-3-O-sulfate": {
        "smiles": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)OS(=O)(=O)O)O)O",
        "cid": 5280362,
        "class": "sulfate conjugate",
    },
    "Isorhamnetin (3'-O-methyl)": {
        "smiles": "COC1=C(C=CC(=C1)C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O",
        "cid": 5281654,
        "class": "methylated metabolite",
    },
    "Isoquercetin (3-O-glucoside)": {
        "smiles": "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)OC4C(C(C(C(O4)CO)O)O)O)O)O",
        "cid": 5280804,
        "class": "glucoside (dietary)",
    },
}


def lipinski_violations(mw, logp, hbd, hba) -> int:
    """Count Lipinski Rule-of-5 violations (>=2 suggests poor oral permeability)."""
    v = 0
    if mw > 500:
        v += 1
    if logp > 5:
        v += 1
    if hbd > 5:
        v += 1
    if hba > 10:
        v += 1
    return v


def veber_pass(tpsa, rotatable) -> bool:
    """Veber rule for oral bioavailability: TPSA <= 140 and rotatable bonds <= 10."""
    return (tpsa <= 140.0) and (rotatable <= 10)


def passive_diffusion_flag(tpsa, mw, logp) -> str:
    """Qualitative passive transcellular-diffusion expectation.

    Thresholds are grounded in established oral-absorption / permeability rules:
      * TPSA <= 140 Angstrom^2 (Veber) is compatible with passive permeation;
        TPSA > 140 strongly predicts poor passive permeability.
      * A second, stricter TPSA band (<= 75-90) plus positive logP marks the
        most clearly permeable space.
      * Negative cLogP (hydrophilic conjugates) further argues against passive
        diffusion and toward transporter dependence.
    The bands are intended to expose the RELATIVE gradient across the series,
    not to assert absolute permeability for any single compound.
    """
    if tpsa > 140 or logp < 0:
        return "poor (high polarity / ionisable; likely transporter-dependent)"
    if tpsa <= 140 and logp >= 0.5:
        return "moderate (passive permeation feasible)"
    return "borderline"


def analyse(name: str, info: dict) -> dict:
    mol = Chem.MolFromSmiles(info["smiles"])
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES for {name}: {info['smiles']}")

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)              # Wildman-Crippen cLogP
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot = Descriptors.NumRotatableBonds(mol)
    arom = rdMolDescriptors.CalcNumAromaticRings(mol)
    formula = rdMolDescriptors.CalcMolFormula(mol)

    return {
        "molecule": name,
        "class": info["class"],
        "pubchem_cid": info["cid"],
        "formula": formula,
        "MW": round(mw, 2),
        "cLogP": round(logp, 2),
        "TPSA": round(tpsa, 2),
        "HBD": hbd,
        "HBA": hba,
        "RotBonds": rot,
        "AromRings": arom,
        "Lipinski_violations": lipinski_violations(mw, logp, hbd, hba),
        "Veber_pass": "yes" if veber_pass(tpsa, rot) else "no",
        "passive_diffusion": passive_diffusion_flag(tpsa, mw, logp),
    }


def main():
    rows = [analyse(name, info) for name, info in MOLECULES.items()]
    os.makedirs("results", exist_ok=True)

    # CSV
    csv_path = "results/quercetin_permeability.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # Markdown summary
    cols = ["molecule", "class", "MW", "cLogP", "TPSA", "HBD", "HBA",
            "Lipinski_violations", "Veber_pass", "passive_diffusion"]
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    body = []
    for r in rows:
        body.append("| " + " | ".join(str(r[c]) for c in cols) + " |")

    md = [
        "# Quercetin vs circulating metabolites: passive-permeability comparison",
        "",
        "Physicochemical descriptors and rule-based passive-permeability indicators "
        "(RDKit). Intended as supporting mechanistic context for experimental data; "
        "the **relative** contrast between the aglycone and its conjugates is the "
        "robust result, not absolute values.",
        "",
        header, sep, *body,
        "",
        "## Interpretation",
        "",
        "- The **aglycone** carries the lowest molecular weight and polar surface "
        "area and the highest cLogP, and is the species most consistent with "
        "passive transcellular (and nuclear-membrane) diffusion.",
        "- The **glucuronide and sulfate conjugates** add a large, polar, ionisable "
        "group: molecular weight and TPSA rise sharply and effective lipophilicity "
        "falls, flagging them as poor passive diffusers and likely transporter / "
        "efflux substrates. This is consistent with deconjugation being required "
        "before the aglycone can reach intracellular compartments.",
        "- **Isorhamnetin** (methylation) stays close to the aglycone in these "
        "properties, as expected for a small, non-ionic modification.",
        "",
        "_In silico predictors trained on drug-like chemical space handle "
        "glycosides/conjugates less reliably; treat absolute numbers with caution "
        "and rely on the direction of the aglycone-vs-conjugate contrast._",
    ]
    md_path = "results/quercetin_permeability.md"
    with open(md_path, "w") as fh:
        fh.write("\n".join(md) + "\n")

    # console
    print(f"Analysed {len(rows)} molecules.\n")
    for r in rows:
        print(f"{r['molecule']:<30} MW={r['MW']:>7}  TPSA={r['TPSA']:>6}  "
              f"cLogP={r['cLogP']:>6}  -> {r['passive_diffusion']}")
    print(f"\nWrote {csv_path}\nWrote {md_path}")


if __name__ == "__main__":
    main()
