# Quercetin permeability analysis (aglycone vs circulating metabolites)

In silico physicochemical and ADMET permeability comparison of **quercetin** and
its major **circulating human metabolites**, supporting the interpretation that
the free **aglycone** — not its conjugated plasma forms — is the species capable
of passive cellular and nuclear entry.

## Context

This analysis accompanies the manuscript:

> Juan C.G., Šimunić-Briški N., Volpe-Zanutto F., et al. *Quercetin activates the
> SIRT6–Nrf2 axis during oxidative stress, modulating DNA repair and ageing-
> associated markers in healthy men.* medRxiv (preprint).

It was added during manuscript **revision** to address a reviewer point on
quercetin bioavailability and on why the measured aglycone fraction is the
mechanistically relevant species. It is supporting context for the experimental
(human, mass-spectrometry, 3D imaging) data, not a standalone result.

## What it does

`quercetin_permeability.py` computes physicochemical and passive-permeability
descriptors for quercetin and four circulating metabolites (quercetin-3-O-
glucuronide, quercetin-3-O-sulfate, isorhamnetin, isoquercetin) using RDKit, and
writes a descriptor table (CSV + Markdown). `make_figure.py` produces the
supplementary permeability figure (TPSA vs cLogP, permeation-favourable region
shaded).

Predicted gastrointestinal absorption values were additionally obtained from
SwissADME (Daina et al., 2017) and are reported in the manuscript table.

## Molecules (SMILES from PubChem)

| Compound | PubChem CID |
| --- | --- |
| Quercetin (aglycone) | 5280343 |
| Quercetin-3-O-glucuronide | 12004528 |
| Quercetin-3-O-sulfate | 5280362 |
| Isorhamnetin (3'-O-methyl) | 5281654 |
| Isoquercetin (3-O-glucoside) | 5280804 |

## Run

```bash
pip install rdkit matplotlib
python quercetin_permeability.py     # descriptor table -> results/
python make_figure.py                # supplementary figure -> results/
```

## Key result

The aglycone and isorhamnetin fall within the passive-permeation-favourable range
(TPSA <= 140 A^2, cLogP >= 0; high predicted GI absorption), whereas the
glucuronide, sulfate, and glucoside conjugates that predominate in plasma fall
outside it (low predicted GI absorption; multiple Lipinski/Veber violations).

## Important note on interpretation

These are in silico predictions interpreted as **relative** indicators across the
metabolite series, not absolute permeability measures — tools trained largely on
drug-like chemical space handle ionisable conjugates approximately. The defensible
conclusion is the direction of the aglycone-vs-conjugate contrast. Even the
aglycone shows only moderate predicted permeability in absolute terms; the
analysis explains a relative accessibility difference, consistent with the
measured human data.

Environment: RDKit v2026.03.2; SwissADME accessed 30 May 2026.
