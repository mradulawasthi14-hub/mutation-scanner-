"""
Mutation Scanner: what happens if EVERY letter of a gene mutates?

For each position in a coding sequence, try the 3 other bases and classify the result:
  silent   - same amino acid (the protein doesn't change)
  missense - different amino acid
  nonsense - creates a STOP codon (protein cut short)
  stop-loss - the original STOP codon is destroyed (protein runs on)

Usage:
    python mutation_scanner.py                      # uses data/example.fasta
    python mutation_scanner.py path/to/your_gene.fasta
"""
import csv, os, sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASES = "TCAG"
AMINO = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {a + b + c: AMINO[i] for i, (a, b, c) in
               enumerate((x, y, z) for x in BASES for y in BASES for z in BASES)}


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
            elif line:
                seq.append(line.upper().replace("U", "T"))
    return name or "sequence", "".join(seq)


def classify(old_aa, new_aa):
    if old_aa == new_aa:
        return "silent"
    if new_aa == "*":
        return "nonsense"
    if old_aa == "*":
        return "stop-loss"
    return "missense"


def scan(seq):
    rows = []
    for pos in range(len(seq)):
        codon_start = pos - pos % 3
        codon = seq[codon_start:codon_start + 3]
        if len(codon) < 3:
            break
        old_aa = CODON_TABLE[codon]
        for base in "ACGT":
            if base == seq[pos]:
                continue
            new_codon = codon[:pos % 3] + base + codon[pos % 3 + 1:]
            new_aa = CODON_TABLE[new_codon]
            rows.append({
                "position": pos + 1, "codon_number": codon_start // 3 + 1,
                "ref_base": seq[pos], "alt_base": base,
                "ref_codon": codon, "alt_codon": new_codon,
                "ref_aa": old_aa, "alt_aa": new_aa,
                "effect": classify(old_aa, new_aa),
            })
    return rows


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join("data", "example.fasta")
    name, seq = read_fasta(path)
    if len(seq) % 3:
        print(f"Warning: length {len(seq)} is not a multiple of 3; trailing bases ignored.")
    rows = scan(seq)

    os.makedirs("results", exist_ok=True)
    with open("results/mutations.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=rows[0].keys())
        w.writeheader(); w.writerows(rows)

    counts = Counter(r["effect"] for r in rows)
    total = len(rows)
    print(f"{name}: {len(seq)} bases, {total} possible single-base mutations")
    for k in ("silent", "missense", "nonsense", "stop-loss"):
        print(f"  {k:<10} {counts[k]:>5}  ({100 * counts[k] / total:.1f}%)")

    # figure 1: overall pie/bar, figure 2: per-codon danger score
    colors = {"silent": "#4c9f70", "missense": "#e0a030", "nonsense": "#c0392b", "stop-loss": "#6c5ce7"}
    n_codons = len(seq) // 3
    per_codon = {k: [0] * n_codons for k in colors}
    for r in rows:
        per_codon[r["effect"]][r["codon_number"] - 1] += 1

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2), gridspec_kw={"width_ratios": [1, 2.4]})
    ks = [k for k in colors if counts[k]]
    ax[0].bar(ks, [counts[k] for k in ks], color=[colors[k] for k in ks])
    ax[0].set(title="All possible point mutations", ylabel="count")
    bottom = [0] * n_codons
    x = range(1, n_codons + 1)
    for k in colors:
        ax[1].bar(x, per_codon[k], bottom=bottom, color=colors[k], label=k, width=1.0)
        bottom = [b + v for b, v in zip(bottom, per_codon[k])]
    ax[1].set(title="Mutation outcomes per codon (9 possible each)", xlabel="Codon number", ylabel="mutations")
    ax[1].legend(ncol=4, fontsize=8)
    fig.suptitle(name, y=1.02)
    plt.tight_layout()
    plt.savefig("results/mutation_landscape.png", dpi=200, bbox_inches="tight")
    print("Saved results/mutations.csv and results/mutation_landscape.png")


if __name__ == "__main__":
    main()
