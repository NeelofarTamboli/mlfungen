"""Compare step-0 usage estimates between two sets of McSplicer output directories.

Used to compare the normal (strand-unaware) run with the strand-split run, e.g.
    python -m src.validation.compare_runs --a out_all --b out_plus out_minus

Prints how many sites are in both, the largest and mean absolute differences, the sites
that differ most, and the sites present in only one run.
"""

from __future__ import annotations

import argparse

from src.validation.check_mcsplicer_output import read_outputs


def merged(dirs):
    out = {}
    for d in dirs:
        out.update({k: v[0] for k, v in read_outputs(d).items()})
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", nargs="+", required=True)
    ap.add_argument("--b", nargs="+", required=True)
    ap.add_argument("--top", type=int, default=6)
    a = ap.parse_args()
    A, B = merged(a.a), merged(a.b)
    both = sorted(set(A) & set(B))
    diffs = sorted(((abs(A[k] - B[k]), k) for k in both), reverse=True)
    print(f"sites: A {len(A)}, B {len(B)}, both {len(both)}")
    print(f"|A-B|: max {diffs[0][0]:.4f}, mean {sum(d for d, _ in diffs) / len(diffs):.5f}, >0.01 in {sum(d > 0.01 for d, _ in diffs)} sites")
    for d, k in diffs[: a.top]:
        print(f"  {k[0]} {k[1]} {k[3]} {k[2]}: A {A[k]:.4f} B {B[k]:.4f}")
    print("only in A:", [(k, round(A[k], 4)) for k in sorted(set(A) - set(B))])
    print("only in B:", [(k, round(B[k], 4)) for k in sorted(set(B) - set(A))])


if __name__ == "__main__":
    main()
