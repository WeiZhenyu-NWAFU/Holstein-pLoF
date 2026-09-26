#!/usr/bin/env python3
"""Build publication-ready fixed-window summaries from jointly normalized iHS.

The input is the headered selscan norm output with columns:
LocusID Position Freq_Derived iHH1 iHH0 iHS_Raw iHS_Norm Is_Sig

The script streams the SNP file once, summarizes several window sizes, and
computes empirical tail probabilities after stratifying windows by SNP count.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--windows", nargs="+", type=int, default=[50000, 100000, 200000])
    parser.add_argument("--extreme-threshold", type=float, default=2.0)
    parser.add_argument("--min-snps", type=int, default=20)
    parser.add_argument("--snp-bins", type=int, default=10)
    parser.add_argument("--max-records", type=int, default=None, help="Testing only")
    return parser.parse_args()


def chrom_key(chrom: str) -> tuple[int, str]:
    value = chrom.removeprefix("chr").removeprefix("Chr")
    try:
        return int(value), chrom
    except ValueError:
        return 10**9, chrom


def conservative_tail_p(sorted_values: list[float], value: float) -> float:
    """P(X >= value), including every tie at the observed value."""
    return (len(sorted_values) - bisect.bisect_left(sorted_values, value)) / len(sorted_values)


def assign_density_bins(rows: list[dict], n_bins: int) -> None:
    ordered = sorted(range(len(rows)), key=lambda i: (rows[i]["N_SNP"], i))
    n = len(ordered)
    for rank, row_i in enumerate(ordered):
        rows[row_i]["SNPCountBin"] = min(n_bins, (rank * n_bins) // n + 1)


def add_empirical_statistics(rows: list[dict], min_snps: int, n_bins: int) -> dict:
    valid = [row for row in rows if row["N_SNP"] >= min_snps]
    if not valid:
        return {"valid_windows": 0, "top_1pct": 0, "top_5pct": 0}

    assign_density_bins(valid, n_bins)
    by_bin: dict[int, list[float]] = defaultdict(list)
    all_fractions = sorted(row["Frac_Extreme"] for row in valid)
    for row in valid:
        by_bin[row["SNPCountBin"]].append(row["Frac_Extreme"])
    for values in by_bin.values():
        values.sort()

    for row in rows:
        if row["N_SNP"] < min_snps:
            row.update(
                SNPCountBin="NA",
                DensityStratifiedTailP="NA",
                GenomewideTailP="NA",
                Top1pct=0,
                Top5pct=0,
            )
            continue
        tail_p = conservative_tail_p(by_bin[row["SNPCountBin"]], row["Frac_Extreme"])
        genome_p = conservative_tail_p(all_fractions, row["Frac_Extreme"])
        row["DensityStratifiedTailP"] = tail_p
        row["GenomewideTailP"] = genome_p
        row["Top1pct"] = int(tail_p <= 0.01)
        row["Top5pct"] = int(tail_p <= 0.05)

    return {
        "valid_windows": len(valid),
        "top_1pct": sum(row["Top1pct"] for row in rows),
        "top_5pct": sum(row["Top5pct"] for row in rows),
    }


def main() -> None:
    args = parse_args()
    windows = sorted(set(args.windows))
    if any(window <= 0 for window in windows):
        raise ValueError("Window sizes must be positive")
    if args.snp_bins <= 0:
        raise ValueError("--snp-bins must be positive")

    args.outdir.mkdir(parents=True, exist_ok=True)
    summaries: dict[int, dict[tuple[str, int], list]] = {window: {} for window in windows}
    seen_chromosomes: set[str] = set()
    processed = 0
    skipped = 0

    with args.input.open("rt", encoding="utf-8") as stream:
        header = stream.readline().strip().split()
        expected = ["LocusID", "Position", "Freq_Derived", "iHH1", "iHH0", "iHS_Raw", "iHS_Norm", "Is_Sig"]
        if header != expected:
            raise ValueError(f"Unexpected header: {header}")

        for line in stream:
            fields = line.split()
            if len(fields) < 7:
                skipped += 1
                continue
            locus = fields[0]
            chrom = locus.rsplit("_", 1)[0]
            try:
                pos = int(fields[1])
                norm = float(fields[6])
            except ValueError:
                skipped += 1
                continue
            if pos < 1 or not math.isfinite(norm):
                skipped += 1
                continue

            seen_chromosomes.add(chrom)
            abs_norm = abs(norm)
            is_extreme = int(abs_norm >= args.extreme_threshold)
            for window in windows:
                start = ((pos - 1) // window) * window + 1
                key = (chrom, start)
                stat = summaries[window].get(key)
                if stat is None:
                    # N, N_extreme, sum_abs, max_abs, max_locus, max_signed
                    summaries[window][key] = [1, is_extreme, abs_norm, abs_norm, locus, norm]
                else:
                    stat[0] += 1
                    stat[1] += is_extreme
                    stat[2] += abs_norm
                    if abs_norm > stat[3]:
                        stat[3] = abs_norm
                        stat[4] = locus
                        stat[5] = norm

            processed += 1
            if args.max_records is not None and processed >= args.max_records:
                break

    expected_chromosomes = {str(i) for i in range(1, 30)}
    if args.max_records is None and seen_chromosomes != expected_chromosomes:
        raise ValueError(
            f"Expected chromosomes 1-29; missing={sorted(expected_chromosomes-seen_chromosomes)}, "
            f"unexpected={sorted(seen_chromosomes-expected_chromosomes)}"
        )

    run_summary = {
        "input": str(args.input),
        "processed_snps": processed,
        "skipped_lines": skipped,
        "chromosomes": sorted(seen_chromosomes, key=chrom_key),
        "extreme_threshold_abs_iHS": args.extreme_threshold,
        "minimum_snps_for_ranking": args.min_snps,
        "snp_count_quantile_bins": args.snp_bins,
        "window_sizes_bp": windows,
        "testing_max_records": args.max_records,
        "outputs": {},
    }

    fieldnames = [
        "Chr", "WinStart_1based", "WinEnd_1based", "WindowSize_bp", "N_SNP",
        "N_Extreme", "Frac_Extreme", "MeanAbs_iHS", "MaxAbs_iHS",
        "MaxAbs_LocusID", "MaxAbs_iHS_Signed", "SNPCountBin",
        "DensityStratifiedTailP", "GenomewideTailP", "Top1pct", "Top5pct",
    ]

    for window in windows:
        rows = []
        for (chrom, start), stat in summaries[window].items():
            n_snp, n_extreme, sum_abs, max_abs, max_locus, max_signed = stat
            rows.append({
                "Chr": chrom,
                "WinStart_1based": start,
                "WinEnd_1based": start + window - 1,
                "WindowSize_bp": window,
                "N_SNP": n_snp,
                "N_Extreme": n_extreme,
                "Frac_Extreme": n_extreme / n_snp,
                "MeanAbs_iHS": sum_abs / n_snp,
                "MaxAbs_iHS": max_abs,
                "MaxAbs_LocusID": max_locus,
                "MaxAbs_iHS_Signed": max_signed,
            })
        rows.sort(key=lambda row: (chrom_key(str(row["Chr"])), row["WinStart_1based"]))
        empirical = add_empirical_statistics(rows, args.min_snps, args.snp_bins)

        output = args.outdir / f"joint_norm_iHS.{window // 1000}kb.windows.tsv"
        with output.open("wt", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

        run_summary["outputs"][str(window)] = {
            "file": str(output),
            "observed_windows": len(rows),
            **empirical,
        }

    with (args.outdir / "run_summary.json").open("wt", encoding="utf-8") as stream:
        json.dump(run_summary, stream, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
