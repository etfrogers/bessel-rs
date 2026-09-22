#!/usr/bin/env python3
import json
import os
import subprocess

def get_commits_from_script():
    commits = []
    with open("bench_commits.sh") as f:
        in_commits = False
        for line in f:
            if "COMMITS=(" in line:
                in_commits = True
                continue
            if in_commits:
                if line.strip().startswith(")"):
                    break
                # Extract "commit" # comment
                parts = line.split("#")
                name_part = parts[0].strip().strip('"').strip("'")
                desc_part = parts[1].strip() if len(parts) > 1 else ""
                if name_part:
                    commits.append((name_part, desc_part))
    return commits

def format_time(t_ns):
    if t_ns is None:
        return "      N/A"
    if t_ns >= 1e6:
        return f"{t_ns/1e6:6.2f} ms"
    if t_ns >= 1e3:
        return f"{t_ns/1e3:6.2f} µs"
    return f"{t_ns:6.1f} ns"

def main():
    commits = get_commits_from_script()
    workloads = ["zeros_workload", "dense_r_workload", "sequence_workload"]

    print("=" * 115)
    print(" CHRONOLOGICAL BENCHMARK TIMELINE (Criterion Baselines)")
    print("=" * 115)

    for wl in workloads:
        print(f"\nWorkload: {wl}")
        print("-" * 115)
        print(f"{'#':2s}  {'Baseline':20s}  {'Time':10s}  {'vs v0.4':9s}  {'vs v1.0':9s}  {'vs Prev':9s}  {'Description':35s}")
        print("-" * 115)

        v04_time = None
        v1_time = None
        # Find abr-v0.4.0 and v1.0-release times for reference
        for i, (commit, _) in enumerate(commits):
            prefix_name = f"{i:02d}_{commit}"
            for candidate in [prefix_name, commit]:
                p = f"target/criterion/{wl}/{candidate}/estimates.json"
                if os.path.exists(p):
                    if "abr-v0.4.0" in commit and v04_time is None:
                        with open(p) as f:
                            v04_time = json.load(f)["mean"]["point_estimate"]
                    elif "v1.0-release" in commit and v1_time is None:
                        with open(p) as f:
                            v1_time = json.load(f)["mean"]["point_estimate"]
                    break

        prev_time = None
        for i, (commit, desc) in enumerate(commits):
            prefix_name = f"{i:02d}_{commit}"
            target_p = None
            for candidate in [prefix_name, commit]:
                p = f"target/criterion/{wl}/{candidate}/estimates.json"
                if os.path.exists(p):
                    target_p = p
                    break

            if target_p:
                with open(target_p) as f:
                    d = json.load(f)
                t_ns = d["mean"]["point_estimate"]
                t_str = format_time(t_ns)

                vs_v04 = f"{v04_time/t_ns:5.2f}x" if v04_time else "  -  "
                vs_v1 = f"{v1_time/t_ns:5.2f}x" if v1_time else "  -  "
                vs_prev = f"{prev_time/t_ns:5.2f}x" if prev_time else "  -  "

                # Highlight significant gains
                star = "★" if vs_v04 and v04_time and (v04_time/t_ns) >= 1.05 else " "
                print(f"{i:02d}  {commit:20s}  {t_str:10s}  {vs_v04:9s}  {vs_v1:9s}  {vs_prev:9s}  {desc[:35]:35s} {star}")
                prev_time = t_ns
            else:
                print(f"{i:02d}  {commit:20s}  {'[not run]':10s}  {'-':9s}  {'-':9s}  {'-':9s}  {desc[:35]:35s}")

    print("\n" + "=" * 115)

if __name__ == "__main__":
    main()
