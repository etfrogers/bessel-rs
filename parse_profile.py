#!/usr/bin/env python3
"""
Profile Parser and Text Reporter for Samply/Firefox profiles in bessel-rs.
Generates a human- and AI-readable text report from profile.json.
"""

import json
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict

def parse_profile(profile_path="profile.json", out_report="profile_summary.txt", out_folded="profile.folded"):
    if not os.path.exists(profile_path):
        print(f"Error: {profile_path} not found.")
        sys.exit(1)

    print(f"Loading {profile_path}...")
    with open(profile_path) as f:
        data = json.load(f)

    # Find the main benchmark thread (thread with the most samples or named 'workload')
    target_thread = None
    max_samples = 0
    for t in data.get("threads", []):
        s_count = len(t.get("samples", {}).get("time", []))
        if s_count > max_samples:
            max_samples = s_count
            target_thread = t

    if not target_thread:
        print("Error: No thread with samples found.")
        sys.exit(1)

    print(f"Analyzing thread '{target_thread.get('name')}' with {max_samples} samples...")

    frame_table = target_thread["frameTable"]
    func_table = target_thread["funcTable"]
    res_table = target_thread["resourceTable"]
    stack_table = target_thread["stackTable"]
    samples = target_thread["samples"]
    string_array = target_thread["stringArray"]

    # Map libraries
    libs = data.get("libs", [])
    workloads_lib_idx = None
    workloads_bin_path = None
    for idx, lib in enumerate(libs):
        if "workload" in lib.get("name", "").lower():
            workloads_lib_idx = idx
            workloads_bin_path = lib.get("path")
            break

    # Collect unique addresses from workloads binary
    addr_to_hex = {}
    if workloads_bin_path and os.path.exists(workloads_bin_path):
        for i in range(frame_table["length"]):
            func_idx = frame_table["func"][i]
            res_idx = func_table["resource"][func_idx]
            if res_idx is not None and res_idx != -1:
                lib_idx = res_table["lib"][res_idx]
                if lib_idx == workloads_lib_idx:
                    addr = frame_table["address"][i]
                    if addr is not None and addr not in addr_to_hex:
                        addr_to_hex[addr] = f"0x{0x100000000 + addr:x}"

    symbol_map = {}
    if addr_to_hex:
        unique_addrs = list(addr_to_hex.keys())
        hex_list = [addr_to_hex[a] for a in unique_addrs]
        print(f"Symbolizing {len(hex_list)} addresses with atos...")
        p1 = subprocess.Popen(
            ["atos", "-o", workloads_bin_path, "-l", "0x100000000"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        atos_out, _ = p1.communicate(input="\n".join(hex_list))

        p2 = subprocess.Popen(
            ["c++filt"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        demangled_out, _ = p2.communicate(input=atos_out)
        for addr, line in zip(unique_addrs, demangled_out.strip().split("\n")):
            symbol_map[addr] = line

    # Resolve frame names
    frame_names = []
    clean_frame_names = []
    for i in range(frame_table["length"]):
        func_idx = frame_table["func"][i]
        res_idx = func_table["resource"][func_idx]
        addr = frame_table["address"][i]
        raw_name = string_array[func_table["name"][func_idx]] if func_table["name"][func_idx] is not None else "unknown"

        if res_idx is not None and res_idx != -1:
            lib_idx = res_table["lib"][res_idx]
            if lib_idx == workloads_lib_idx and addr in symbol_map:
                full_sym = symbol_map[addr]
                frame_names.append(full_sym)
                clean_sym = re.sub(r" \(in workloads-[^\)]+\).*", "", full_sym)
                clean_frame_names.append(clean_sym)
                continue
            lib_name = libs[lib_idx]["name"] if lib_idx < len(libs) else "unknown"
            frame_names.append(f"[{lib_name}] {raw_name}")
            clean_frame_names.append(f"[{lib_name}] {raw_name}")
        else:
            frame_names.append(raw_name)
            clean_frame_names.append(raw_name)

    stack_frames = stack_table["frame"]
    stack_prefixes = stack_table["prefix"]
    sample_stacks = samples["stack"]
    total_samples = len(sample_stacks)

    # Process samples
    folded_counter = Counter()
    self_counts = Counter()
    total_counts = Counter()
    line_counts = Counter()
    workload_samples = Counter()
    workload_self = defaultdict(Counter)

    for s in sample_stacks:
        if s is None:
            continue
        curr = s
        stack_full = []
        stack_clean = []
        while curr is not None:
            f_idx = stack_frames[curr]
            stack_full.append(frame_names[f_idx])
            stack_clean.append(clean_frame_names[f_idx])
            curr = stack_prefixes[curr]

        if not stack_full:
            continue

        # Reverse stack for folded format (root -> leaf)
        stack_clean_reversed = list(reversed(stack_clean))
        folded_line = ";".join(stack_clean_reversed)
        folded_counter[folded_line] += 1

        # Leaf = top of stack (exclusive self time)
        leaf_full = stack_full[0]
        leaf_clean = stack_clean[0]
        self_counts[leaf_clean] += 1
        line_counts[leaf_full] += 1

        # Check which workload this sample belongs to
        wl_found = "other"
        for fn in stack_clean:
            if "run_zeros_workload" in fn:
                wl_found = "zeros_workload"
                break
            elif "run_dense_r_workload" in fn:
                wl_found = "dense_r_workload"
                break
            elif "run_sequence_workload" in fn:
                wl_found = "sequence_workload"
                break
        workload_samples[wl_found] += 1
        workload_self[wl_found][leaf_clean] += 1

        for fn in set(stack_clean):
            total_counts[fn] += 1

    # Write folded stacks text file
    print(f"Writing folded stacks to {out_folded}...")
    with open(out_folded, "w") as f:
        for stack, count in folded_counter.most_common():
            f.write(f"{stack} {count}\n")

    # Generate text summary report
    print(f"Writing text report to {out_report}...")
    lines = []
    lines.append("=" * 80)
    lines.append(f"PROFILE TEXT REPORT: {profile_path}")
    lines.append(f"Total Samples: {total_samples}")
    lines.append("=" * 80)

    lines.append("\n" + "-" * 80)
    lines.append("1. CRITERION WORKLOAD SAMPLE DISTRIBUTION")
    lines.append("-" * 80)
    for wl, count in workload_samples.most_common():
        pct = (count / total_samples) * 100
        lines.append(f"{pct:6.2f}% ({count:5d} samples) : {wl}")

    lines.append("\n" + "-" * 80)
    lines.append("2. TOP 30 FUNCTIONS BY EXCLUSIVE (SELF) CPU TIME")
    lines.append("-" * 80)
    lines.append(f"{'Self %':>8} {'Samples':>8} | Function Name")
    lines.append("-" * 80)
    for fn, count in self_counts.most_common(30):
        pct = (count / total_samples) * 100
        lines.append(f"{pct:7.2f}% {count:8d} | {fn}")

    lines.append("\n" + "-" * 80)
    lines.append("3. TOP 25 HOTTEST SOURCE CODE LINES (Demangled with File & Line)")
    lines.append("-" * 80)
    lines.append(f"{'Self %':>8} {'Samples':>8} | Location")
    lines.append("-" * 80)
    for loc, count in line_counts.most_common(25):
        pct = (count / total_samples) * 100
        # clean location display
        loc_display = re.sub(r" \(in workloads-[^\)]+\)", "", loc)
        lines.append(f"{pct:7.2f}% {count:8d} | {loc_display}")

    for wl in ["zeros_workload", "dense_r_workload", "sequence_workload"]:
        wl_total = workload_samples[wl]
        if wl_total == 0:
            continue
        lines.append("\n" + "-" * 80)
        lines.append(f"4. TOP 10 HOTTEST FUNCTIONS IN {wl.upper()} ({wl_total} samples)")
        lines.append("-" * 80)
        lines.append(f"{'Workload %':>11} {'Samples':>8} | Function Name")
        lines.append("-" * 80)
        for fn, count in workload_self[wl].most_common(10):
            pct = (count / wl_total) * 100
            lines.append(f"{pct:10.2f}% {count:8d} | {fn}")

    lines.append("\n" + "-" * 80)
    lines.append("5. TOP 25 FUNCTIONS BY INCLUSIVE (TOTAL) TIME")
    lines.append("-" * 80)
    lines.append(f"{'Total %':>8} {'Self %':>8} {'Samples':>8} | Function Name")
    lines.append("-" * 80)
    amos_total = [(fn, count) for fn, count in total_counts.items() if "amos_bessel" in fn]
    amos_total.sort(key=lambda x: x[1], reverse=True)
    for fn, count in amos_total[:25]:
        tot_pct = (count / total_samples) * 100
        self_c = self_counts[fn]
        self_pct = (self_c / total_samples) * 100
        lines.append(f"{tot_pct:7.2f}% {self_pct:7.2f}% {self_c:8d} | {fn}")

    report_text = "\n".join(lines)
    with open(out_report, "w") as f:
        f.write(report_text)
    print("Done!")

if __name__ == "__main__":
    parse_profile()
