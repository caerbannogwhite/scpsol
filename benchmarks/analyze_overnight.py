"""Analyze overnight parameter sweep results.

Usage:
    python benchmarks/analyze_overnight.py benchmarks/sweep-overnight-train.csv [--top N]
"""
import csv
import sys
from collections import defaultdict

KNOWN_OPT = {
    'scpb1': 69, 'scpb2': 76, 'scpb3': 80, 'scpb4': 79, 'scpb5': 72,
    'scpnre1': 29, 'scpnre2': 30, 'scpnre3': 27, 'scpnre4': 28, 'scpnre5': 28,
    'scpnrf1': 14, 'scpnrf2': 15, 'scpnrf3': 14, 'scpnrf4': 14, 'scpnrf5': 13,
}

def analyze(csvfile, top_n=15):
    combos = defaultdict(list)
    with open(csvfile) as f:
        for row in csv.DictReader(f):
            if not row.get('status') or not row.get('time'):
                continue
            key = (row['branch'], row['cut_freq'], row['balas_freq'],
                   row['eta'], row['sb'])
            combos[key].append(row)

    instances = set()
    for entries in combos.values():
        for e in entries:
            instances.add(e['instance'])
    n_instances = len(instances)

    results = []
    for key, entries in combos.items():
        branch, cf, bf, eta, sb = key
        n = len(entries)
        total_time = 0
        total_nodes = 0
        total_gap = 0
        n_optimal = 0
        n_timeout = 0
        n_correct = 0
        for e in entries:
            t = float(e['time']) if e['time'] else 600
            total_time += t
            total_nodes += int(e['nodes']) if e['nodes'] else 0
            gap = float(e['gap']) if e['gap'] else 100
            total_gap += gap
            if e['status'] == 'Optimal':
                n_optimal += 1
            elif e['status'] == 'TimeLimit':
                n_timeout += 1
            inst = e['instance']
            if inst in KNOWN_OPT and e['primal']:
                if abs(float(e['primal']) - KNOWN_OPT[inst]) < 0.5:
                    n_correct += 1

        avg_gap = total_gap / n if n > 0 else 100
        label = f"{branch[:3]} cf={cf} bf={bf}"
        if branch == 'reliability':
            label += f" e={eta} s={sb}"

        results.append({
            'key': key, 'label': label,
            'n': n, 'total_time': total_time, 'total_nodes': total_nodes,
            'avg_gap': avg_gap, 'n_optimal': n_optimal, 'n_timeout': n_timeout,
            'n_correct': n_correct,
        })

    # Sort: most optimal first, then lowest avg gap, then lowest time
    results.sort(key=lambda r: (-r['n_optimal'], r['avg_gap'], r['total_time']))

    print(f"{'#':>3} {'Config':<35} | {'n':>3} {'opt':>3} {'to':>3} | "
          f"{'avg_gap':>8} {'tot_time':>9} {'tot_nodes':>10}")
    print("-" * 95)
    for i, r in enumerate(results[:top_n]):
        print(f"{i+1:>3} {r['label']:<35} | {r['n']:>3} {r['n_optimal']:>3} {r['n_timeout']:>3} | "
              f"{r['avg_gap']:>7.4f}% {r['total_time']:>8.1f}s {r['total_nodes']:>10}")

    # Print all results sorted by time (for those with max optimal count)
    max_opt = results[0]['n_optimal'] if results else 0
    print(f"\n=== All combos with {max_opt} optimal (sorted by avg gap, then time) ===")
    top_results = [r for r in results if r['n_optimal'] == max_opt]
    top_results.sort(key=lambda r: (r['avg_gap'], r['total_time']))
    for i, r in enumerate(top_results[:top_n]):
        print(f"  {i+1:>2}. {r['label']:<35} gap={r['avg_gap']:.4f}% "
              f"time={r['total_time']:.1f}s nodes={r['total_nodes']}")

    print(f"\n=== Top 5 recommendations ===")
    for i, r in enumerate(top_results[:5]):
        branch, cf, bf, eta, sb = r['key']
        cmd = f"--branch {branch} --cut-frequency {cf} --balas-frequency {bf}"
        if branch == 'reliability':
            cmd += f" --reliability-eta {eta} --reliability-sb {sb}"
        print(f"  {i+1}. {cmd}")
        print(f"     opt={r['n_optimal']}/{r['n']} gap={r['avg_gap']:.4f}% "
              f"time={r['total_time']:.1f}s nodes={r['total_nodes']}")

if __name__ == '__main__':
    csvfile = sys.argv[1] if len(sys.argv) > 1 else 'benchmarks/sweep-overnight-train.csv'
    top_n = 15
    if '--top' in sys.argv:
        idx = sys.argv.index('--top')
        top_n = int(sys.argv[idx + 1])
    analyze(csvfile, top_n)
