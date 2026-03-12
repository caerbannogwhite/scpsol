import csv
import sys
from collections import defaultdict

KNOWN_OPT = {
    'scp51': 253, 'scp52': 302, 'scp53': 226, 'scp54': 242, 'scp55': 211,
    'scp56': 213, 'scp57': 293, 'scp58': 288, 'scp59': 279, 'scp510': 265,
    'scpa1': 253, 'scpa2': 252, 'scpa3': 232, 'scpa4': 234, 'scpa5': 236,
    'scpb1': 69, 'scpb2': 76, 'scpb3': 80, 'scpb4': 79, 'scpb5': 72,
}

def analyze(csvfile):
    combos = defaultdict(list)
    with open(csvfile) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if not row['status'] or not row['primal'] or not row['time']:
                continue
            key = (row['cut_freq'], row['balas_freq'])
            combos[key].append({
                'instance': row['instance'],
                'primal': float(row['primal']),
                'time': float(row['time']),
                'nodes': int(row['nodes']),
            })

    instances = set()
    for entries in combos.values():
        for e in entries:
            instances.add(e['instance'])
    n_instances = len(instances)

    results = []
    for (cf, bf), entries in sorted(combos.items()):
        n = len(entries)
        total_time = sum(e['time'] for e in entries)
        total_nodes = sum(e['nodes'] for e in entries)
        n_optimal = sum(1 for e in entries if e['instance'] in KNOWN_OPT
                        and abs(e['primal'] - KNOWN_OPT[e['instance']]) < 0.5)
        n_subopt = sum(1 for e in entries if e['instance'] in KNOWN_OPT
                       and e['primal'] > KNOWN_OPT[e['instance']] + 0.5)
        subopt_instances = [e['instance'] for e in entries if e['instance'] in KNOWN_OPT
                            and e['primal'] > KNOWN_OPT[e['instance']] + 0.5]
        results.append({
            'cf': cf, 'bf': bf,
            'n': n, 'total_time': total_time, 'total_nodes': total_nodes,
            'n_optimal': n_optimal, 'n_subopt': n_subopt,
            'subopt': subopt_instances,
        })

    # Sort by: all optimal first, then by total time
    results.sort(key=lambda r: (-r['n_optimal'], r['n_subopt'], r['total_time']))

    print(f"{'cf':>4} {'bf':>4} | {'runs':>4} {'opt':>3} {'sub':>3} | {'time':>8} {'nodes':>7} | suboptimal")
    print("-" * 75)
    for r in results:
        subopt_str = ', '.join(r['subopt']) if r['subopt'] else ''
        print(f"{r['cf']:>4} {r['bf']:>4} | {r['n']:>4} {r['n_optimal']:>3} {r['n_subopt']:>3} | "
              f"{r['total_time']:>8.1f} {r['total_nodes']:>7} | {subopt_str}")

    # Top 5
    print(f"\n=== Top 5 (by correctness, then time) ===")
    for i, r in enumerate(results[:5]):
        print(f"  {i+1}. cf={r['cf']} bf={r['bf']}  time={r['total_time']:.1f}s  nodes={r['total_nodes']}  "
              f"opt={r['n_optimal']}/{r['n']}")

if __name__ == '__main__':
    analyze(sys.argv[1])
