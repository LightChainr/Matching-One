#!/usr/bin/env python3
"""Freeze existing moments and existing64 tails; never produce trajectories."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parent
OLD = '375cd3a12b2b7a87d79148a59f62b95898f9e471'
OLD_PATH = 'results/p334-exact-score-quartet-moments'
EXISTING = 'experiments/p334-mechanism-response-20260831'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, default=ROOT.parent.parent)
    args = parser.parse_args()
    dest = ROOT/'inputs'
    dest.mkdir(exist_ok=False)
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=args.repo, text=True).strip()
    files = []

    def copy(commit, source, target, expected=None):
        data = subprocess.check_output(['git', 'show', f'{commit}:{source}'], cwd=args.repo)
        sha = hashlib.sha256(data).hexdigest()
        if expected is not None:
            assert sha == expected, (source, sha, expected)
        path = dest/target
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
        files.append({'commit': commit, 'source_path': source, 'local_path': target,
                      'sha256': sha, 'bytes': len(data)})
        return data

    meta = json.loads(copy(OLD, OLD_PATH+'/metadata.json', 'old8/metadata.json'))
    for n in (325, 425):
        copy(OLD, f'{OLD_PATH}/N{n}.npz', f'old8/N{n}.npz', meta['sizes'][str(n)]['sha256'])
    receipt = json.loads(copy(head, EXISTING+'/extension/run_receipt.json', 'new64/run_receipt.json'))
    assert receipt['prefixes'] == 3053 and receipt['new_tail_paths'] == 781568
    for entry in receipt['files']:
        copy(head, EXISTING+'/extension/'+entry['name'], 'new64/'+entry['name'], entry['sha256'])
    copy(head, EXISTING+'/results-extension/score.json', 'existing_local_mean_rank.json')
    copy(head, EXISTING+'/src/cell00_extension.cpp', 'provenance/cell00_extension.cpp')
    copy(OLD, 'scripts/p334_exact_score_quartet_moments.py', 'provenance/old8_reader.py')
    copy('dc4bb041', 'notes/p334-birth-covariance-hierarchy.md', 'provenance/completed_hierarchy.md')
    copy('03603388', 'notes/p334-nested-covariance-response.md', 'provenance/completed_identities.md')
    manifest = {'schema': 'p334.conditional-shape-inputs.v1', 'new_tail_source_head': head,
                'old8_archive_commit': OLD, 'new_sampling': 0,
                'plan_sha256_before_scoring': hashlib.sha256((ROOT/'PLAN.md').read_bytes()).hexdigest(),
                'files': files}
    (dest/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'files': len(files), 'bytes': sum(f['bytes'] for f in files),
                      'new_tail_source_head': head, 'new_samples': 0}))


if __name__ == '__main__':
    main()
