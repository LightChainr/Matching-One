#!/usr/bin/env python3
"""Single score of the fixed negative target on the existing independent archive."""
import sys
sys.dont_write_bytecode = True
import csv
import gzip
from itertools import zip_longest
import math
from pathlib import Path
import platform
from runner import ROOT, FIXED, OLD_HEADER, NEW_HEADER, digest, mapped, read, require, save, utc

METRICS = ('deltaY', 'deltaE', 'deltaZ2', 'earlyY', 'lateY')
M = FIXED['N'] - FIXED['b']


def load_blocks():
    archive = read(ROOT / 'inputs/archive/run.json')
    replay = read(ROOT / 'data/run.json')
    require(replay['status'] == 'completed' and replay['new_random_prefixes'] == 0, 'incomplete replay')
    require(replay['source_manifest_sha256'] == digest(ROOT / 'inputs/archive/run.json'), 'source manifest differs')
    aa = sorted(archive['batches'], key=lambda r: r['batch'])
    rr = sorted(replay['batches'], key=lambda r: r['batch'])
    require([r['batch'] for r in aa] == [r['batch'] for r in rr] == list(range(14)), 'batch IDs differ')
    blocks, alignment = [], []
    for old, new in zip(aa, rr):
        batch = old['batch']
        require(old['seed'] == new['seed'] and old['dependency_group'] == new['dependency_group'], 'batch dependency mismatch')
        require(new['original_file'] == old['file'], 'source file mismatch')
        op, np = ROOT / 'inputs/archive' / old['file'], ROOT / 'data' / new['file']
        require(digest(op) == old['sha256'] and digest(np) == new['sha256'], 'input hash differs')
        cells, rank_counts, rank1, eligible = {}, [0, 0, 0], [0, 0], [0, 0]
        rows, candidates, rank0_mapped = 0, 0, 0
        with gzip.open(op, 'rt') as f, gzip.open(np, 'rt') as g:
            a, z = csv.reader(f), csv.reader(g)
            require(next(a) == OLD_HEADER and next(z) == NEW_HEADER, 'column mismatch')
            for line, (o, n) in enumerate(zip_longest(a, z), 2):
                label = 'batch %d line %d' % (batch, line)
                require(o is not None and n is not None and len(o) == 7 and len(n) == 7, label + ': row mismatch')
                v = list(map(int, n))
                require(mapped(o) == v[:5], label + ': first5 mismatch')
                j1, rank, dx, dy, c, e, candidate = v
                rows += 1
                rank_counts[rank] += 1
                if rank != 1:
                    require(v[2:] == [0, 0, -1, -1, -1], label + ': sentinel mismatch')
                    rank0_mapped += int(rank == 0)
                    continue
                require(1 <= j1 <= FIXED['b'] and 0 <= c <= M, label + ': rank1 bounds')
                require(math.gcd(dx, dy) == 1 and (dx > 0 or (dx == 0 and dy > 0)), label + ': direction convention')
                s = M - c
                require(0 <= e <= s * (s - 1) // 2 and candidate >= 0, label + ': pair count')
                group = int(j1 > FIXED['a'])
                rank1[group] += 1
                candidates += candidate
                if s == 0:
                    continue
                eligible[group] += 1
                cell = cells.setdefault((dx, dy, c), [[0, 0], [0, 0]])
                cell[group][0] += 1
                cell[group][1] += e
        require(rows == old['samples'] == new['samples'] == 10000, 'batch size differs')
        blocks.append(dict(batch=batch, cells=cells, rank_counts=rank_counts, rank1=rank1, eligible=eligible))
        alignment.append(dict(batch=batch, rows=rows, all_mapped_first5_match=True, rank0_J1_mapped_to_zero=rank0_mapped,
                              seed=old['seed'], dependency_group=old['dependency_group'], sum_nonlocal_candidates=candidates,
                              original_file=old['file'], original_sha256=old['sha256'],
                              replay_file=new['file'], replay_sha256=new['sha256']))
    return blocks, alignment


def estimate(blocks, omit=None):
    cells, risk, eligible, ranks = {}, [0, 0], [0, 0], [0, 0, 0]
    for block in blocks:
        if block['batch'] == omit:
            continue
        for j in range(2):
            risk[j] += block['rank1'][j]
            eligible[j] += block['eligible'][j]
        for j in range(3):
            ranks[j] += block['rank_counts'][j]
        for key, groups in block['cells'].items():
            cell = cells.setdefault(key, [[0, 0], [0, 0]])
            for j in range(2):
                for k in range(2):
                    cell[j][k] += groups[j][k]
    terms, weights, rows = [[] for _ in METRICS], [], []
    supported = [0, 0]
    counts = dict(common=0, early_only=0, late_only=0)
    for (dx, dy, c), groups in sorted(cells.items()):
        (ne, se), (nl, sl) = groups
        s = M - c
        common = bool(ne and nl)
        w = ne * nl / (ne + nl) if common else 0.0
        counts['common' if common else ('early_only' if ne else 'late_only')] += 1
        row = dict(direction=[dx, dy], c=c, safe_sites=s, common_support=common, overlap_weight=w,
                   early=dict(n=ne, sum_e=se, meanY=2 * se / (s * ne) if ne else None),
                   late=dict(n=nl, sum_e=sl, meanY=2 * sl / (s * nl) if nl else None))
        if common:
            # Integer numerator before division; no fabricated probe observations.
            de = (se * nl - sl * ne) / (ne * nl)
            dy_value = 2 * de / s
            vector = [dy_value, de, -2 * de / (M * (M - 1)), row['early']['meanY'], row['late']['meanY']]
            weights.append(w)
            for term, value in zip(terms, vector):
                term.append(w * value)
            row.update(dict(zip(METRICS[:3], vector[:3])))
            supported[0] += ne
            supported[1] += nl
        rows.append(row)
    total_w = math.fsum(weights)
    require(total_w > 0, 'no common support for deletion ' + str(omit))
    point = [math.fsum(term) / total_w for term in terms]
    point[2] = -2 * point[1] / (M * (M - 1))
    for row in rows:
        row['normalized_weight'] = row['overlap_weight'] / total_w
    coverage = {}
    for group, j in [('early', 0), ('late', 1), ('combined', None)]:
        r, e, s = (sum(values) if j is None else values[j] for values in (risk, eligible, supported))
        coverage[group] = dict(rank1_prefixes=r, eligible_prefixes=e, supported_prefixes=s,
                               unsupported_eligible_prefixes=e-s, no_safe_site_prefixes=r-e,
                               supported_fraction_of_eligible=s/e if e else None,
                               supported_fraction_of_rank1=s/r if r else None)
    support = dict(cell_counts=counts, overlap_weight_sum=total_w, coverage=coverage, rank_counts=ranks)
    return point, support, rows


def main():
    require(not (ROOT / 'results/result.json').exists(), 'score already exists; no repeat scoring')
    blocks, alignment = load_blocks()
    point, support, rows = estimate(blocks)
    deletions = []
    for batch in range(14):
        vector, info, cells = estimate(blocks, omit=batch)
        deletions.append(dict(omitted_batch=batch, point_vector=vector, support=info,
             common_cells=[dict(direction=r['direction'], c=r['c'], nE=r['early']['n'], nL=r['late']['n'],
                                weight=r['overlap_weight'], normalized_weight=r['normalized_weight'])
                           for r in cells if r['common_support']]))
    means = [math.fsum(d['point_vector'][i] for d in deletions) / 14 for i in range(5)]
    covariance = [[13 / 14 * math.fsum((d['point_vector'][i]-means[i]) * (d['point_vector'][j]-means[j])
                                      for d in deletions) for j in range(5)] for i in range(5)]
    se = [math.sqrt(covariance[i][i]) for i in range(5)]
    reference = read(ROOT / 'inputs/reference-result.json')
    compare = {}
    for i, name in enumerate(METRICS):
        ref_point, ref_se = reference['point_estimates']['RB_' + name], reference['standard_errors']['RB_' + name]
        uncertainty = math.hypot(ref_se, se[i])
        compare[name] = dict(reference=ref_point, reference_SE=ref_se, archive=point[i], archive_SE=se[i],
                             archive_minus_reference=point[i]-ref_point, independent_difference_SE=uncertainty,
                             difference_over_SE=(point[i]-ref_point)/uncertainty)
    result = dict(status='scored', created_utc=utc(), fixed=FIXED, m=M, metric_order=METRICS,
                  point_estimates=dict(zip(METRICS, point)), standard_errors=dict(zip(METRICS, se)),
                  primary=dict(**support, cells=rows),
                  joint_uncertainty=dict(method='aligned whole-batch delete-one; exact (D,c) support and overlap weights recomputed',
                       covariance_formula='(B-1)/B sum (theta_minus_batch - mean)(theta_minus_batch - mean)^T',
                       point_estimator='pooled plug-in; no jackknife bias correction', covariance=covariance,
                       metric_order=METRICS, delete_one_mean=means, delete_one=deletions),
                  aligned_batches=alignment,
                  batch_sufficient_statistics=[dict(batch=b['batch'], rank_counts=b['rank_counts'], rank1=b['rank1'], eligible=b['eligible'],
                      cells=[dict(direction=list(k[:2]), c=k[2], early=dict(n=v[0][0], sum_e=v[0][1]),
                                  late=dict(n=v[1][0], sum_e=v[1][1])) for k, v in sorted(b['cells'].items())]) for b in blocks],
                  independent_block_comparison=compare,
                  prediction=dict(fixed_sign='negative', observed_sign='negative' if point[0] < 0 else 'nonnegative',
                                  primary_delta_over_SE=point[0]/se[0]),
                  provenance=dict(source_record=read(ROOT / 'source-record.json'),
                      scorer_sha256=digest(Path(__file__)), runner_sha256=digest(ROOT / 'runner.py'),
                      reference_result_sha256=digest(ROOT / 'inputs/reference-result.json'),
                      replay_manifest_sha256=digest(ROOT / 'data/run.json'), command=[sys.executable, *sys.argv],
                      python=sys.version, platform=platform.platform()),
                  scope=['Fixed target readout of previously unmeasured e in an existing independent random block.',
                         'Other results from this archive were previously inspected. Not new prospective sampling or untouched holdout.',
                         'No new samples, no block pooling, no cutoff/feature/model search, no artificial probe observations.',
                         'deltaZ2 = -2 deltaE/[m(m-1)] is derived from the same configurations, not independent evidence.',
                         'Finite L512 count-clock readout; no continuum, fixed-label-time or full Markov conclusion.'])
    save(ROOT / 'results/result.json', result)
    lines = ['# 独立既有档案的固定目标读出', '', '早组减晚组；±为14个对齐整批删除的标准误。', '',
             '| 指标 | 本档案 | SE |', '|---|---:|---:|']
    for name, value, error in zip(METRICS, point, se):
        lines.append('| %s | %.15g | %.15g |' % (name, value, error))
    cmp = compare['deltaY']
    lines += ['', '前块 deltaY = %.15g ± %.15g；本块减前块 = %.15g ± %.15g（sqrt(SE1²+SE2²)）。' %
              (cmp['reference'], cmp['reference_SE'], cmp['archive_minus_reference'], cmp['independent_difference_SE']),
              '本块 deltaY/SE = %.6f；两块未合并。' % (point[0]/se[0]), '',
              '140000行按约定映射后，前五字段全部逐行对应。共同支持：%d/%d（%.6f%%）；共同格数%d。' %
              (support['coverage']['combined']['supported_prefixes'], support['rank_counts'][1],
               100*support['coverage']['combined']['supported_fraction_of_rank1'], support['cell_counts']['common']), '',
              '按exact(D,c)格直接累加整数n、sum_e；Ybar=2sum_e/[(N-b-c)n]；w=nE*nL/(nE+nL)。',
              '每次删除整批重新求共同支持和权重；保存5×5完整协方差、14次删除向量/支持/权重、整数格统计和来源依赖。', '',
              '这是另一独立随机块上此前未测过e的固定目标读出。该档案的其它结果早已读过，不能称全新前瞻采样或原封未看的holdout。',
              '没有新采样、旧census或全套测试；32行成本检查为旧batch0的重复，科学权重为0。',
              'deltaZ2为同几何公式的派生量；结果仅属于有限L512插入时钟，不作连续极限或完整Markov判决。', '']
    with (ROOT / 'results/RESULT.md').open('x') as f:
        f.write('\n'.join(lines))
    print({k: result[k] for k in ('point_estimates', 'standard_errors', 'prediction')})
    print(compare['deltaY'])


if __name__ == '__main__':
    main()
