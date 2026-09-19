"""Validate quoted #775 tables only. This script cannot enumerate a new width."""
from pathlib import Path
from fractions import Fraction
import json
import math
import resource
from time import perf_counter

ROOT = Path(__file__).resolve().parent

def main():
    start = perf_counter()
    data = json.loads((ROOT/'quoted-input.json').read_text())
    checks = []
    for s in data['systems']:
        L=s['L']; n=L*L; c=s['C']
        assert len(c)==3 and all(len(row)==n+1 for row in c)
        assert all(type(v) is int and v>=0 for row in c for v in row)
        residuals=[sum(row[k] for row in c)-math.comb(n,k) for k in range(n+1)]
        assert not any(residuals), (L,residuals)
        assert sum(map(sum,c)) == 1 << n
        checks.append({'L':L,'integer_columns_checked':n+1,'nonzero_residuals':0,
                       'quoted_rank_totals':list(map(sum,c))})
    # This is scenario arithmetic from quoted telemetry, NOT a run or a TM forecast.
    cost = Fraction('24661.427')*(1 << 13)
    result={'issue':816,'status':'BLOCKED_BEFORE_ENUMERATION','L7_produced':False,'L8_started':False,
            'blocker':'Required #775 rank/winding source not obtained; the cited TM is a design, not a retrieved validated executable.',
            'new_configurations_enumerated':0,'transfer_states':None,'transfer_peak_rss':None,'transfer_wall_seconds':None,
            'quoted_table_checks':checks,'integer_columns_checked':sum(x['integer_columns_checked'] for x in checks),
            'brute_force_work':{'L7_masks':1<<49,'L8_masks':1<<64,'L7_to_L6_factor':1<<13,
               'quoted_L6_seconds':'24661.427','constant_rate_L7_seconds':str(cost),
               'constant_rate_L7_days':str(cost/86400),'caveat':'Illustrative same-rate extrapolation only; no measured TM state count or runtime.'},
            'validation_only_wall_seconds':round(perf_counter()-start,6),
            'validation_only_peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'scope':'Only integer revalidation of quoted L3-L6 tables. No rank kernel, fit, new-width run, or full CI.'}
    (ROOT/'result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__ == '__main__':
    main()
