"""Integer 4/8 boundary-convention regression, NOT a torus jump2 census."""
from pathlib import Path
import json
from time import perf_counter


def crosses(mask, w, h, eight, horizontal):
    steps = [(1,0),(-1,0),(0,1),(0,-1)]
    if eight:
        steps += [(1,1),(1,-1),(-1,1),(-1,-1)]
    start = [y*w for y in range(h)] if horizontal else list(range(w))
    seen = {i for i in start if mask >> i & 1}
    todo = list(seen)
    while todo:
        u = todo.pop(); x, y = u % w, u // w
        if (horizontal and x == w-1) or (not horizontal and y == h-1):
            return True
        for dx,dy in steps:
            a,b=x+dx,y+dy
            if 0 <= a < w and 0 <= b < h:
                v=b*w+a
                if mask >> v & 1 and v not in seen:
                    seen.add(v); todo.append(v)
    return False


def main():
    t=perf_counter(); checks=[]
    for w,h in [(2,2),(3,3),(4,4)]:
        total=1 << (w*h); bad=[]
        for m in range(total):
            black=crosses(m,w,h,False,True)
            white=crosses((total-1)^m,w,h,True,False)
            if black == white:
                bad.append(m)
        assert not bad, (w,h,bad[:1])
        checks.append({'w':w,'h':h,'configurations':total,'violations':len(bad)})
    # With 4/4 adjacency the 2x2 diagonal coloring has neither crossing.
    assert not crosses(9,2,2,False,True)
    assert not crosses(6,2,2,False,False)
    result={
        'issue':815,'deterministic_lemma':'proved for the explicitly declared lattice annuli',
        'acceptance':'partial: requested all-jump2 L<=5 census not executed',
        'rectangle_4_8_checks':checks,'wrong_4_4_control':{'black_mask':9,'white_mask':6},
        'all_jump2_configurations_checked':0,
        'blocker':'Named #775/#769 C kernel and raw witness files were not obtained from the checked branch sources; no replacement rank kernel was written.',
        'wall_seconds':round(perf_counter()-t,6)
    }
    Path(__file__).with_name('result.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__ == '__main__':
    main()
