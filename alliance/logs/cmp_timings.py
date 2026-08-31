"""Compare per-recipe timings between two apero_processing logs.

Parses the "Timings:" summary block emitted by apero_processing, which
lists lines of the form::

    ID = 0    Time = 14.400
              apero_preprocess_spirou.py ... --shortname=PP --parallel=True

Aggregates total/mean wall time per recipe shortname and reports the
difference between the two logs.
"""
import re
import sys
from collections import defaultdict

ANSI = re.compile(r'\x1b\[[0-9;]*m')
TIME_RE = re.compile(r'\bID\s*=\s*(\d+)\s+Time\s*=\s*([0-9.]+)')
SHORT_RE = re.compile(r'--shortname=(\S+)')
RECIPE_RE = re.compile(r'(apero_\w+\.py)')


def clean(line):
    """Strip ANSI colour codes and the log prefix."""
    return ANSI.sub('', line).rstrip('\n').rstrip('\r')


def parse(path):
    """Return {shortname: [total_time, count]} for one log file."""
    agg = defaultdict(lambda: [0.0, 0])
    seen_ids = set()
    with open(path, 'rb') as fobj:
        raw = fobj.read().decode('utf-8', 'replace')
    lines = raw.split('\n')
    nlines = len(lines)
    for idx in range(nlines):
        line = clean(lines[idx])
        match = TIME_RE.search(line)
        if match is None:
            continue
        runid = int(match.group(1))
        secs = float(match.group(2))
        # look ahead for the runstring (usually the very next line)
        name = None
        for jdx in range(idx + 1, min(idx + 4, nlines)):
            nxt = clean(lines[jdx])
            if TIME_RE.search(nxt):
                break
            smatch = SHORT_RE.search(nxt)
            if smatch is not None:
                name = smatch.group(1)
                break
            rmatch = RECIPE_RE.search(nxt)
            if rmatch is not None:
                name = rmatch.group(1)
                break
        if name is None:
            name = 'UNKNOWN'
        key = (runid, name)
        if key in seen_ids:
            continue
        seen_ids.add(key)
        agg[name][0] += secs
        agg[name][1] += 1
    return agg


def main():
    path_a, path_b = sys.argv[1], sys.argv[2]
    agg_a = parse(path_a)
    agg_b = parse(path_b)

    tot_a = sum(v[0] for v in agg_a.values())
    tot_b = sum(v[0] for v in agg_b.values())
    cnt_a = sum(v[1] for v in agg_a.values())
    cnt_b = sum(v[1] for v in agg_b.values())

    print('%-16s %22s %22s %14s' % ('', 'v0.7', 'v0.8', 'delta'))
    print('%-16s %8s %6s %6s %8s %6s %6s %14s'
          % ('shortname', 'total', 'n', 'mean', 'total', 'n', 'mean',
             'total_delta'))
    print('-' * 88)

    names = sorted(set(agg_a) | set(agg_b),
                   key=lambda k: -(agg_b.get(k, [0, 0])[0]
                                   - agg_a.get(k, [0, 0])[0]))
    for name in names:
        ta, na = agg_a.get(name, [0.0, 0])
        tb, nb = agg_b.get(name, [0.0, 0])
        ma = ta / na if na else 0.0
        mb = tb / nb if nb else 0.0
        print('%-16s %8.0f %6d %6.1f %8.0f %6d %6.1f %+14.0f'
              % (name, ta, na, ma, tb, nb, mb, tb - ta))

    print('-' * 88)
    print('%-16s %8.0f %6d %6s %8.0f %6d %6s %+14.0f'
          % ('TOTAL', tot_a, cnt_a, '', tot_b, cnt_b, '', tot_b - tot_a))


if __name__ == '__main__':
    main()
