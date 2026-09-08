"""Frozen analytic contract checks, runnable against immutable base or candidate.

Expected powers of ten follow the concentration definition, without calling the
candidate or duplicating its log implementation. Absolute p tolerance 1e-12 and
nM relative tolerance 1e-12 cover binary64 arithmetic for these simple cases.
This is author testing, not a blind independent evaluation.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import pathlib
import platform
import sys


def run(module_path):
    spec = importlib.util.spec_from_file_location('assay_under_test', module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    cases = []

    def check(name, args, expected):
        try:
            actual = module.normalize(*args)
            ok = actual.get('status') == expected['status']
            for key, value in expected.items():
                if key in ('p_activity', 'value_nm'):
                    got = actual.get(key)
                    ok = ok and isinstance(got, (int, float)) and math.isfinite(got)
                    if ok:
                        ok = math.isclose(got, value, rel_tol=1e-12 if key == 'value_nm' else 0,
                                          abs_tol=1e-12 if key == 'p_activity' else 0)
                else:
                    ok = ok and actual.get(key) == value
            if expected['status'] == 'invalid':
                ok = ok and bool(actual.get('reason')) and not any(
                    key in actual for key in ('value_nm', 'p_activity'))
            cases.append({'case': name, 'pass': bool(ok), 'expected': expected, 'actual': actual})
        except Exception as error:
            cases.append({'case': name, 'pass': False, 'expected': expected,
                          'exception': type(error).__name__ + ': ' + str(error)})

    # Same concentration expressed through every required unit, exact p = 6.
    for unit, value in [('nM', 1000), ('uM', 1), ('µM', 1), ('μM', 1),
                        ('mM', .001), ('M', .000001)]:
        check('equivalent_unit_' + unit, (value, unit),
              {'status': 'ok', 'endpoint': 'IC50', 'value_nm': 1000,
               'p_activity': 6, 'relation_nm': '=', 'p_relation': '='})
    for rel, canon, inverse in [('=', '=', '='), ('<', '<', '>'), ('<=', '<=', '>='),
                                ('>', '>', '<'), ('>=', '>=', '<='),
                                ('≤', '<=', '>='), ('≥', '>=', '<=')]:
        check('inequality_' + rel, ('10', 'nM', rel),
              {'status': 'ok', 'value_nm': 10, 'p_activity': 8,
               'relation_nm': canon, 'p_relation': inverse})
    for endpoint in ['IC50', 'Ki', 'Kd']:
        check('endpoint_' + endpoint, (100, 'nM', '=', endpoint),
              {'status': 'ok', 'endpoint': endpoint, 'value_nm': 100, 'p_activity': 7})
    for value, unit, nm, p in [(1e-300, 'M', 1e-291, 300),
                               (1e-20, 'nM', 1e-20, 29),
                               (1e100, 'nM', 1e100, -91)]:
        check('finite_range_' + str(value), (value, unit),
              {'status': 'ok', 'value_nm': nm, 'p_activity': p})
    for name, args in [
        ('zero', (0, 'nM')), ('negative', (-1, 'nM')),
        ('nan', ('nan', 'nM')), ('positive_inf', ('inf', 'nM')),
        ('negative_inf', ('-inf', 'nM')), ('missing', (None, 'nM')),
        ('empty', ('', 'nM')), ('boolean', (True, 'nM')),
        ('nonnumeric', ('potent', 'nM')), ('unknown_unit', (10, 'pM')),
        ('unit_missing', (10, None)), ('unit_unhashable', (10, [])),
        ('unknown_relation', (10, 'nM', '~')), ('relation_missing', (10, 'nM', None)),
        ('relation_unhashable', (10, 'nM', [])),
        ('unknown_endpoint', (10, 'nM', '=', 'EC50')),
        ('endpoint_missing', (10, 'nM', '=', None)),
        ('endpoint_unhashable', (10, 'nM', '=', [])),
        ('nm_overflow', (1e308, 'M')), ('float_overflow', (10**400, 'nM'))]:
        check(name, args, {'status': 'invalid'})
    return {'module': str(pathlib.Path(module.__file__).resolve()),
            'sha256': hashlib.sha256(pathlib.Path(module.__file__).read_bytes()).hexdigest(),
            'python': sys.executable, 'python_version': platform.python_version(),
            'collected': len(cases), 'passed': sum(c['pass'] for c in cases),
            'failed': sum(not c['pass'] for c in cases), 'skipped': 0, 'xfail': 0,
            'cases': cases}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--module', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = run(args.module)
    # Baseline NaN/Infinity is represented in JSON as a diagnostic string.
    def safe(value):
        if isinstance(value, float) and not math.isfinite(value):
            return str(value)
        if isinstance(value, dict):
            return {k: safe(v) for k, v in value.items()}
        if isinstance(value, list):
            return [safe(v) for v in value]
        return value
    pathlib.Path(args.output).write_text(json.dumps(safe(result), indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['module', 'collected', 'passed', 'failed', 'skipped', 'xfail']}))
    sys.exit(1 if result['failed'] else 0)
