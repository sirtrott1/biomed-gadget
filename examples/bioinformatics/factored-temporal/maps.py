"""Independent expansion of complete factored certificates using SymPy rationals."""
import copy
import sympy as S


def integer(n):
    if type(n) is not int or n < 0:
        raise ValueError('invalid map coordinate')
    return n


def sparse(m):
    rows, cols = integer(m['rows']), integer(m['columns'])
    entries = {}
    for e in m['entries']:
        i, j = integer(e['row']), integer(e['column'])
        if i >= rows or j >= cols or (i, j) in entries:
            raise ValueError('duplicate or out-of-range coefficient')
        value = S.Rational(e['value'])
        if not value:
            raise ValueError('explicit zero in sparse matrix')
        entries[i, j] = value
    return S.MutableSparseMatrix(rows, cols, entries)


class Expansion:
    def __init__(self, archive):
        if archive['format'] != 'gadget-factored-maps-v1':
            raise ValueError('unsupported factored map format')
        self.nodes = archive['nodes']
        self.cache = {}

    def node(self, i):
        i = integer(i)
        if i >= len(self.nodes):
            raise ValueError('foreign expression reference')
        if i in self.cache:
            return self.cache[i]
        n = self.nodes[i]
        rows, cols = integer(n['rows']), integer(n['columns'])

        def previous(j):
            if integer(j) >= i:
                raise ValueError('factored map is cyclic or forward-referenced')
            return self.node(j)

        kind = n['kind']
        if kind == 'matrix':
            value = sparse(n)
        elif kind == 'identity':
            if rows != cols:
                raise ValueError('rectangular identity')
            value = S.eye(rows).as_mutable().todok()
            value = S.MutableSparseMatrix(rows, rows, value)
        elif kind == 'transpose':
            value = previous(n['input']).T
        elif kind == 'scale':
            value = S.Rational(n['coefficient']) * previous(n['input'])
        elif kind == 'sum':
            if not n['inputs']:
                raise ValueError('untyped empty sum')
            value = S.zeros(rows, cols)
            for j in n['inputs']:
                value += previous(j)
        elif kind == 'compose':
            if not n['application_order']:
                raise ValueError('untyped empty composition')
            value = S.eye(cols)
            for j in n['application_order']:
                value = previous(j) * value
        elif kind == 'solve':
            a = previous(n['system'])
            if a.rows != a.cols:
                raise ValueError('rectangular solve system')
            if a.is_diagonal():
                if any(a[j, j] == 0 for j in range(a.rows)):
                    raise ValueError('singular solve')
                value = S.diag(*[1 / a[j, j] for j in range(a.rows)])
            else:
                value = a.inv()
            if a * value != S.eye(rows):
                raise ValueError('nonzero complete inverse residual')
        else:
            raise ValueError('unsupported factored node')
        if value.shape != (rows, cols):
            raise ValueError('expression dimension mismatch')
        self.cache[i] = S.MutableSparseMatrix(value)
        return self.cache[i]

    def map(self, descriptor):
        value = self.node(descriptor['expression'])
        if value.shape != (descriptor['rows'], descriptor['columns']):
            raise ValueError('descriptor dimension mismatch')
        return value

    def all(self):
        for i in range(len(self.nodes)):
            self.node(i)


def encode(m):
    return dict(rows=int(m.rows), columns=int(m.cols),
                entries=[dict(row=int(i), column=int(j), value=str(v))
                         for (i, j), v in sorted(m.todok().items()) if v])


def expand_result(result):
    """Return the original enumerated result shape, preserving all other evidence."""
    def walk(value, expansion=None):
        if isinstance(value, dict):
            if 'factored_maps' in value:
                expansion = Expansion(value['factored_maps'])
                expansion.all()
            if 'expression' in value:
                if expansion is None:
                    raise ValueError('unbound factored descriptor')
                return encode(expansion.map(value))
            return {k: walk(v, expansion) for k, v in value.items() if k != 'factored_maps'}
        if isinstance(value, list):
            return [walk(v, expansion) for v in value]
        return copy.deepcopy(value)
    return walk(result)


def check_source_bindings(archive, endpoints):
    """Check complete metric/boundary/projection recipes against retained source."""
    exp = Expansion(archive)
    exp.all()
    systems = []
    for endpoint in endpoints:
        d = endpoint['declaration']
        grades = d['source_grades']
        boundaries = [S.zeros(0, len(grades[0]))]
        metrics = []
        for k, cells in enumerate(grades):
            spec = d['metrics'][k]
            metrics.append(S.eye(len(cells)) if spec['kind'] == 'identity' else
                           S.diag(*[S.Rational(v) for v in spec['values']]))
            if k:
                lo = {c['key']: i for i, c in enumerate(grades[k-1])}
                up = {c['key']: i for i, c in enumerate(cells)}
                b = S.zeros(len(lo), len(up))
                for o in d['source_occurrences']:
                    if o['grade'] == k:
                        b[lo[o['lower_key']], up[o['parent_key']]] += S.Rational(o['coefficient'])
                boundaries.append(b)
        adjoints = [S.zeros(len(grades[0]), 0)]
        for k in range(1, len(grades)):
            adjoints.append(metrics[k].inv() * boundaries[k].T * metrics[k-1])
        systems.append((boundaries, metrics, adjoints))
    checks = 0
    for i, n in enumerate(exp.nodes):
        if 'endpoint' not in n:
            continue
        t, k, role = n['endpoint'], n['grade'], n['role']
        bs, ms, ads = systems[t]
        actual = exp.node(i)
        h = ads[k] * bs[k]
        if k+1 < len(bs):
            h += bs[k+1] * ads[k+1]
        if role == 'metric':
            expected = ms[k]
        elif role == 'boundary':
            expected = bs[k]
        elif role == 'adjoint':
            expected = ads[k]
        elif role == 'hodge':
            expected = h
        elif role.startswith('green:'):
            # Compare the full equation rather than solving it again.
            assert (S.eye(h.rows) + S.Rational(role[6:]) * h) * actual == S.eye(h.rows)
            checks += actual.rows * actual.cols
            continue
        elif role.startswith('sector:'):
            sector = int(role[7:])
            frames = [ads[k], bs[k+1] if k+1 < len(bs) else S.zeros(ms[k].rows, 0)]
            ps = []
            for f in frames:
                independent = f.columnspace()
                if independent:
                    f = S.Matrix.hstack(*independent)
                    ps.append(f * (f.T * ms[k] * f).inv() * f.T * ms[k])
                else:
                    ps.append(S.zeros(ms[k].rows))
            expected = ps[sector] if sector < 2 else S.eye(ms[k].rows) - ps[0] - ps[1]
        else:
            raise ValueError('unrecognized source binding')
        assert actual == expected, (t, k, role)
        checks += actual.rows * actual.cols
    return checks
