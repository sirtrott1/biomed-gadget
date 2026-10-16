#!/usr/bin/env python3
"""Compare certified physical fits with a floating point constrained optimizer."""
import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import LinearConstraint, minimize


def matrix(declaration):
    out = np.zeros((declaration['rows'], declaration['columns']))
    for entry in declaration['entries']:
        out[entry['row'], entry['column']] = float(Fraction(entry['value']))
    return out


def compare(root):
    physical = json.loads((root/'certificate.json').read_text())['result']['physical_shape']
    records = []
    for fit in physical['results']:
        preparation = physical['preparations'][fit['preparation']]
        family = physical['families'][preparation['family']]
        a, g = matrix(family['normals']), matrix(family['gram'])
        b = np.array([float(Fraction(v)) for v in family['rhs']])
        v0 = np.array([float(Fraction(v)) for v in fit['unconstrained_coordinates']])
        exact = np.array([float(Fraction(v)) for v in fit['coordinates']])
        if not len(v0):
            records.append(dict(state=fit['state'], response=fit['response'], model=fit['model'], dimension=0))
            continue
        # Unit changes for the comparison optimizer, with its full objective and
        # all inequalities retained. The exact source/certificate is unchanged.
        units = max(float(np.max(np.abs(v0))), float(np.max(np.abs(exact))), np.finfo(float).tiny)
        gscale = float(np.max(np.abs(g)))
        form = g/gscale
        center = v0/units
        objective = lambda z: float((z-center) @ form @ (z-center))
        gradient = lambda z: 2*form @ (z-center)
        rowscale = np.maximum(np.max(np.abs(a*units), axis=1), np.abs(b)) if len(b) else np.array([])
        rowscale[rowscale == 0] = 1
        constraints = [LinearConstraint(a*units/rowscale[:, None], -np.inf, b/rowscale)] if len(b) else []
        answer = minimize(objective, np.zeros(len(v0)), jac=gradient, constraints=constraints,
                          method='SLSQP', options=dict(ftol=1e-12, maxiter=10000))
        field_error = float(np.max(np.abs(matrix(family['frame']) @ (answer.x*units-exact))))
        relative_error = float(np.max(np.abs(answer.x-exact/units)))
        objective_error = abs(objective(answer.x)-objective(exact/units))
        if not answer.success or relative_error > 2e-5 or objective_error > 2e-8:
            raise ValueError((fit['response'], fit['model'], answer.message, relative_error, objective_error))
        records.append(dict(state=fit['state'], response=fit['response'], model=fit['model'],
                            iterations=int(answer.nit), maximum_original_field_error=field_error,
                            relative_coordinate_error=relative_error, scaled_objective_error=objective_error))
    return dict(status='PASS', fits=len(records), reference='SciPy SLSQP on the same full constrained quadratic objectives',
                scope='Numerical comparison; exact feasibility and optimality are verified separately', records=records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(compare(args.run.resolve()), indent=2))
