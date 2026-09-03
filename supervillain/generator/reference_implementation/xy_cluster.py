"""Deliberately slow, full-bond-graph oracle for constrained XY clusters.

Unlike the production growth kernel, this implementation first constructs the
entire random bond graph.  Bond energies are evaluated from angle differences,
without the production kernel's spin-projection identity.  The seed's connected
component is then reflected.  Inputs follow the low-level kernel contract;
the public ``ClusterUpdate`` validates fields and random tapes.
"""

import math

import numpy as np


_TWO_PI = 2.0 * math.pi


def _distance(left, right):
    delta = float(left) - float(right)
    if -math.pi <= delta < math.pi:
        return abs(delta)
    delta %= _TWO_PI
    return abs(delta - _TWO_PI if delta >= math.pi else delta)


def reference_cluster(phi, plus, minus, seed_site, reflection_angle, uniforms,
                      kappa, epsilon):
    """Reflect the seed's component in place; return size and growth tests.

    ``uniforms[mu, x]`` belongs to the unique positive bond ``x -> x+mu``.
    Thus the two parallel bonds of a periodic length-two direction retain
    separate activation draws.  The reported test count refers to boundary
    bonds encountered while growing the seed component, not the oracle's
    additional whole-graph work.

    Both possible single-end reflections are checked for strict admissibility,
    symmetrically at floating-point boundaries.  If rounding makes the final
    reflected field invalid, fail before mutating the supplied field rather
    than silently introduce a state-dependent rejection probability.
    """
    volume = len(phi)
    dimension = plus.shape[0]
    reflection_angle = float(reflection_angle) % _TWO_PI
    reflected = np.empty(volume, dtype=np.float64)
    for site in range(volume):
        angle = (2.0 * reflection_angle + math.pi - float(phi[site])) % _TWO_PI
        reflected[site] = 0.0 if angle >= _TWO_PI else angle

    active = np.zeros((dimension, volume), dtype=bool)
    for mu in range(dimension):
        for site in range(volume):
            neighbor = int(plus[mu, site])
            if (_distance(reflected[site], phi[neighbor]) >= epsilon
                    or _distance(phi[site], reflected[neighbor]) >= epsilon):
                active[mu, site] = True
            else:
                original_angle = float(phi[site]) - float(phi[neighbor])
                reflected_angle = reflected[site] - float(phi[neighbor])
                # cos(a)-cos(b) without subtracting nearly equal cosines.
                # That cancellation can erase a small positive bond cost
                # which becomes order one at a large, finite kappa.
                difference = (-2.0 * math.sin((original_angle + reflected_angle) / 2.0)
                              * math.sin((original_angle - reflected_angle) / 2.0))
                increase = kappa * difference
                probability = -math.expm1(-max(0.0, increase))
                active[mu, site] = float(uniforms[mu, site]) < probability

    component = {int(seed_site)}
    queue = [int(seed_site)]
    tests = 0
    for site in queue:
        for mu in range(dimension):
            # This order only determines the work count; the full graph's
            # connected component is independent of traversal order.
            for neighbor, origin in ((int(plus[mu, site]), site),
                                     (int(minus[mu, site]), int(minus[mu, site]))):
                if neighbor in component:
                    continue
                tests += 1
                if active[mu, origin]:
                    component.add(neighbor)
                    queue.append(neighbor)

    candidate = np.array(phi, dtype=np.float64, copy=True)
    for site in component:
        candidate[site] = reflected[site]
    for mu in range(dimension):
        for site in range(volume):
            if _distance(candidate[site], candidate[int(plus[mu, site])]) >= epsilon:
                raise FloatingPointError(
                    "Cluster reflection lost strict admissibility through rounding"
                )
    phi[:] = candidate
    return len(component), tests
