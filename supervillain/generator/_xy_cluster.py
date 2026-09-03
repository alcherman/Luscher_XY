"""Serial Numba Wolff kernel for the constrained XY model.

Cluster growth is sequential.  Independent chains, rather than sites within
one growing cluster, are the natural unit of parallel work.
"""

import math

import numba
import numpy as np

from supervillain.generator._xy_kernels import _canonical, _wrap


@numba.njit(cache=True)
def cluster_serial(phi, plus, minus, seed_site, reflection_angle, uniforms,
                   kappa, epsilon):
    """Grow one cluster, mutate flat ``phi``, and return size and bond tests.

    Callers supply validated canonical float64 fields, ferromagnetic kappa,
    and uniforms indexed by ``(direction, positive-bond origin)``.  Thus the
    two parallel bonds in a periodic extent-two direction remain independent.
    All cluster decisions use the original field.  No mutation takes place
    until the candidate passes a full, strict admissibility check.
    """
    volume = phi.size
    alpha = _canonical(reflection_angle)
    reflected = np.empty(volume, dtype=np.float64)
    projections = np.empty(volume, dtype=np.float64)
    for site in range(volume):
        reflected[site] = _canonical(2.0 * alpha + math.pi - phi[site])
        projections[site] = math.cos(phi[site] - alpha)

    included = np.zeros(volume, dtype=np.bool_)
    queue = np.empty(volume, dtype=np.int64)
    queue[0] = seed_site
    included[seed_site] = True
    head = 0
    size = 1
    bond_tests = 0
    while head < size:
        site = queue[head]
        head += 1
        for mu in range(plus.shape[0]):
            for orientation in range(2):
                if orientation == 0:
                    neighbor = plus[mu, site]
                    origin = site
                else:
                    neighbor = minus[mu, site]
                    origin = neighbor
                if included[neighbor]:
                    continue
                bond_tests += 1

                # Mathematically these two distances agree.  Test both to
                # avoid choosing a cluster boundary by traversal direction
                # when floating-point rounding differs near strict epsilon.
                forced = (
                    abs(_wrap(phi[neighbor] - reflected[site])) >= epsilon
                    or abs(_wrap(reflected[neighbor] - phi[site])) >= epsilon
                )
                product = projections[site] * projections[neighbor]
                exponent = 2.0 * (kappa * max(0.0, product))
                probability = -math.expm1(-exponent)
                if forced or uniforms[mu, origin] < probability:
                    included[neighbor] = True
                    queue[size] = neighbor
                    size += 1

    # Exact reflections preserve internal bonds.  Finite-precision angle
    # storage can nevertheless push a strict-cutoff bond onto its boundary.
    # Fail closed, without a tolerance or an unaccounted-for rejection step.
    for site in range(volume):
        value = reflected[site] if included[site] else phi[site]
        for mu in range(plus.shape[0]):
            neighbor = plus[mu, site]
            other = reflected[neighbor] if included[neighbor] else phi[neighbor]
            if abs(_wrap(other - value)) >= epsilon:
                raise FloatingPointError(
                    "Cluster reflection reached the strict admissibility "
                    "boundary through floating-point rounding; field unchanged."
                )

    for position in range(size):
        site = queue[position]
        phi[site] = reflected[site]
    return size, bond_tests
