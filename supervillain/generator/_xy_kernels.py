"""Numba kernels for constrained XY single-site Metropolis sweeps.

The arrays passed here are deliberately plain, flat NumPy arrays.  Keeping
the lattice bookkeeping in :mod:`supervillain.generator.xy` makes the kernels
usable as a small, independently testable numerical core.
"""

import math

import numba
from numba import prange


_TWO_PI = 2.0 * math.pi


@numba.njit(cache=True)
def _wrap(angle):
    """Map an angle to [-pi, pi), preserving already-principal tiny links."""
    if -math.pi <= angle and angle < math.pi:
        return angle
    remainder = angle % _TWO_PI
    return remainder - _TWO_PI if remainder >= math.pi else remainder


@numba.njit(cache=True)
def _canonical(angle):
    """Map an angle to [0, 2*pi)."""
    remainder = angle % _TWO_PI
    # Some libm implementations round a negative sub-ulp input to 2*pi.
    # Zero is the nearest point on the circle and preserves the half-open
    # spin domain without introducing a spurious branch-edge link.
    return 0.0 if remainder >= _TWO_PI else remainder


def _make_sweep(parallel):
    """Build a serial or colour-parallel sweep with the oracle signature."""

    @numba.njit(cache=True, parallel=parallel)
    def sweep(phi, plus, minus, color_sites, offsets, proposals, uniforms,
              kappa, epsilon):
        accepted_by_color = 0

        # Every colour has no nearest-neighbour pairs, so updates within a
        # colour are independent.  The outer loop is intentionally serial:
        # it is the barrier separating colours.
        for color in range(offsets.shape[0] - 1):
            start = offsets[color]
            stop = offsets[color + 1]
            accepted_here = 0

            for position in (prange(start, stop) if parallel else range(start, stop)):
                site = color_sites[position]
                old_value = phi[site]
                new_value = _canonical(old_value + proposals[site])

                admissible = True
                delta_s = 0.0
                for mu in range(plus.shape[0]):
                    forward = plus[mu, site]
                    backward = minus[mu, site]

                    old_forward = _wrap(phi[forward] - old_value)
                    new_forward = _wrap(phi[forward] - new_value)
                    old_backward = _wrap(old_value - phi[backward])
                    new_backward = _wrap(new_value - phi[backward])

                    if abs(new_forward) >= epsilon or abs(new_backward) >= epsilon:
                        admissible = False
                        break

                    # Both the incoming and outgoing oriented bond are
                    # included.  At N=2 they have the same endpoint but are
                    # still distinct bonds and therefore are counted twice.
                    delta_s += -kappa * (
                        math.cos(new_forward) - math.cos(old_forward)
                        + math.cos(new_backward) - math.cos(old_backward)
                    )

                accepted = False
                if admissible:
                    if delta_s <= 0.0:
                        accepted = True
                    else:
                        accepted = uniforms[site] < math.exp(-delta_s)

                if accepted:
                    phi[site] = new_value
                    accepted_here += 1

            accepted_by_color += accepted_here

        return accepted_by_color

    return sweep


# These share the pure-Python oracle's positional signature and mutate ``phi``
# in place, returning the number of accepted site proposals.
sweep_serial = _make_sweep(parallel=False)
sweep_parallel = _make_sweep(parallel=True)
