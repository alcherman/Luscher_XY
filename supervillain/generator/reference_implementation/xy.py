"""Slow, independent reference routines for compact XY local updates.

These routines deliberately use scalar Python loops and recompute the entire
action after each proposal.  They are an oracle for optimized implementations,
not a production sampler.
"""

import math

import numpy as np


_TWO_PI = 2.0 * math.pi


def _validate_parameters(kappa, epsilon):
    if (isinstance(kappa, (bool, np.bool_))
            or not isinstance(kappa, (int, float, np.integer, np.floating))
            or not math.isfinite(float(kappa))):
        raise ValueError("kappa must be a finite real number.")
    if (isinstance(epsilon, (bool, np.bool_))
            or not isinstance(epsilon, (int, float, np.integer, np.floating))
            or not math.isfinite(float(epsilon))
            or float(epsilon) <= 0.0 or float(epsilon) > math.pi):
        raise ValueError("epsilon must be finite and satisfy 0 < epsilon <= pi.")


def _principal_difference(delta):
    """Return ``delta`` on [-pi, pi) without erasing tiny values."""
    if -math.pi <= delta < math.pi:
        return delta
    wrapped = delta % _TWO_PI
    return wrapped - _TWO_PI if wrapped >= math.pi else wrapped


def _canonical_angle(angle):
    """Map an angle to [0, 2*pi), including modulo-rounding edge cases."""
    wrapped = angle % _TWO_PI
    return 0.0 if wrapped >= _TWO_PI else wrapped


def _validate_spatial_phi(phi, *, require_finite=True):
    array = np.asarray(phi)
    if array.ndim < 1:
        raise ValueError("phi must be a D-dimensional spatial ndarray with D >= 1.")
    if not np.issubdtype(array.dtype, np.number) or np.iscomplexobj(array):
        raise TypeError("phi must be a real-valued numeric ndarray.")
    if any(side < 2 for side in array.shape):
        raise ValueError("every spatial extent of phi must be at least 2.")
    if require_finite and not np.isfinite(array).all():
        raise ValueError("phi must contain only finite values.")
    return array


def reference_action(phi, kappa, epsilon=math.pi):
    r"""Evaluate the constrained compact XY action with explicit loops.

    ``phi`` must be a spatial array with one axis per lattice direction.  This
    function intentionally does not accept the leading component axis used by
    :class:`supervillain.lattice.Form`.
    """
    _validate_parameters(kappa, epsilon)
    field = _validate_spatial_phi(phi, require_finite=False)
    if not np.isfinite(field).all():
        return math.inf
    shape = field.shape
    total = 0.0
    for site in np.ndindex(shape):
        old = float(field[site])
        for direction in range(field.ndim):
            neighbor = list(site)
            neighbor[direction] = (neighbor[direction] + 1) % shape[direction]
            delta = _principal_difference(float(field[tuple(neighbor)]) - old)
            if not math.isfinite(delta) or abs(delta) >= float(epsilon):
                return math.inf
            total += -float(kappa) * math.cos(delta)
    return total


def _shape_from_neighbor_table(plus, minus):
    plus_array = np.asarray(plus)
    minus_array = np.asarray(minus)
    if plus_array.ndim != 2 or plus_array.shape != minus_array.shape:
        raise ValueError("plus and minus must have the same shape (D, V).")
    D, volume = plus_array.shape
    if D < 1 or volume < 2:
        raise ValueError("neighbor tables require D >= 1 and V >= 2.")
    N = int(round(volume ** (1.0 / D)))
    if N < 2 or N ** D != volume:
        raise ValueError("neighbor-table volume must equal N**D for integer N >= 2.")
    if (not np.issubdtype(plus_array.dtype, np.integer)
            or not np.issubdtype(minus_array.dtype, np.integer)
            or np.any(plus_array < 0) or np.any(plus_array >= volume)
            or np.any(minus_array < 0) or np.any(minus_array >= volume)):
        raise ValueError("neighbor tables must contain site indices in [0, V).")
    return D, N, volume


def reference_sweep(phi, plus, minus, color_sites, offsets, proposals, uniforms,
                    kappa, epsilon=math.pi):
    r"""Perform one deliberately slow checkerboard Metropolis sweep in place.

    ``phi`` is a contiguous flat array of length ``V=N**D``.  ``color_sites``
    concatenates the sites for each independent color, while ``offsets`` holds
    their start positions and one final ``V``.  Proposal and uniform arrays
    are indexed by *site*, so checkerboard order does not change their random
    association.  A proposal is wrapped to ``[0, 2*pi)`` before evaluating the
    full action.
    """
    _validate_parameters(kappa, epsilon)
    D, N, volume = _shape_from_neighbor_table(plus, minus)
    field = np.asarray(phi)
    if field.ndim != 1 or field.shape[0] != volume:
        raise ValueError("phi must be a flat array with length V.")
    if not field.flags.c_contiguous:
        raise ValueError("phi must be C-contiguous.")
    if not np.issubdtype(field.dtype, np.floating):
        raise TypeError("phi must have a floating-point dtype for in-place updates.")
    if not np.isfinite(field).all():
        raise ValueError("phi must contain only finite values.")

    colors = np.asarray(color_sites)
    boundaries = np.asarray(offsets)
    delta = np.asarray(proposals)
    draws = np.asarray(uniforms)
    if (colors.ndim != 1 or colors.shape[0] != volume
            or not np.issubdtype(colors.dtype, np.integer)
            or np.any(colors < 0) or np.any(colors >= volume)
            or np.unique(colors).size != volume):
        raise ValueError("color_sites must be a permutation of the V site indices.")
    if (boundaries.ndim != 1 or boundaries.size < 2
            or not np.issubdtype(boundaries.dtype, np.integer)
            or boundaries[0] != 0 or boundaries[-1] != volume
            or np.any(boundaries[1:] < boundaries[:-1])):
        raise ValueError("offsets must be nondecreasing color boundaries from 0 to V.")
    if delta.shape != (volume,) or draws.shape != (volume,):
        raise ValueError("proposals and uniforms must each have length V.")
    if (not np.issubdtype(delta.dtype, np.number) or np.iscomplexobj(delta)
            or not np.isfinite(delta).all()):
        raise ValueError("proposals must be finite real numbers.")
    if (not np.issubdtype(draws.dtype, np.number) or np.iscomplexobj(draws)
            or not np.isfinite(draws).all() or np.any(draws < 0.0)
            or np.any(draws >= 1.0)):
        raise ValueError("uniforms must be finite real numbers in [0, 1).")

    spatial = field.reshape((N,) * D)
    current = reference_action(spatial, kappa, epsilon)
    accepted = 0
    for color in range(boundaries.size - 1):
        for position in range(int(boundaries[color]), int(boundaries[color + 1])):
            site = int(colors[position])
            previous = float(field[site])
            field[site] = _canonical_angle(previous + float(delta[site]))
            proposed = reference_action(spatial, kappa, epsilon)
            if not math.isfinite(proposed):
                probability = 0.0
            elif not math.isfinite(current):
                probability = 1.0
            else:
                change = proposed - current
                probability = 1.0 if change <= 0.0 else math.exp(-change)
            if float(draws[site]) < probability:
                current = proposed
                accepted += 1
            else:
                field[site] = previous
    return accepted
