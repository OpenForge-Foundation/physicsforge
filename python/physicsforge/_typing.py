"""Shared type aliases for the physicsforge public API."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Union

import numpy as np
from numpy.typing import ArrayLike, NDArray

if TYPE_CHECKING:
    from physicsforge.vector import Vector3

#: A float64 NumPy array.
FloatArray = NDArray[np.float64]

#: Return type of every broadcasting law: a Python ``float`` when all inputs are
#: scalars, otherwise a ``float64`` array with the broadcast shape of the inputs.
FloatOrArray = float | FloatArray

#: Anything convertible to a 3-vector: a :class:`~physicsforge.vector.Vector3`,
#: a length-3 sequence, or a shape ``(3,)`` array.
VectorLike = Union["Vector3", Sequence[float], NDArray[np.floating]]

__all__ = ["ArrayLike", "FloatArray", "FloatOrArray", "VectorLike"]
