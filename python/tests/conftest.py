"""Shared pytest configuration for the physicsforge Python test suite."""

import numpy as np
import pytest


@pytest.fixture
def rng() -> np.random.Generator:
    """Deterministic random generator so failures are reproducible."""
    return np.random.default_rng(20260101)
