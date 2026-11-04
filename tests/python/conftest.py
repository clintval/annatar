from pathlib import Path

import pytest

FIXTURE = Path(__file__).parent / "out" / "fixture.h5ad"


@pytest.fixture(scope="session")
def fixture_path() -> Path:
    if not FIXTURE.exists():
        pytest.fail(f"missing {FIXTURE}; run `pixi run fixture-h5ad` first")
    return FIXTURE
