import numpy as np
import pytest
import scipy.sparse as sp

import scrna_toolkit as sct

ALL_FORMATS = ["csr", "csc", "coo", "bsr", "lil", "dok", "dia"]


@pytest.fixture(autouse=True)
def _no_leftover_guard():
    yield
    sct.remove_dense_guard()


def _matrix(fmt: str, shape=(30, 40), array: bool = False):
    m = sp.random(*shape, density=0.1, format=fmt, random_state=0)
    return getattr(sp, f"{fmt}_array")(m) if array else m


@pytest.mark.parametrize("array", [False, True])
@pytest.mark.parametrize("fmt", ALL_FORMATS)
def test_guard_blocks_large_and_allows_small(fmt, array):
    big, small = _matrix(fmt, array=array), _matrix(fmt, (5, 5), array=array)
    expected = big.toarray()

    with sct.dense_guard(max_elements=100):
        with pytest.raises(sct.DenseConversionError, match="30 x 40"):
            big.toarray()
        if not array:
            with pytest.raises(sct.DenseConversionError):
                big.todense()
        assert small.toarray().shape == (5, 5)

    np.testing.assert_array_equal(big.toarray(), expected)  # removed on exit


def test_guard_error_is_a_memory_error():
    with sct.dense_guard(max_elements=10), pytest.raises(MemoryError):
        sp.random(10, 10, density=0.5, format="csr").toarray()


def test_nested_guards_restore_the_outer_limit():
    m = _matrix("csr")  # 1200 elements
    assert sct.dense_guard_limit() is None
    with sct.dense_guard(max_elements=5000):
        m.toarray()
        with sct.dense_guard(max_elements=100):
            assert sct.dense_guard_limit() == 100
            with pytest.raises(sct.DenseConversionError):
                m.toarray()
        assert sct.dense_guard_limit() == 5000
        m.toarray()
    assert sct.dense_guard_limit() is None


def test_install_and_remove_dense_guard():
    m = _matrix("csc")
    sct.install_dense_guard(max_elements=100)
    sct.install_dense_guard(max_elements=50)  # re-install only changes the limit
    with pytest.raises(sct.DenseConversionError, match="limit 50 elements"):
        m.toarray()
    sct.remove_dense_guard()
    sct.remove_dense_guard()  # idempotent
    assert m.toarray().shape == (30, 40)
    with pytest.raises(ValueError, match="positive"):
        sct.install_dense_guard(0)


def test_guard_also_stops_densifying_inside_the_toolkit(adata):
    sct.normalize_total(adata)
    sct.log1p(adata)
    with sct.dense_guard(max_elements=1000), pytest.raises(sct.DenseConversionError):
        sct.scale(adata)
