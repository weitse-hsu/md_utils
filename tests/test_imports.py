import md_utils


def test_package_importable():
    assert md_utils.__version__ is not None
