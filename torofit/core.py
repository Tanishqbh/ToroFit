"""Core utilities for the ToroFit package."""


def version_info() -> dict:
    """Return version and package metadata.

    Returns
    -------
    dict
        A dictionary containing ``name`` and ``version`` keys.

    Examples
    --------
    >>> from torofit.core import version_info
    >>> info = version_info()
    >>> info["name"]
    'torofit'
    """
    from torofit import __version__

    return {"name": "torofit", "version": __version__}
