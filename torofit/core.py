"""Core utilities for the ToroFit package."""


def greet(name: str = "World") -> str:
    """Return a greeting message.

    Parameters
    ----------
    name : str, optional
        The name to greet. Defaults to ``"World"``.

    Returns
    -------
    str
        A greeting string.

    Examples
    --------
    >>> from torofit import greet
    >>> greet("Athlete")
    'Hello, Athlete! Welcome to ToroFit.'
    """
    return f"Hello, {name}! Welcome to ToroFit."
