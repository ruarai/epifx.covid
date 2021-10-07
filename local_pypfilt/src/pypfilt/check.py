"""Provide convenience functions for checking invariants."""

def _expect_cols(ncol, exp, raise_exc=True):
    if ncol == exp:
        return True
    elif raise_exc:
        raise ValueError('Expected {} columns but found {}'.format(exp, ncol))
    else:
        return False


def is_entire_matrix(params, hist, raise_exc=True):
    """
    Check whether the history matrix includes all columns (including, e.g.,
    the particle weights and parent indices).

    :param params: The simulation parameters.
    :param hist: The history matrix.
    :param raise_exc: Whether to raise an exception if the check fails.

    :returns: ``True`` if the check is successful. If the check fails, either
        a ``ValueError`` exception is raised (if ``raise_exc == True``) or
        ``False`` is returned (if ``raise_exc == False``).
    """
    exp = params['hist']['state_cols'] + params['hist']['extra_cols']
    return _expect_cols(hist.shape[-1], exp, raise_exc)


def is_only_statevec(params, hist, raise_exc=True):
    """
    Check whether the history matrix contains only the particle state vector
    columns.

    :param params: The simulation parameters.
    :param hist: The history matrix.
    :param raise_exc: Whether to raise an exception if the check fails.

    :returns: ``True`` if the check is successful. If the check fails, either
        a ``ValueError`` exception is raised (if ``raise_exc == True``) or
        ``False`` is returned (if ``raise_exc == False``).
    """
    exp = params['hist']['state_cols']
    return _expect_cols(hist.shape[-1], exp, raise_exc)
