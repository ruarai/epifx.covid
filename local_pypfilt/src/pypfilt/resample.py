"""Various methods for resampling particles."""

import logging
import numpy as np

from . import context


def post_regularise(ctx, px, new_px):
    """
    Sample model parameter values from a continuous approximation of the
    optimal filter, assuming that it has a smooth density.

    This is the post-regularised particle filter (post-RPF). For details, see
    chapter 12 of Doucet et al., Sequential Monte Carlo Methods in Practice,
    Springer, 2001.
    `doi:10.1007/978-1-4757-3437-9_12
    <https://doi.org/10.1007/978-1-4757-3437-9_12>`_

    :param params: The simulation parameters.
    :param px: The particles, prior to resampling.
    :param new_px: The particles after resampling directly from the discrete
        distribution (``px``). This matrix will be **destructively updated**
        with model parameter values samples from the regularisation kernel.
    """

    from . import stats

    logger = logging.getLogger(__name__)

    rnd = ctx.component['random']['resample']
    count = px.shape[0]
    # Only resample parameters that can be sampled continuously.
    details = ctx.component['model'].describe()
    p_info = [(ix, name)
              for ix, (name, smooth, _vmin, _vmax) in enumerate(details)
              if smooth]
    p_ixs = np.array([info[0] for info in p_info])
    #print('Resample can smooth: {}'.format([info[1] for info in p_info]))

    if len(p_ixs) == 0:
        logger.debug("Post-RPF: no parameters to resample")
        return

    # Check for parameters that are constant (or nearly so) for all particles.
    # These parameters must be ignored or the covariance matrix will not be
    # positive definite, and the Cholesky decomposition will fail.
    p_range = np.ptp(px[:, p_ixs], axis=0)
    toln = ctx.params['resample']['reg_toln']
    good = p_range >= toln
    if not np.all(good):
        bad = np.logical_not(good)
        msg = "Post-RPF found {} constant parameter(s) at {}".format(
            sum(bad), p_ixs[np.where(bad)])
        logger.debug(msg)
        #print(p_range)
        #print(good)
        p_ixs = p_ixs[good]
        if len(p_ixs) == 0:
            logger.debug("Post-RPF: no non-constant parameters to resample")
            return

    # Use a bandwidth that is half that of the optimal bandwidth for a
    # Gaussian kernel (when the underlying density is Gaussian with unit
    # covariance), to handle multi-model densities.
    npar = len(p_ixs)
    h = 0.5 * (4 / (count * (npar + 2))) ** (1 / (npar + 4))

    #print(p_ixs)
    #print(px[:5, p_ixs])

    # Calculate the Cholesky decomposition of the parameter covariance
    # matrix V, which is used to transform independent normal samples
    # into multivariate normal samples with covariance matrix V.
    try:
        cov_mat = stats.cov_wt(px[:, p_ixs], px[:, -2])
        a_mat = np.linalg.cholesky(cov_mat)
    except np.linalg.LinAlgError as e:
        # When the covariance matrix is not positive definite, print the name
        # and range of each parameter, and the covariance matrix itself.
        names = [name for (ix, name) in p_info if ix in p_ixs]
        mins = np.min(px[:, p_ixs], axis=0)
        maxs = np.max(px[:, p_ixs], axis=0)
        means = np.mean(px[:, p_ixs], axis=0)
        mat_lines = str(cov_mat).splitlines()
        mat_sep = "\n      "
        mat_disp = mat_sep.join(["Covariance matrix:"] + mat_lines)
        logger = logging.getLogger(__name__)
        logger.warn("Post-RPF Cholesky decomposition: {}".format(e))
        logger.warn("Post-RPF parameters: {}".format(", ".join(names)))
        logger.warn("Minimum values: {}".format(mins))
        logger.warn("Maximum values: {}".format(maxs))
        logger.warn("Mean values:    {}".format(means))
        logger.warn(mat_disp)
        if ctx.params['resample']['regularise_or_fail']:
            raise
        else:
            return

    # Sample the multivariate normal with covariance V and mean of zero.
    std_samples = rnd.normal(size=(npar, count))
    scaled_samples = np.transpose(np.dot(a_mat, h * std_samples))

    #print(cov_mat)
    #print(a_mat)
    #print(std_samples[:, :5])
    #print(scaled_samples[:5])
    #print()

    # Add the sampled noise and clip to respect parameter bounds.
    new_px[:, p_ixs] = np.clip(
        new_px[:, p_ixs] + scaled_samples,
        ctx.params['model']['param_min'][None, p_ixs],
        ctx.params['model']['param_max'][None, p_ixs])


def resample(ctx, px):
    """Resample a particle population.

    :param params: The simulation parameters.
    :param px: An array of particle state vectors.

    The supported resampling methods are:

    - ``'basic'``:         uniform random numbers from [0, 1].
    - ``'stratified'``:    uniform random numbers from [j / m, (j + 1) / m).
    - ``'deterministic'``: select (j - a) / m for some fixed a.

    Where m is the number of particles and j = 0, ..., m - 1.

    These algorithms are described in G Kitagawa, J Comp Graph Stat
    5(1):1-25, 1996.
    `doi:10.2307/1390750 <https://doi.org/10.2307/1390750>`_
    """
    # Sort the particle indices according to weight (in descending order), so
    # that we can determine the original index of each resampled particle.
    # Use the merge sort algorithm because it is stable (thus preserving the
    # behaviour of Python's built-in `sorted` function).
    sorted_ix = np.argsort(- px[:, -2], kind='mergesort')
    # Sort the weights in descending order.
    sorted_ws = px[sorted_ix, -2]
    # Calculate the upper bounds for each interval.
    bounds = np.cumsum(sorted_ws)
    # Generate the random samples using the specified resampling method.
    count = px.shape[0]
    method = ctx.params['resample']['method']
    rnd = ctx.component['random']['resample']
    if method == 'basic':
        choices = np.sort(rnd.uniform(size=count))
    elif method == 'stratified':
        choices = (rnd.uniform(size=count) + np.arange(count)) / count
    elif method == 'deterministic':
        choices = (rnd.uniform() + np.arange(count)) / count
    else:
        # This is an error.
        raise ValueError("Invalid resampling method '{}'".format(method))
    # Resample the particles.
    new_px = np.copy(px)
    # Since the intervals and random samples are both monotonic increasing, we
    # only need step through the samples and record the current interval.
    bix = 0
    for (j, rand_val) in enumerate(choices):
        while bounds[bix] < rand_val:
            bix += 1
        new_px[j, 0:-2] = px[sorted_ix[bix]][0:-2]
        new_px[j, -1] = sorted_ix[bix]
    # Renormalise the weights.
    new_px[:, -2] = 1.0 / count
    # Sample model parameter values from a regularised kernel, if requested.
    if ctx.params['resample']['regularisation']:
        post_regularise(ctx, px, new_px)
    # Copy the resampled particles back into the original array.
    px[:, :] = new_px[:, :]


def resample_weights(weights, rnd, method='deterministic'):
    """
    Resample a particle weight array.

    :param np.ndarray weights: The particle weights.
    :param rnd: A random number generator.
    :param method: The resampling method: ``'basic'``, ``'stratified'``, or
        ``'deterministic'`` (default).
    :returns: A ``(sample_ixs, weight)`` tuple, where ``sample_ixs`` are the
        indices of the resampled particles and ``weight`` is the new weight
        for each particle (a single float).
    """
    # Resample to obtain a random sample from the correct density.
    params = {
        'resample': {
            'method': method,
            'regularisation': False,
        },
    }
    component = {
        'random': {
            'resample': rnd,
        },
    }

    ctx = context.Scaffold(params=params, component=component)

    px = np.array([weights, np.zeros(weights.shape)]).T
    resample(ctx, px)
    sample_ixs = px[:, -1].astype(int)
    new_weight = px[0, -2]
    return (sample_ixs, new_weight)
