#!/usr/bin/env python
"""
An example of using the ``pypfilt`` package to estimate the state of a
two-species system described by the Lotka-Volterra equations.
"""

import pypfilt
import pypfilt.summary
import numpy as np
import scipy.integrate
import scipy.stats
import h5py
import pypfilt.plot
import logging
import sys
import pkgutil


class LotkaVolterra(pypfilt.Model):
    """An implementation of the (continuous) Lotka-Volterra equations."""

    def init(self, ctx, vec):
        """Initialise a matrix of state vectors."""
        # Select x(0), y(0), and the parameters according to the priors.
        rnd = ctx.component['random']['model']
        size = vec[..., 0].shape
        vec[..., 0] = ctx.params['model']['prior']['x'](rnd, size)
        vec[..., 1] = ctx.params['model']['prior']['y'](rnd, size)
        vec[..., 2] = ctx.params['model']['prior']['alpha'](rnd, size)
        vec[..., 3] = ctx.params['model']['prior']['beta'](rnd, size)
        vec[..., 4] = ctx.params['model']['prior']['gamma'](rnd, size)
        vec[..., 5] = ctx.params['model']['prior']['delta'](rnd, size)

    def state_size(self):
        """Return the size of the state vector."""
        return 6

    def d_dt(self, xt, t):
        """Calculate the derivatives of x(t) and y(t)."""
        # Restore the 2D shape of the flattened state matrix.
        xt = xt.reshape((-1, 6))
        x, y = xt[..., 0], xt[..., 1]
        d_dt = np.zeros(xt.shape)
        # Calculate dx/dt and dy/dt.
        d_dt[..., 0] = xt[..., 2] * x - xt[..., 3] * x * y
        d_dt[..., 1] = xt[..., 4] * x * y - xt[..., 5] * y
        # Flatten the 2D derivatives matrix.
        return d_dt.reshape(-1)

    def update(self, ctx, t, dt, is_fs, prev, curr):
        """Perform a single time-step."""
        # Use scalar time, so that ``t + dt`` is well-defined.
        t = ctx.component['time'].to_scalar(t)
        # The state matrix must be flattened for odeint.
        xt = scipy.integrate.odeint(self.d_dt, prev.reshape(-1),
                                    [t, t + dt])[1]
        # Restore the 2D shape of the flattened state matrix.
        curr[:] = xt.reshape(curr.shape)

    def describe(self):
        """Describe each component of the state vector."""
        return [
            # Restrict x(t), y(t) to [0, 10^5], don't allow regularisation.
            ('x', False, 0, 1e5),
            ('y', False, 0, 1e5),
            # Restrict parameters to [0, 2], allow regularisation.
            ('alpha', True, 0, 2),
            ('beta', True, 0, 2),
            ('gamma', True, 0, 2),
            ('delta', True, 0, 2),
        ]


class ObsModel(pypfilt.Obs):
    def __init__(self, obs_unit, obs_period):
        self.unit = obs_unit
        self.period = obs_period

    def log_llhd(self, params, op, time, obs, curr, hist):
        # NOTE: the expected observations are x(t) and y(t).
        # Calculate the log-likelihood of each observation in turn.
        unit = obs['unit']
        if unit == 'x':
            x_t = curr[..., 0]
            x_dist = scipy.stats.norm(loc=x_t, scale=op['sdev'])
            return x_dist.logpdf(obs['value'])
        elif unit == 'y':
            y_t = curr[..., 1]
            y_dist = scipy.stats.norm(loc=y_t, scale=op['sdev'])
            return y_dist.logpdf(obs['value'])
        else:
            raise ValueError('invalid observation unit: {}'.format(unit))

    def simulate(self, params, op, time, period, expect, rng=None):
        if rng is None:
            return scipy.stats.norm(loc=expect, scale=op['sdev']).rvs()
        else:
            return rng.normal(loc=expect, scale=op['sdev'])

    def expect(self, ctx, op, time, period, prev, curr):
        if self.unit == 'x':
            expect = curr[..., 0]
        elif self.unit == 'y':
            expect = curr[..., 1]
        else:
            raise ValueError('invalid observation unit: {}'.format(self.unit))

        return expect

    def quantiles(self, params, op, time, mu, wt, probs):
        # The minimum interval width before we decide that a value is
        # sufficiently accurate.
        tolerance = 0.00001
        scale = op['sdev']
        normal = scipy.stats.norm(loc=mu, scale=scale)

        def cdf(y):
            """Calculate the CDF of the weighted sum over all particles."""
            return np.dot(wt, normal.cdf(y))

        def bisect(a, b):
            """
            Return the midpoint of the interval [a, b], or ``None`` if the
            minimum tolerance has been reached.
            """
            if b > a + tolerance:
                return (a + b) / 2
            else:
                return None

        # Find appropriate lower and upper bounds for y_i.
        pr_min = np.min(probs)
        pr_max = np.max(probs)
        y0_lower = scipy.stats.norm(loc=np.min(mu), scale=scale).ppf(pr_min)
        y0_upper = scipy.stats.norm(loc=np.max(mu), scale=scale).ppf(pr_max)

        return pypfilt.obs.bisect_cdf(probs, cdf, bisect, y0_lower, y0_upper)

    def from_file(self, filename, time_scale):
        cols = [time_scale.column('date'), ('value', np.float)]
        df = pypfilt.io.read_table(filename, cols)
        obs_list = [{'date': row['date'],
                     'value': row['value'],
                     'unit': self.unit,
                     'period': self.period,
                     'source': filename}
                    for row in df]
        return (obs_list, df)


def default_priors():
    """Define default model prior distributions."""
    return {
        'x': lambda r, size=None: r.uniform(0.5, 1.5, size=size),
        'y': lambda r, size=None: r.uniform(0.2, 0.4, size=size),
        'alpha': lambda r, size=None: r.uniform(0.6, 0.8, size=size),
        'beta': lambda r, size=None: r.uniform(1.2, 1.4, size=size),
        'gamma': lambda r, size=None: r.uniform(0.9, 1.1, size=size),
        'delta': lambda r, size=None: r.uniform(0.9, 1.1, size=size),
    }


def make_params(px_count, seed, obs_sdev, max_days=14):
    """Define the default simulation parameters for this model."""
    model = LotkaVolterra()
    time_scale = pypfilt.Scalar()
    params = pypfilt.default_params(model, time_scale,
                                    max_days=max_days,
                                    px_count=px_count,
                                    prng_seed=seed)
    # Use one time-step per unit time, odeint will interpolate as needed.
    params['steps_per_unit'] = 1
    # Calculate statistics from the start of the simulation period.
    params['summary']['from_first_day'] = True
    # Define default model prior distributions.
    params['model']['prior'] = default_priors()

    # Define the observation model parameters and likelihood functions.
    params['obs'] = {
        'x': {'sdev': obs_sdev},
        'y': {'sdev': obs_sdev},
    }
    params['component']['obs'] = {
        'x': ObsModel(obs_unit='x', obs_period=0),
        'y': ObsModel(obs_unit='y', obs_period=0),
    }
    # Write output to the working directory.
    params['out_dir'] = '.'
    params['tmp_dir'] = '.'
    return params


def make_observations(params, obs_tables=True):
    # Record the original prior distributions and particle count.
    original_priors = params['model']['prior']
    px_count = params['hist']['px_count']

    # Define the ground truth and construct the corresponding priors.
    x0 = 0.9
    y0 = 0.25
    alpha = 2/3
    beta = 4/3
    gamma = 1
    delta = 1
    params['model']['prior'] = {
        'x': lambda r, size=None: x0 * np.ones(size),
        'y': lambda r, size=None: y0 * np.ones(size),
        'alpha': lambda r, size=None: alpha * np.ones(size),
        'beta': lambda r, size=None: beta * np.ones(size),
        'gamma': lambda r, size=None: gamma * np.ones(size),
        'delta': lambda r, size=None: delta * np.ones(size),
    }

    # Simulate the observations from this model.
    params['hist']['px_count'] = 1
    sim_obs = pypfilt.simulate_from_model(params)

    # Restore the original prior distributions and particle count.
    params['model']['prior'] = original_priors
    params['hist']['px_count'] = px_count

    # Convert each row in the simulated observations table into an observation
    # dictionary. Note that this involves converting the observation dates
    # from their serialised form.
    time = params['component']['time']
    obs = []
    for row in sim_obs:
        obs.append({
            'date': time.from_dtype(row['date']),
            'period': 0,
            'unit': row['unit'],
            'value': row['value'],
            'source': 'make_observations()',
        })

    if obs_tables:
        params['data']['obs']['x'] = sim_obs[sim_obs['unit'] == 'x']
        params['data']['obs']['x'] = sim_obs[sim_obs['unit'] == 'y']

    return obs


def save_scalar_observations(sim_obs):
    """Save simulated observations to disk."""
    x_tbl = sim_obs[sim_obs['unit'] == 'x'][['date', 'value']]
    y_tbl = sim_obs[sim_obs['unit'] == 'y'][['date', 'value']]
    x_tbl = x_tbl[x_tbl['date'] > 0]
    y_tbl = y_tbl[y_tbl['date'] > 0]
    np.savetxt('predation-counts-x.ssv', x_tbl, fmt='%d %f',
               header='date value', comments='')
    np.savetxt('predation-counts-y.ssv', y_tbl, fmt='%d %f',
               header='date value', comments='')


def forecast(data_file):
    """Run a suite of forecasts against generated observations."""
    logger = logging.getLogger(__name__)
    logger.info('Preparing the forecast simulations')

    # Define the simulation period and forecasting times.
    t0 = 0.0
    t1 = 15.0
    fs_times = [1.0, 3.0, 5.0, 7.0, 9.0]
    params = make_params(px_count=1000, seed=42, obs_sdev=0.2)
    params['time']['start'] = t0
    params['time']['until'] = t1

    # Generate noisy observations.
    obs = make_observations(params, t1)

    # Define the summary tables to be saved to disk.
    summary = pypfilt.summary.HDF5(params, obs)
    params['component']['summary'] = summary
    params['component']['summary_monitor'] = {
        'expected_obs': pypfilt.summary.ExpectedObsMonitor(),
    }
    params['component']['summary_table'] = {
        'model_cints': pypfilt.summary.ModelCIs(probs=[0, 50, 95]),
        'obs': pypfilt.summary.Obs(),
        'forecasts': pypfilt.summary.PredictiveCIs('expected_obs'),
    }

    # Run the forecast simulations.
    pypfilt.forecast(params, [obs], fs_times, data_file)


def plot_forecasts(state_cints, x_obs, y_obs, pdf_file=None, png_file=None):
    """Plot the population predictions at each forecasting date."""
    logger = logging.getLogger(__name__)
    with pypfilt.plot.apply_style():
        plot = pypfilt.plot.Grid(
            state_cints, 'Time', 'Population Size (1,000s)',
            ('fs_date', 'Forecast @ t = {:0.0f}'),
            ('unit', lambda s: '{}(t)'.format(s)))
        plot.expand_x_lims('date')
        plot.expand_y_lims('ymax')

        for (ax, df) in plot.subplots():
            ax.axhline(y=0, xmin=0, xmax=1,
                       linewidth=1, linestyle='--', color='k')
            hs = pypfilt.plot.cred_ints(ax, df, 'date', 'prob')
            if df['unit'][0] == 'x':
                df_obs = x_obs
            else:
                df_obs = y_obs
            past_obs = df_obs[df_obs['date'] <= df['fs_date'][0]]
            future_obs = df_obs[df_obs['date'] > df['fs_date'][0]]
            hs.extend(pypfilt.plot.observations(ax, past_obs,
                                                label='Past observations'))
            hs.extend(pypfilt.plot.observations(ax, future_obs,
                                                future=True,
                                                label='Future observations'))
            plot.add_to_legend(hs)

            # Adjust the axis limits and the number of ticks.
            ax.set_xlim(left=0)
            ax.locator_params(axis='x', nbins=4)
            ax.set_ylim(bottom=-0.2)
            ax.locator_params(axis='y', nbins=4)

        plot.legend(loc='upper center', ncol=5)

        if pdf_file:
            logger.info('Plotting to {}'.format(pdf_file))
            plot.save(pdf_file, format='pdf', width=10, height=5)
        if png_file:
            logger.info('Plotting to {}'.format(png_file))
            plot.save(png_file, format='png', width=10, height=5)


def plot_params(param_cints, pdf_file=None, png_file=None):
    """Plot the parameter posteriors over the estimation run."""
    logger = logging.getLogger(__name__)
    with pypfilt.plot.apply_style():
        plot = pypfilt.plot.Wrap(
            param_cints, 'Time', 'Value',
            ('name', lambda s: '$\\{}$'.format(s)),
            nr=1)
        plot.expand_y_lims('ymax')

        for (ax, df) in plot.subplots(dy=-0.025):
            hs = pypfilt.plot.cred_ints(ax, df, 'date', 'prob')
            if df['name'][0] == 'alpha':
                y_true = 2/3
            elif df['name'][0] == 'beta':
                y_true = 4/3
            elif df['name'][0] == 'gamma':
                y_true = 1
            elif df['name'][0] == 'delta':
                y_true = 1
            hs.append(ax.axhline(y=y_true, xmin=0, xmax=1, label='True value',
                                 linewidth=1, linestyle='--', color='k'))
            plot.add_to_legend(hs)

        plot.legend(loc='upper center', ncol=5, borderaxespad=0)

        if pdf_file:
            logger.info('Plotting to {}'.format(pdf_file))
            plot.save(pdf_file, format='pdf', width=10, height=3)
        if png_file:
            logger.info('Plotting to {}'.format(png_file))
            plot.save(png_file, format='png', width=10, height=3)


def plot(data_file, png=True, pdf=True):
    """
    Save the plots produced by :func:`plot_params` and :func:`plot_forecasts`.

    This will save the plots to files whose names begin with
    "predation_params" and "predation_forecasts".

    :param png: Whether to save plots as PNG files.
    :param pdf: Whether to save plots as PDF files.
    """
    logger = logging.getLogger(__name__)
    logger.info('Loading outputs from {}'.format(data_file))

    # Use the 'Agg' backend so that plots can be generated non-interactively.
    import matplotlib
    matplotlib.use('Agg')

    # File names for the generated plots.
    fs_pdf = 'predation_forecasts.pdf'
    fs_png = 'predation_forecasts.png'
    pp_pdf = 'predation_params.pdf'
    pp_png = 'predation_params.png'

    # Read in the model credible intervals and the observations.
    with h5py.File(data_file, 'r') as f:
        cints = f['/data/model_cints'][()]
        forecasts = f['/data/forecasts'][()]
        obs = f['/data/obs'][()]

    # Convert serialised values into more convenient data types.
    convs = pypfilt.summary.default_converters(pypfilt.Scalar())
    cints = pypfilt.summary.convert_cols(cints, convs)
    forecasts = pypfilt.summary.convert_cols(forecasts, convs)
    obs = pypfilt.summary.convert_cols(obs, convs)

    # Separate the observations of the two populations.
    x_obs = obs[obs['unit'] == 'x']
    y_obs = obs[obs['unit'] == 'y']

    # Separate the credible intervals for the population sizes from the
    # credible intervals for the model parameters.
    var_mask = np.logical_or(cints['name'] == 'x',
                             cints['name'] == 'y')
    param_cints = cints[np.logical_not(var_mask)]

    # Only retain forecasts, ignore results from the estimation pass, if any.
    fs_mask = forecasts['fs_date'] < max(forecasts['date'])
    forecasts = forecasts[fs_mask]

    # Only keep the model parameter posteriors from the estimation run.
    est_mask = param_cints['fs_date'] == max(param_cints['date'])
    param_cints = param_cints[est_mask]

    # Plot the population forecasts.
    pdf_file = fs_pdf if pdf else None
    png_file = fs_png if png else None
    plot_forecasts(forecasts, x_obs, y_obs, pdf_file, png_file)

    # Plot the model parameter posterior distributions.
    pdf_file = pp_pdf if pdf else None
    png_file = pp_png if png else None
    plot_params(param_cints, pdf_file, png_file)


def __example_data(filename):
    return pkgutil.get_data('pypfilt.examples', filename).decode()


def example_toml_data():
    """Return the contents of the example file "predation.toml"."""
    return __example_data('predation.toml')


def example_obs_x_data():
    """Return the contents of the example file "predation-counts-x.ssv"."""
    return __example_data('predation-counts-x.ssv')


def example_obs_y_data():
    """Return the contents of the example file "predation-counts-y.ssv"."""
    return __example_data('predation-counts-y.ssv')


def write_example_files():
    """
    Save the following example files to the working directory:

    * The forecast scenario file "predation.toml";
    * The observations file "predation-counts-x.ssv"; and
    * The observations file "predation-counts-y.ssv".
    """
    toml_file = 'predation.toml'
    obs_x_file = 'predation-counts-x.ssv'
    obs_y_file = 'predation-counts-y.ssv'

    toml_data = example_toml_data()
    with open(toml_file, 'w') as f:
        f.write(toml_data)

    obs_x_data = example_obs_x_data()
    with open(obs_x_file, 'w') as f:
        f.write(obs_x_data)

    obs_y_data = example_obs_y_data()
    with open(obs_y_file, 'w') as f:
        f.write(obs_y_data)


def main(args=None):
    logging.basicConfig(level=logging.INFO)
    data_file = 'predation.hdf5'
    forecast(data_file)
    plot(data_file, pdf=False)


if __name__ == '__main__':
    sys.exit(main())
