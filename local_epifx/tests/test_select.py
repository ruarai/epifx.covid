import logging
import os
import warnings
import numpy as np
import pypfilt
import pypfilt.sweep
import epifx
import epifx.select


def test_select():
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("epifx.select")
    logger.setLevel(logging.DEBUG)
    out_file = 'test_select_samples.ssv'

    # Define the PRNG seed for the selection process.
    seed = 2020

    # Draw proposals from the prior distributions defined in the simulation
    # parameters.
    proposal = epifx.select.DefaultProposal()

    # Peak times and sizes based on weekly numbers of seasonal influenza case
    # notifications for metropolitan Melbourne over 2012-2017.
    peak_times = np.array([[136, 177, 176, 182, 187, 193]])
    peak_sizes = np.array([360, 417, 691, 1329, 975, 2036])
    target = epifx.select.TargetPeakMVN(peak_sizes, peak_times)

    # Write an empty observations file, otherwise an exception will be raised
    # when attempting to load the observations.
    obs_file = 'no-observations.ssv'
    with open(obs_file, 'w') as f:
        f.write('date count\n')

    config = pypfilt.config.from_string(config_str())
    with warnings.catch_warnings():
        # NOTE: pypfilt.sweep() will read the observations date file, which
        # makes numpy.loadtxt() produce a warning about an empty data file.
        # We can suppress this with warnings.filterwarnings().
        warnings.filterwarnings('ignore',
                                message='loadtxt: Empty input file:',
                                category=UserWarning)
        forecasts = list(pypfilt.sweep.forecasts(config))
    params = forecasts[0].params

    # Select particles according to the peak size and time target.
    vec = epifx.select.select(params, proposal, target, seed)

    # Retrieve the parameter columns from these particles.
    sample_cols = params['component']['model'].sample_columns()
    column_names = list(sample_cols.keys())
    column_ixs = np.array([sample_cols[n] for n in column_names])
    tbl = vec[:, column_ixs]

    # Save the sampled parameter values.
    logger.debug("Saving samples to {}".format(out_file))
    np.savetxt(out_file, tbl, header=' '.join(column_names), comments='')

    # Remove the empty observations file.
    os.remove(obs_file)
    # Leave the sampled parameter file for now.
    # os.remove(out_file)


def config_str():
    """Define forecast scenarios for these test cases."""
    return """
    [components]
    model = "epifx.det.SEIR"
    time = "pypfilt.Datetime"
    summary = "epifx.summary.make"

    [parameters]
    particles = 2000
    prng_seed = 3001
    steps_per_unit = 1
    last_n_periods = 1
    data_dir = "."
    tmp_dir = "."
    out_dir = "."
    json_dir = "."
    max_days = 7
    resample.threshold = 0.25
    resample.regularisation = true
    time.start = "2020-03-01"
    time.until = "2020-12-31"
    summary.from_first_day = true
    summary.only_forecasts = true
    summary.metadata.packages = [ "epifx" ]
    fresh_cache = true
    remove_cache = true

    [model.bounds]
    R0 = { min = 1.2, max = 1.6 }
    sigma = { min = 0.1, max = 10.0 }
    gamma = { min = 0.1, max = 10.0 }
    eta = { min = 1.0, max = 1.0 }
    alpha = { min = -0.2, max = 0.0 }
    t0 = { min = 0, max = 60 }

    [model.priors]
    R0 = { function = "uniform", args.low = 1.2, args.high = 1.6 }
    sigma = { function = "inverse_uniform", args.low = 0.5, args.high = 4.0 }
    gamma = { function = "inverse_uniform", args.low = 0.5, args.high = 4.0 }
    eta = { function = "uniform", args.low = 1.0, args.high = 1.0 }
    alpha = { function = "uniform", args.low = 0.0, args.high = 0.0 }
    t0 = { function = "uniform", args.low = 0, args.high = 60 }

    [summary.init]
    default = true

    [scenario.Melbourne]
    name = "Melbourne"
    parameters.model.population_size = 5_191_000

    [observations.Notifications]
    model = "epifx.obs.PopnCounts"
    file = "no-observations.ssv"
    init.obs_period = 7
    file_args.time_col = "date"
    file_args.value_col = "count"
    parameters.bg_obs = 70
    parameters.bg_var = 80
    parameters.pr_obs = 0.00275
    parameters.disp = 100
    format = { bg_obs = "3.0f", bg_var = "03.0f", pr_obs = "0.5f", disp = "03.0f" }
    name = { bg_obs = "bg", bg_var = "bgvar", pr_obs = "pr", disp = "disp" }
    """
