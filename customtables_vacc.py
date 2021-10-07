import numpy as np
import pypfilt.resample
from pypfilt.summary import Table


class SEEIIR(Table):
    """
    Record the state (S, E1, E2, I1, I2, R) of each particle at each day.
    """
    def dtype(self, ctx, obs_list, name):
        self.__ctx = ctx
        self.__params = ctx.params
        self.__popn_size = ctx.params['model']['population_size']
        self.__rnd = np.random.default_rng(ctx.params.get('prng_seed'))
        self.__resample = ctx.params.get_chained(
            ['summary', 'tables', name, 'resample'],
            default=True)
        self.__sample_ixs = None
        fs_date = ctx.component['time'].dtype('fs_date')
        date = ctx.component['time'].dtype('date')
        ix = ('ix', np.int32)
        wt = ('weight', np.float_)
        s_u = ('S_U', np.int32)
        e1_u = ('E1_U', np.int32)
        e2_u = ('E2_U', np.int32)
        i1_u = ('I1_U', np.int32)
        i2_u = ('I2_U', np.int32)
        r_u = ('R_U', np.int32)

        s_v = ('S_V', np.int32)
        e1_v = ('E1_V', np.int32)
        e2_v = ('E2_V', np.int32)
        i1_v = ('I1_V', np.int32)
        i2_v = ('I2_V', np.int32)
        r_v = ('R_V', np.int32)

        sigma = ('sigma', np.float64)
        gamma = ('gamma', np.float64)
        reff = ('R_eff', np.float64)

        adjustment = ('adjustment', np.float64)
        mean_Ei = ('mean_Ei', np.float64)
        mean_Et = ('mean_Et', np.float64)

        return [fs_date, date, ix, wt,
                s_u, e1_u, e2_u, i1_u, i2_u, r_u,
                s_v, e1_v, e2_v, i1_v, i2_v, r_v,
                reff, sigma, gamma,
                adjustment,
                mean_Ei, mean_Et]

    def n_rows(self, start_date, end_date, n_days, n_sys, forecasting):
        # Need one row for each particle.
        self.__sample_ixs = None
        num_px = self.__params['hist']['px_count']
        if not self.__resample:
            self.__sample_ixs = np.array(range(num_px))
        return num_px * n_days

    def add_rows(self, hist, weights, fs_date, dates, obs_types, insert_fn):
        if self.__sample_ixs is None:
            # Resample the particles once, before adding any rows for each
            # forecast, so that the particle weights are uniform.
            weight_ix = dates[0][1]
            weights_in = weights[weight_ix, :]
            (sample_ixs, _weight) = pypfilt.resample.resample_weights(
                weights_in, self.__rnd)
            self.__sample_ixs = sample_ixs

        fs_date_enc = self.__ctx.component['time'].to_dtype(fs_date)
        for date, date_ix, hist_ix in dates:
            date_enc = self.__ctx.component['time'].to_dtype(date)
            for (sample_ix, px) in enumerate(self.__sample_ixs):
                insert_fn((fs_date_enc, date_enc, sample_ix,
                           weights[date_ix, px],    # weight

                           hist[hist_ix, px, 0],    # Su(t)
                           hist[hist_ix, px, 1],    # E1u(t)
                           hist[hist_ix, px, 2],    # E2u(t)
                           hist[hist_ix, px, 3],    # I1u(t)
                           hist[hist_ix, px, 4],    # I2u(t)
                           hist[hist_ix, px, 5],    # Ru(t)

                           hist[hist_ix, px, 6],    # Sv(t)
                           hist[hist_ix, px, 7],    # E1v(t)
                           hist[hist_ix, px, 8],    # E2v(t)
                           hist[hist_ix, px, 9],    # I1v(t)
                           hist[hist_ix, px, 10],    # I2v(t)
                           hist[hist_ix, px, 11],    # Rv(t)

                           hist[hist_ix, px, 17],    # R_eff(t)
                           hist[hist_ix, px, 13],    # sigma
                           hist[hist_ix, px, 14],    # gamma
                           hist[hist_ix, px, 18],    # adjustment
                           hist[hist_ix, px, 19],    # mean_Ei
                           hist[hist_ix, px, 20]))   # mean_Et


class ExposureAndPresymptomaticEnsemble(Table):
    """
    Record the cumulative number of exposures and the prevalence of
    pre-symptomatic individuals for each particle.

    Note that this currently only supports the ``epifx.stoch.SEEIIR`` model.
    """

    def dtype(self, ctx, obs_list, name):
        self.__ctx = ctx
        self.__params = ctx.params
        self.__popn_size = ctx.params['model']['population_size']
        self.__rnd = np.random.default_rng(ctx.params.get('prng_seed'))
        self.__sample_ixs = None
        fs_date = ctx.component['time'].dtype('fs_date')
        date = ctx.component['time'].dtype('date')
        ix = ('ix', np.int32)
        cum_exposures = ('cum_exposures', np.int32)
        pre_symptomatic = ('pre_symptomatic', np.int32)
        return [fs_date, date, ix, cum_exposures, pre_symptomatic]

    def n_rows(self, start_date, end_date, n_days, n_sys, forecasting):
        # Need one row for each particle.
        if not forecasting:
            return 0
        self.__sample_ixs = None
        num_px = self.__params['hist']['px_count']
        return num_px * n_days

    def add_rows(self, hist, weights, fs_date, dates, obs_types, insert_fn):
        if self.__sample_ixs is None:
            # Resample the particles once, before adding any rows for each
            # forecast, so that the particle weights are uniform.
            weight_ix = dates[0][1]
            weights_in = weights[weight_ix, :]
            (sample_ixs, _weight) = pypfilt.resample.resample_weights(
                weights_in, self.__rnd)
            self.__sample_ixs = sample_ixs

        fs_date_enc = self.__ctx.component['time'].to_dtype(fs_date)
        for date, _ix, hist_ix in dates:
            date_enc = self.__ctx.component['time'].to_dtype(date)
            for (sample_ix, px) in enumerate(self.__sample_ixs):
                # The cumulative number of exposures is N - S(t).
                cum_exps = self.__popn_size - hist[hist_ix, px, 0]
                # Pre-symptomatic prevalence is E1(t) + E2(t) + I1(t).
                pre_symp_prev = np.sum(hist[hist_ix, px, 1:4])
                insert_fn((fs_date_enc, date_enc, sample_ix,
                           cum_exps, pre_symp_prev))
