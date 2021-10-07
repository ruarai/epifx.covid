.. _gs-example-obs:

Observations
============

For simplicity, we assume that both the prey and predator populations ---
:math:`x(t)` and :math:`y(t)` --- are directly observed, and that the
observation error is distributed normally with zero mean and a known standard
deviation :math:`\sigma_\mathrm{obs}`.

   .. math::

      \mathbf{x_t} &= [x, y, \alpha, \beta, \delta, \gamma]^T \\
      \mathcal{L}(x(t) \mid \mathbf{x_t}) &\sim \mathcal{N}(\mu=x,
          \sigma=\sigma_\mathrm{obs}) \\
      \mathcal{L}(y(t) \mid \mathbf{x_t}) &\sim \mathcal{N}(\mu=y,
          \sigma=\sigma_\mathrm{obs})

These observation models are implemented by
:class:`pypfilt.examples.predation.ObsModel`, and are used to update our
beliefs about how well each particle :math:`\mathbf{x_t}` explains each of the
provided observations (shown below).

.. https://sphinx-rtd-theme.readthedocs.io/en/stable/demo/lists_tables.html

.. hlist::
   :columns: 2

   - .. literalinclude:: /../src/pypfilt/examples/predation-counts-x.ssv
        :language: text
        :caption: Example observations of :math:`x(t)`.
        :name: obs-x

   - .. literalinclude:: /../src/pypfilt/examples/predation-counts-y.ssv
        :language: text
        :caption: Example observations of :math:`y(t)`.
        :name: obs-y
