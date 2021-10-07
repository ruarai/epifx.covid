.. _gs-equations:

Lotka-Volterra equations
========================

Here we will show how to generate forecasts for the (continuous)
`Lotka-Volterra equations`_, which describe the dynamics of biological systems
in which two species interact (one predator, one prey).

.. math::

   \frac{dx}{dt} &= \alpha x - \beta xy \\
   \frac{dy}{dt} &= \delta xy - \gamma y

==============  =========================================================
Symbol          Meaning
==============  =========================================================
:math:`x(t)`    The size of the prey population (1,000s).
:math:`y(t)`    The size of the predator population (1,000s).
:math:`\alpha`  Exponential prey growth rate in the absence of predators.
:math:`\beta`   The rate at which prey suffer from predation.
:math:`\delta`  The predator growth rate, driven by predation.
:math:`\gamma`  Exponential decay rate of the predator population.
==============  =========================================================

All of the state variables and parameters are stored in the state vector:

.. math::

   \mathbf{x_t} = [x, y, \alpha, \beta, \delta, \gamma]^T

.. _Lotka-Volterra equations: https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations
