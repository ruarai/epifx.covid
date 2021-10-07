import copy
import numpy as np
import types


class Context:
    """
    A simulation context, which encapsulates the simulation parameters,
    simulation components, data tables, event handlers, etc.

    :param params: The simulation parameters.

    :Examples:

    >>> import pypfilt
    >>> import pypfilt.examples.predation
    >>> from pypfilt.context import Context
    >>> pypfilt.examples.predation.write_example_files()
    >>> config_file = 'predation.toml'
    >>> for forecast in pypfilt.forecasts_iter(config_file):
    ...     ctx = Context(forecast.params)
    ...     # Inspect parameters.
    ...     print(ctx.params['prng_seed'])
    ...     # Inspect components.
    ...     print(ctx.component['time'].__class__)
    ...     # Inspect data tables.
    ...     print(ctx.data['obs']['x'].shape)
    ...     # Inspect event handlers.
    ...     print(ctx.hook)
    42
    <class 'pypfilt.time.Scalar'>
    (15,)
    {'log_llhd': []}
    """
    def __init__(self, params):
        model = params['component']['model']
        # Ensure that prior distributions have been defined.
        if 'prior' not in params['model'] or not params['model']['prior']:
            raise ValueError('Model prior distribution(s) are undefined')

        params['epoch'] = params['time']['start']
        params['dt'] = 1.0 / params['steps_per_unit']

        params['size'] = params['hist']['px_count']
        # Record the number of state columns and extra columns.
        params['hist']['state_cols'] = model.state_size()
        # Allocate an array that enumerates the particles.
        params['px_range'] = np.arange(params['hist']['px_count'])

        # NOTE: we copy each of the lists/dictionaries so that any
        # modifications made during the simulation have no effect on the
        # input parameters dictionary.

        # Simulation components.
        component = params['component']
        self.component = component.copy()

        # Event hooks.
        hooks = params['hooks']
        self.hook = hooks.copy()

        # Data tables.
        data = params['data']
        self.data = data.copy()

        # Construct the random number generators.
        for name in params['random']:
            self.component['random'][name] = np.random.default_rng(
                params['prng_seed'])

        # Simulation parameters.
        del params['component']
        del params['hooks']
        del params['data']
        self.params = Params(params)
        params['component'] = component
        params['hooks'] = hooks
        params['data'] = data

    def install_hook(self, hook_name, hook_fn):
        """
        Register a function that should be called in response to an event.

        :param hook_name: The event name.
        :param hook_fn: The event-handler function.
        """
        if hook_name not in self.hook:
            self.hook[hook_name] = [hook_fn]
        elif hook_fn not in self.hook[hook_name]:
            self.hook[hook_name].append(hook_fn)

    def call_hooks(self, hook_name, *args, **kwargs):
        """
        Call all event-handler functions associated with an event name.

        :param hook_name: The event name.
        :param \\*args: Positional arguments to pass to the event handlers.
        :param \\**kwargs: Keyword arguments to pass to the event handlers.
        """
        if hook_name in self.hook:
            for hook_fn in self.hook[hook_name]:
                hook_fn(*args, **kwargs)


def into_nonmut(value):
    # Convert known mutable types into non-mutable equivalents.
    if isinstance(value, dict):
        return Params(value)
    elif isinstance(value, list):
        # NOTE: recursively descend into the list items.
        return tuple(into_nonmut(item) for item in value)
    elif isinstance(value, set):
        # NOTE: recursively descend into the set items.
        return frozenset(into_nonmut(item) for item in value)
    elif isinstance(value, (types.FunctionType, types.LambdaType)):
        # NOTE: no need to make a copy of functions.
        return value
    elif isinstance(value, (float, int, bool, str, tuple, type(None))):
        # NOTE: no need to make copies of basic types.
        return value
    elif isinstance(value, np.ndarray):
        # NOTE: NumPy arrays provide their own copy() implementation.
        return value.copy()
    else:
        # Make a deep copy of all other values.
        # print('Making deep copy of {}'.format(type(value)))
        return copy.deepcopy(value)


class Params:
    """
    A non-mutable dictionary wrapper that ensures the simulation parameters
    remain unchanged during each simulation.
    """
    def __init__(self, params):
        self.__params = {}
        for (key, value) in params.items():
            self.__params[key] = into_nonmut(value)

    def __len__(self):
        return len(self.__params)

    def __getitem__(self, key):
        return self.__params[key]

    def __setitem__(self, key, value):
        raise ValueError('Cannot set parameter "{}"'.format(key))

    def __delitem__(self, key):
        raise ValueError('Cannot delete parameter "{}"'.format(key))

    def __missing__(self, key):
        return self.__params.__missing__(key)

    def __iter__(self):
        return self.__params.__iter__()

    def __contains__(self, item):
        return item in self.__params

    def __str__(self):
        return str(self.__params)

    def __repr__(self):
        return self.__class__.__name__ + "(" + repr(self.__params) + ")"

    def copy(self):
        return Params(self.__params.copy())

    def items(self):
        return self.__params.items()

    def keys(self):
        return self.__params.keys()

    def values(self):
        return self.__params.values()

    def get(self, key, default=None):
        """
        Return the value associated with key, if present; otherwise return
        ``default``.
        """
        return self.__params.get(key, default)

    def get_chained(self, keys, default=None):
        """
        Return the value associated with a sequence of keys, if all keys are
        present, otherwise return ``default``.

        :Examples:

        >>> from pypfilt.context import Params
        >>> x = {'a': 1, 'b': {'c': 2}}
        >>> params = Params(x)
        >>> params.get_chained(['b', 'c'])
        2
        >>> params.get_chained(['b', 'd']) is None
        True
        """
        result = self
        for key in keys:
            # NOTE: if we encounter a non-dictionary type, return the default.
            if not isinstance(result, Params):
                return default
            try:
                result = result[key]
            except KeyError:
                return default
        return result


class Scaffold:
    """
    Convenience class for defining minimal context objects.

    This is useful when, e.g., writing test cases for functions that require a
    context object.

    :Examples:

    >>> import numpy as np
    >>> import pypfilt.resample
    >>> from pypfilt.context import Scaffold
    >>> params = {'resample': {'method': 'basic', 'regularisation': False}}
    >>> component = {'random': {'resample': np.random.default_rng()}}
    >>> ctx = Scaffold(params=params, component=component)
    >>> weights = np.array([0.50, 0.25, 0.1, 0.1, 0.02, 0.02, 0.01])
    >>> ixs = np.zeros(weights.shape)
    >>> attempts = 10
    >>> for i in range(attempts):
    ...     x = np.array([weights, ixs]).T
    ...     pypfilt.resample.resample(ctx, x)
    ...     if any(x[:, 1] == 0):
    ...         # The first particle (weight 0.5) was selected at least once.
    ...         break
    """
    def __init__(self, **entries):
        self.__dict__.update(entries)
