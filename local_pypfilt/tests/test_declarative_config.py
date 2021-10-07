"""Test cases for the declarative configuration module."""

import pypfilt.config
import pkgutil


def test_flat_overrides():
    """Ensure that one flat dictionary correctly overrides another."""
    base = {'a': 1, 'b': 2, 'c': 3}
    update = {'a': 5, 'd': 7}
    pypfilt.config.override_dict(base, update)
    assert base == {'a': 5, 'b': 2, 'c': 3, 'd': 7}


def test_nested_overrides_1():
    """Ensure that one nested dictionary correctly overrides another."""
    base = {'a': 1, 'nested': {'b': 2, 'c': 3, 'd': {'x': 1}}}
    update = {'nested': {'b': 3, 'd': {'y': 5}}}
    pypfilt.config.override_dict(base, update)
    assert base == {'a': 1, 'nested': {'b': 3, 'c': 3, 'd': {'x': 1, 'y': 5}}}


def test_nested_overrides_2():
    """Ensure that one nested dictionary correctly overrides another."""
    base = {'a': 1, 'nested': {'b': 2}, 'replaced': 5}
    update = {'nested': {'b': 3}, 'added': {'f': 10}, 'replaced': {'x': 1}}
    pypfilt.config.override_dict(base, update)
    assert base == {'a': 1, 'nested': {'b': 3}, 'added': {'f': 10},
                    'replaced': {'x': 1}}


def test_load_example_config():
    """Ensure that the provide configuration can be parsed."""
    toml_data = pkgutil.get_data('pypfilt.examples', 'predation.toml')
    config = pypfilt.config.from_string(toml_data.decode())
    assert isinstance(config, pypfilt.config.Config)
