try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

_PARENT_MODULE = ".".join(__name__.split('.')[:-1])

from . import vep


def read_vep_basic_schema():
    import yaml

    with pkg_resources.open_binary(vep.__name__, 'basic_schema.yaml') as fp:
        vep_basic_schema = yaml.safe_load(fp)

    return vep_basic_schema


def read_vep_basic_args():
    import yaml

    with pkg_resources.open_binary(vep.__name__, 'basic_args.yaml') as fp:
        vep_basic_args = yaml.safe_load(fp)

    return vep_basic_args


def read_vep_predefined_custom_annot():
    import yaml

    with pkg_resources.open_binary(vep.__name__, 'predefined_custom_annot.yaml') as fp:
        vep_predefined_custom_annot = yaml.safe_load(fp)

    return vep_predefined_custom_annot


_PKG_RESOURCES = {
    "vep/basic_schema": read_vep_basic_schema,
    "vep/basic_args": read_vep_basic_args,
    "vep/predefined_custom_annot": read_vep_predefined_custom_annot,
}


def get(uri):
    fn = _PKG_RESOURCES.get(uri)

    if not fn:
        raise ValueError("Unknown URI: %s" % uri)

    return fn()
