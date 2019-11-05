import importlib


def load_perturbation_function(version: str):
    PACKAGE = "srf_generation.source_parameter_generation.uncertainties.versions"
    foo = importlib.import_module("." + version, PACKAGE)
    return foo.generate_source_params
