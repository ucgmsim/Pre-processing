import importlib
import os


def __load_all__():
    PACKAGE = "SrfGen.source_parameter_generation.uncertainties.versions"
    perturbators = {}
    versions = {}
    list_modules = [
        f.replace(".py", "")
        for f in os.listdir(os.path.dirname(__file__))
        if f.endswith(".py")
    ]
    list_modules.remove("__init__")
    for module_name in list_modules:
        print("Load module ' ", module_name, "' :")
        foo = importlib.import_module("." + module_name, PACKAGE)
        if foo.VERSION in perturbators.keys():
            raise ImportError(
                "Version {} already exists in the versions module. Files {} and {} have the same version. Please change one of the files.".format(
                    foo.VERSION, foo.__file__, versions[foo.VERSION]
                )
            )
        perturbators[foo.VERSION] = foo.generate_source_params
        versions[foo.VERSION] = foo.__file__
    return perturbators


PERTURBATORS = __load_all__()
