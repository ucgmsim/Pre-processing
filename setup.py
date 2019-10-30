from setuptools import setup, find_namespace_packages

setup(
    name="SRF generation",
    version="19.9.1",
    packages=find_namespace_packages(include=["srf_generation.*"]),
    url="https://github.com/ucgmsim/Pre-processing",
    description="Srf generation code",
    install_requires=["numpy>=1.14.3", "scipy>=1.1.0"],
    scripts=[
        "srf_generation/source_parameter_generation/gcmt_to_srfgen.py",
        "srf_generation/source_parameter_generation/srfgenparams_to_srf.py"
     ],
)