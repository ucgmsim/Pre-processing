try:
    from setuptools import setup, find_namespace_packages, find_packages
except ImportError:
    print(
        "An error occurred while trying to import setuptools. "
        "Try running pip with the additional flag '-r requirements.txt'"
    )
    from sys import exit
    exit(1)

setup(
    name="SRF generation",
    version="19.9.1",
    packages=find_namespace_packages(include=["srf_generation.*"]),
    url="https://github.com/ucgmsim/Pre-processing",
    description="Srf generation code",
    install_requires=["numpy>=1.14.3", "scipy>=1.1.0"],
    include_package_data=True,
    scripts=[
        "srf_generation/source_parameter_generation/gcmt_to_realisation.py",
        "srf_generation/source_parameter_generation/generate_realisations_from_gcmt.py",
        "srf_generation/input_file_generation/realisation_to_srf.py",
        "srf_generation/input_file_generation/generate_srf_from_realisations.py",
    ],
)
