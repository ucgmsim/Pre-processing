from setuptools import setup, find_packages

setup(
    name="SrfGen",
    version="19.9.1",
    packages=find_packages(),
    url="https://github.com/ucgmsim/Pre-processing",
    description="Srf generation code",
    install_requires=["numpy>=1.14.3", "scipy>=1.1.0"],
)
