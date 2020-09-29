from setuptools import setup
from setuptools.command.test import test as _test
from os import path
import sys

here = path.abspath(path.dirname(__file__))


class PyTest(_test):
    user_options = [("pytest-args=", "a", "Arguments to pass to pytest")]

    def initialize_options(self):
        _test.initialize_options(self)
        self.pytest_args = []

    def run_tests(self):
        import shlex
        import pytest

        if self.pytest_args:
            args = shlex.split(self.pytest_args)
        else:
            args = None
            errno = pytest.main(args)
            sys.exit(errno)


with open(path.join(here, "README.org"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="neutrinomass",
    description="Code behind arXiv:20xx.xxxxx.",
    # long_description=long_description,
    author="John Gargalionis",
    author_email="garj@student.unimelb.edu.au",
    version="1.0.1",
    url="https://github.com/johngarg/neutrinomass",
    license="MIT",
    packages=[
        "neutrinomass",
        "neutrinomass.tensormethod",
        "neutrinomass.completions",
        "neutrinomass.database",
        "neutrinomass.analysis",
        "neutrinomass.utils",
    ],
    install_requires=[
        "matplotlib>=3.0.2",
        "pytest>=3.8.0",
        "basisgen>=1.0.1",
        "networkx>=2.2",
        "sympy==1.2",
        "alive_progress>=1.5.1",
        "matchpy>=0.5.2",
        "pandas>=1.1.2",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    tests_require=["pytest"],
    package_data={},
    cmdclass={"test": PyTest},
    include_package_data=True,
    python_requires=">=3.6",
)
