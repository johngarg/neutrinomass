from setuptools import setup, find_packages
from setuptools.command.test import test as _test
from os import path

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
    name="mv",
    description="Example code behind arXiv:190x.xxxxx.",
    long_description=long_description,
    author="John Gargalionis",
    author_email="garj@student.unimelb.edu.au",
    license="BSD-2",
    packages=["mv", "mv.tensormethod", "mv.completions", "mv.database", "mv.analysis"],
    install_requires=["numpy", "sympy"],
    tests_require=["pytest"],
    package_data={},
    cmdclass={"test": PyTest},
)
