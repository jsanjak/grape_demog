from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
import sys
import setuptools
import os,glob
import fwdpy11

__version__ = '0.0'
class get_pybind_include(object):
    """Helper class to determine the pybind11 include path

    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __init__(self, user=False):
        self.user = user

    def __str__(self):
        import pybind11
        return pybind11.get_include(self.user)

INCLUDES=[
    # Path to pybind11 headers
    get_pybind_include(),
    get_pybind_include(user=True),
    # Path to fwdpy11 and fwdpp headers
    fwdpy11.get_includes(),
    fwdpy11.get_fwdpp_includes(),
    os.path.join(sys.prefix, 'include')
]

LIBRARY_DIRS=[
    os.path.join(sys.prefix, 'lib')
    ]

ext_modules = [
    Extension(
        'clonal.cl_evolve',
        ['clonal/src/cl_evolve.cpp'],
        library_dirs=LIBRARY_DIRS,
        include_dirs=INCLUDES,
        libraries=['gsl','gslcblas'],
        language='c++'
    )
]

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler):
    """Return the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-stdlib=libc++', '-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if ct == 'unix':
            opts.append('-DVERSION_INFO="%s"' % self.distribution.get_version())
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('/DVERSION_INFO=\\"%s\\"' % self.distribution.get_version())
        for ext in self.extensions:
            ext.extra_compile_args = opts
        build_ext.build_extensions(self)

setup(name="clonal",
      version=__version__,
      description="A clonal evolution function extension for fwdpy11",
      url="https://github.com/ThorntonLab/clonalevolve",
      author="Jaleal Sanjak",
      author_email="jsanjak@uci.edu",
      license="GPLv3",
      ext_modules=ext_modules,
      cmdclass={'build_ext': BuildExt},
      packages=['clonal'],
      requires=['pybind11','numpy','fwdpy11'],
      #install_requires=['pybind11>=1.7','numpy','fwdpy11'],
      zip_safe=False)

