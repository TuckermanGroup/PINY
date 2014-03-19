import sys
import glob

# PINY directories
dir_piny = '..'
dir_piny_include = dir_piny + '/include'
dir_libpiny = dir_piny + '/lib'


# global Cython directives
cython_directives = {
    'wraparound':False,
    'boundscheck': False,
    'cdivision': True,
    'profile': False,
    'embedsignature': True
}


#
# check requirements
#

missing_packages = []

try:
    from distutils.core import setup
    from distutils.extension import Extension
except ImportError:
    missing_packages.append('distutils')

try:
    import numpy as np
except ImportError:
    missing_packages.append('numpy')

try:
    from Cython.Distutils import build_ext
except ImportError:
    missing_packages.append('Cython')

if len(missing_packages) > 0:
    print 'Required packages not found:'
    for p in missing_packages:
        print p
    sys.exit(1)


#
# prepare extensions
#

ext_names = glob.glob('piny/*.pyx')

extensions = [
    Extension(
        name[:-4].replace('/', '.'),
        [name],
        include_dirs = [np.get_include(), dir_piny_include],
        library_dirs = [dir_libpiny],
        libraries = ['piny'],
        extra_compile_args=['-fopenmp'],
        extra_link_args=['-fopenmp'])
    for name in ext_names]

for e in extensions:
    e.cython_directives = cython_directives


#
# build
#

setup(
    name = 'PINY',
    version = '',
    description = '',
    author = 'Ondrej Marsalek',
    author_email = 'ondrej.marsalek@gmail.com',
    url = '',

    packages = ['pypiny'],
    ext_modules = extensions,
    cmdclass = {'build_ext': build_ext},
)
