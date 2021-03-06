# from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages

module1 = Extension('PyReCoDe',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['../../include','../../src'],
                    libraries = ['zlib'],
                    library_dirs = ['../../lib'],
                    sources = ['pyrecode.cpp'])

setup (name = 'PyReCoDe',
       version = '1.0',
       description = '',
       author = 'Abhik Datta',
       author_email = '',
       url = '',
       long_description = 'Contains readers and writers for ReCoDe',
       ext_modules = [module1],
       packages=find_packages()
    )