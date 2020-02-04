# from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages

module1 = Extension('c_recode',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['include','src'],
                    libraries = [],
                    library_dirs = ['lib'],
                    sources = ['src/pyrecode/pyrecode.cpp'],
                    extra_compile_args=["-lz", "-O3"],
                    extra_link_args=["-lz", "-L/usr/lib64/libz.so"])
# "-L/usr/lib64/libz.so"


setup (name = 'pyrecode',
       version = '1.0',
       description = '',
       author = 'Abhik Datta',
       author_email = '',
       url = 'https://github.com/NDLOHGRP/ReCoDe',
       long_description = 'Readers and writers for ReCoDe',
       package_dir={'': 'src'},
       ext_modules = [module1],
       packages=find_packages('src')
    )