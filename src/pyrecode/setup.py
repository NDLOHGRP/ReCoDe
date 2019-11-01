from distutils.core import setup, Extension

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
       long_description = 'This is really just a demo package.',
       ext_modules = [module1]
    )