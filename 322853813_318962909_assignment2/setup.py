from setuptools import Extension, setup

module = Extension("mykmeanssp", sources=['kmeansmodule.c'])
setup(name='mykmeanssp',
      version='1.0.0.0',
      description='Python wrapper',
      ext_modules=[module])