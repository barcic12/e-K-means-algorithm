from setuptools import setup, Extension

setup(
  name="mykmeanssp",
  version="0.1.0",
  install_requires=['invoke'],
  ext_modules=[Extension("mykmeanssp", ["kmeans.c"])]
)