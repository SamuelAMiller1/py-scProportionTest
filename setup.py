from setuptools import setup, find_packages

setup(
    name='py-scProportionTest',
    version='0.1.0',
    url='https://github.com/SamuelAMiller1/py-scProportionTest',
    author='Sam Miller',
    description='Python package equivalent to R package scProportionTest',
    packages=find_packages(),    
    install_requires=['numpy', 'pandas', 'anndata', 'scanpy'],
)
