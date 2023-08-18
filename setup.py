from setuptools import setup, find_packages
from scProportionTest import __version__

setup(
    name='py-scProportionTest',
    version=__version__,
    url='https://github.com/SamuelAMiller1/py-scProportionTest',
    author='Sam Miller',
    packages=find_packages(),    
    install_requires=['numpy', 'pandas', 'statsmodels', 'matplotlib', 'tqdm'],
    description='Python package to evaluate differences in cell type proportions',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown'
)
