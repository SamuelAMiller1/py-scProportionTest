from setuptools import setup, find_packages
from scProportionTest import __version__

with open('README.rst', 'r') as file:
    long_description = file.read()
setup(
    name='py-scProportionTest',
    version=__version__,
    url='https://github.com/SamuelAMiller1/py-scProportionTest',
    author='Sam Miller',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'pandas',
        'statsmodels',
        'matplotlib',
        'tqdm',
        'seaborn',
    ],
    description='Python package to evaluate differences in cell type proportions',
    long_description=long_description,
    long_description_content_type='text/x-rst'
)

