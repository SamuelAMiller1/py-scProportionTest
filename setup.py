from setuptools import setup, find_packages

setup(
    name='py-scProportionTest',
    version='0.1.0',
    url='https://github.com/SamuelAMiller1/py-scProportionTest',
    author='Sam Miller',
    packages=find_packages(),    
    install_requires=['numpy', 'pandas', 'statsmodels', 'matplotlib', 'tqdm'],
    description='Python package to evaluate differences in cell type proportions',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown'
)
