from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name='repro-image',
    version='0.1.0',
    packages=find_packages('src', exclude=['tests', 'tests.*', 'docs']),
    package_dir={'': 'src'},
    url='https://github.com/virtualpregnancy/repro-image.git',
    license=license,
    author='ABI Pregnancy Modelling Group',
    author_email='alys.clark@auckland.ac.nz',
    test_suite='nose.collector',
    tests_require=['nose'],
    description=''
)
