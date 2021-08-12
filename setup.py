from os import path

from setuptools import find_packages, setup

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'readme.md'), encoding='utf-8') as f:
    l_d = f.read()

# Get the version.
version = {}
with open("hydropandas/version.py") as fp:
    exec(fp.read(), version)

setup(
    name='hydropandas',
    version=version['__version__'],
    description='hydropandas module by Artesia',
    long_description=l_d,
    long_description_content_type="text/markdown",
    url='https://github.com/ArtesiaWater/hydropandas',
    author='Artesia',
    license='MIT',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Other Audience',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
    platforms='Windows, Mac OS-X',
    install_requires=['scipy',
                      "pandas>=1.0,<1.2.0",  # see PR #64
                      'tqdm',
                      'zeep',
                      'requests',
                      'geopandas'
                      ],
    packages=find_packages(exclude=[]),
    package_data={"hydropandas": ["data/*.json"]},
    include_package_data=True
)
