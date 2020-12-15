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
    install_requires=['numpy>=1.15',
                      'matplotlib>=3.0',
                      'pandas>=0.24',
                      'scipy>=1.2',
                      'geopandas >= 0.6.0',
                      'zeep >= 3.4.0',
                      'bokeh >= 2.0.2',
                      'branca >= 0.3.1',
                      'folium >= 0.8.2',
                      'pastas >= 0.13.0',
                      'flopy >= 3.3.0',
                      'xarray >= 0.15.0',
                      'netCDF4 >= 1.5.4',
                      'requests >= 2.22.0',
                      'shapely >= 1.7.0',
					  'cython >= 0.29.21',
                      'pyproj >= 2.6.0',
                      'tqdm >= 4.42.1'],
    packages=find_packages(exclude=[]),
    package_data={"hydropandas": ["data/*.json"]},
    include_package_data=True
)
