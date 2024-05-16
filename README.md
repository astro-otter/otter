# OTTER API
### **O**pen mul**T**iwavelength **T**ransient **E**vent **R**epository

A Python API for the OTTER.

[actions-badge]:            https://github.com/astro-otter/otter/workflows/CI/badge.svg
[actions-link]:             https://github.com/astro-otter/otter/actions
[black-badge]:              https://img.shields.io/badge/code%20style-black-000000.svg
[black-link]:               https://github.com/psf/black
[conda-badge]:              https://img.shields.io/conda/vn/conda-forge/hepfile
[conda-link]:               https://github.com/conda-forge/hepfile-feedstock
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/mattbellis/hepfile/discussions
[gitter-badge]:             https://badges.gitter.im/https://github.com/mattbellis/hepfile/community.svg
[gitter-link]:              https://gitter.im/https://github.com/mattbellis/hepfile/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
[pypi-link]:                https://pypi.org/project/astro-otter/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/astro-otter
[pypi-version]:             https://badge.fury.io/py/astro-otter.svg
[rtd-badge]:                https://readthedocs.org/projects/otter/badge/?version=latest
[rtd-link]:                 https://otter.readthedocs.io/en/latest/?badge=latest
[sk-badge]:                 https://scikit-hep.org/assets/images/Scikit--HEP-Project-blue.svg
[ruff-badge]:               https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json
[ruff-link]:                https://github.com/astral-sh/ruff
[codecov-badge]:            https://codecov.io/gh/astro-otter/otter/graph/badge.svg?token=BtCerOdTc0
[codecov-link]:             https://codecov.io/gh/astro-otter/otter

[![Documentation Status](https://readthedocs.org/projects/astro-otter/badge/?version=latest)](https://astro-otter.readthedocs.io/en/latest/?badge=latest)
[![Actions Status][actions-badge]][actions-link]
[![PyPI version][pypi-version]][pypi-link]
[![Linting: Ruff][ruff-badge]][ruff-link]
[![codecov][codecov-badge]][codecov-link]

## Installation
To install the OTTER API use
```
git clone https://github.com/astro-otter/otter.git
cd otter
python -m pip install .
```
This will be changed into the more convenient `python -m pip install astro-otter` at a later date!

For developers, please also enable the pre-commit hooks using
```
pre-commit install
```

## Tutorial
### Connecting to the OTTER
```python
# import the API
from otter import Otter, Transient
```


```python
# connect to the database
# this username and password is just for now and will be updated later!
db = Otter(username='user@otter', password='insecure')
```

### A typical workflow

First use `Otter.getMeta` to query

```python
# can query by ANY name associated with an object
db.getMeta(names=['ASASSN-15oi', 'AT2020opy'])
```




    [{'name': {'default_name': 'ASASSN-15oi', 'alias': [{'value': 'ASASSN-15oi', 'reference': 'ASASSN'}]}, 'coordinate': {'equitorial': [{'ra': '20 39 09.096', 'dec': '-30 45 20.71', 'epoch': 'J2000', 'system': 'ICRS', 'ra_units': 'hourangle', 'dec_units': 'deg', 'reference': ['2021NatAs...5..491H'], 'computed': False, 'default': True, 'uuid': 'a06be641-1601-4737-9a1a-bd25c5dd61e6'}], 'galactic': [{'l': 13.01154485751856, 'b': -35.41877256185317, 'l_units': 'deg', 'b_units': 'deg', 'reference': 'a06be641-1601-4737-9a1a-bd25c5dd61e6', 'computed': True}]}, 'epoch': {'date_discovery': [{'value': 57248.2, 'date_format': 'MJD', 'reference': ['2021NatAs...5..491H'], 'computed': False}]}, 'distance': {'redshift': [{'value': '0.0484', 'reference': ['2021NatAs...5..491H'], 'computed': False, 'default': True}]}, 'classification': [{'object_class': 'TDE', 'confidence': 1, 'reference': ['2021NatAs...5..491H'], 'default': True}]},
     {'name': {'default_name': 'AT2020opy', 'alias': [{'value': 'AT2020opy', 'reference': 'TNS'}]}, 'coordinate': {'equitorial': [{'ra': '15 56 25.728', 'dec': '+23 22 21.15', 'epoch': 'J2000', 'system': 'ICRS', 'ra_units': 'hourangle', 'dec_units': 'deg', 'reference': ['2023MNRAS.518..847G'], 'computed': False, 'default': True, 'uuid': '4f414ede-e0f0-4423-b2cc-f3ad309f0936'}], 'galactic': [{'l': 38.411569358993255, 'b': 48.23253616380999, 'l_units': 'deg', 'b_units': 'deg', 'reference': '4f414ede-e0f0-4423-b2cc-f3ad309f0936', 'computed': True}]}, 'epoch': {'date_discovery': [{'value': 59038.23, 'date_format': 'MJD', 'reference': ['2023MNRAS.518..847G'], 'computed': False}]}, 'distance': {'redshift': [{'value': '0.159', 'reference': ['2023MNRAS.518..847G'], 'computed': False, 'default': True}]}, 'classification': [{'object_class': 'TDE', 'confidence': 1, 'reference': ['2023MNRAS.518..847G'], 'default': True}]}]



We can also do a cone search
```python
from astropy.coordinates import SkyCoord
import astropy.units as u
coord = SkyCoord(239, 23, unit=('deg', 'deg'))
rad = (1*u.deg).to(u.arcsec).value
db.getMeta(coords=coord, radius=rad)
```




    [{'name': {'default_name': 'AT2020opy', 'alias': [{'value': 'AT2020opy', 'reference': 'TNS'}]}, 'coordinate': {'equitorial': [{'ra': '15 56 25.728', 'dec': '+23 22 21.15', 'epoch': 'J2000', 'system': 'ICRS', 'ra_units': 'hourangle', 'dec_units': 'deg', 'reference': ['2023MNRAS.518..847G'], 'computed': False, 'default': True, 'uuid': '4f414ede-e0f0-4423-b2cc-f3ad309f0936'}], 'galactic': [{'l': 38.411569358993255, 'b': 48.23253616380999, 'l_units': 'deg', 'b_units': 'deg', 'reference': '4f414ede-e0f0-4423-b2cc-f3ad309f0936', 'computed': True}]}, 'epoch': {'date_discovery': [{'value': 59038.23, 'date_format': 'MJD', 'reference': ['2023MNRAS.518..847G'], 'computed': False}]}, 'distance': {'redshift': [{'value': '0.159', 'reference': ['2023MNRAS.518..847G'], 'computed': False, 'default': True}]}, 'classification': [{'object_class': 'TDE', 'confidence': 1, 'reference': ['2023MNRAS.518..847G'], 'default': True}]}]



Or search within a redshift range (or just a maximum or minimum)
```python
# can search a redshift range
db.getMeta(minZ=0.5, maxZ=0.9)
```




    [{'name': {'default_name': 'Sw J1112-82', 'alias': [{'value': 'Sw J1112-82', 'reference': 'Swift'}]}, 'coordinate': {'equitorial': [{'ra': '11 11 47.6', 'dec': '-82 38 44.44', 'epoch': 'J2000', 'system': 'ICRS', 'ra_units': 'hourangle', 'dec_units': 'deg', 'reference': ['2017MNRAS.472.4469B'], 'computed': False, 'default': True, 'uuid': '3bb02aa4-eb60-4f09-acdf-d1431a6addc4'}], 'galactic': [{'l': 299.6337165869647, 'b': -20.420594756871665, 'l_units': 'deg', 'b_units': 'deg', 'reference': '3bb02aa4-eb60-4f09-acdf-d1431a6addc4', 'computed': True}]}, 'epoch': {'date_discovery': [{'value': '55729.5', 'date_format': 'MJD', 'reference': ['2017MNRAS.472.4469B'], 'computed': False}]}, 'distance': {'redshift': [{'value': '0.89', 'reference': ['2017MNRAS.472.4469B'], 'computed': False, 'default': True}]}, 'classification': [{'object_class': 'TDE', 'confidence': 1, 'reference': ['2017MNRAS.472.4469B'], 'default': True}]}]



We can even get all objects that have spectra associated with them
```python
# just get objects that have spectra associated with them
db.getMeta(hasSpec=True)
```




    []

This should be empty because at the time of developing this tutorial there were no spectra in
OTTER. Similarly, we can get all objects that have photometry associated with them with
`db.getMeta(hasPhot=True)`.

These outputs may appear like dictionaries but they're actually customized!

Besides the typical dictionary methods the following methods are also implemented for Transient objects.

```python
help(Transient)
```

    Help on class Transient in module otter.transient:

    class Transient(collections.abc.MutableMapping)
     |  Transient(d={}, name=None)
     |
     |  Method resolution order:
     |      Transient
     |      collections.abc.MutableMapping
     |      collections.abc.Mapping
     |      collections.abc.Collection
     |      collections.abc.Sized
     |      collections.abc.Iterable
     |      collections.abc.Container
     |      builtins.object
     |
     |  Methods defined here:
     |
     |  __add__(self, other, strict_merge=True)
     |      Merge this transient object with another transient object
     |
     |      Args:
     |          other [Transient]: A Transient object to merge with
     |          strict_merge [bool]: If True it won't let you merge objects that
     |                               intuitively shouldn't be merged (ie. different
     |                               transient events).
     |
     |  __delitem__(self, keys)
     |
     |  __getitem__(self, keys)
     |
     |  __init__(self, d={}, name=None)
     |      Overwrite the dictionary init
     |
     |      Args:
     |          d [dict]: A transient dictionary
     |
     |  __iter__(self)
     |
     |  __len__(self)
     |
     |  __repr__(self, html=False)
     |      Return repr(self).
     |
     |  __setitem__(self, key, value)
     |
     |  cleanPhotometry(self, flux_unit='mag(AB)', date_unit='MJD')
     |      Ensure the photometry associated with this transient is all in the same units/system/etc
     |
     |  getMeta(self, keys=None)
     |      Get the metadata (no photometry or spectra)
     |
     |      This essentially just wraps on __getitem__ but with some checks
     |
     |      Args:
     |          keys [list[str]] : list of keys
     |
     |  getSkyCoord(self, coord_type='equitorial', idx=0)
     |      Convert the coordinates to an astropy SkyCoord
     |
     |  keys(self)
     |      D.keys() -> a set-like object providing a view on D's keys
     |
     |  plotPhotometry(self, flux_unit='mag(AB)', date_unit='datetime', **kwargs)
     |      Plot the photometry associated with this transient (if any)
     |
     |      Args:
     |          flux_unit [str]: Valid astropy unit string for the flux (y-axis) units.
     |                           Default: 'ABmag'
     |          date_unit [str]: Valid astropy unit string for the date (x-axis) units.
     |                           Default: 'MJD'


Some other advantages of the Transient objects are
```python
t = db.getMeta(minZ=0.5, maxZ=0.9)[0]
print(type(t))
print()
# say you want to get the equitorial coordinates
# you can do it classically
print(t['coordinate']['equitorial'])

# or you can use the hdf5 style
print(t['coordinate/equitorial'])
```

    <class 'otter.transient.Transient'>

    [{'ra': '11 11 47.6', 'dec': '-82 38 44.44', 'epoch': 'J2000', 'system': 'ICRS', 'ra_units': 'hourangle', 'dec_units': 'deg', 'reference': ['2017MNRAS.472.4469B'], 'computed': False, 'default': True, 'uuid': '3bb02aa4-eb60-4f09-acdf-d1431a6addc4'}]
    [{'ra': '11 11 47.6', 'dec': '-82 38 44.44', 'epoch': 'J2000', 'system': 'ICRS', 'ra_units': 'hourangle', 'dec_units': 'deg', 'reference': ['2017MNRAS.472.4469B'], 'computed': False, 'default': True, 'uuid': '3bb02aa4-eb60-4f09-acdf-d1431a6addc4'}]


You can even get multiple fields at once like with astropy Tables or pandas DataFrames
```python
# You can also get multiple fields at once
t[['name/default_name', 'coordinate/equitorial', 'distance']]
```




    {'name/default_name': 'Sw J1112-82', 'coordinate/equitorial': [{'ra': '11 11 47.6', 'dec': '-82 38 44.44', 'epoch': 'J2000', 'system': 'ICRS', 'ra_units': 'hourangle', 'dec_units': 'deg', 'reference': ['2017MNRAS.472.4469B'], 'computed': False, 'default': True, 'uuid': '3bb02aa4-eb60-4f09-acdf-d1431a6addc4'}], 'distance': {'redshift': [{'value': '0.89', 'reference': ['2017MNRAS.472.4469B'], 'computed': False, 'default': True}]}}



You can also add two Transient objects to merge them
```python
t1, t2 = db.query(names=['ASASSN-15oi', 'AT2020opy'])

try:
    t1 + t2
except ValueError as ve:
    print('The following error is actually expected!')
    print('We dont want you to be able to combine any old transients!')
    print()
    print('Error Message:')
    print(ve)
```

    The following error is actually expected!
    We dont want you to be able to combine any old transients!

    Error Message:
    These two transients are not within 5 arcseconds! They probably do not belong together! If they do You can set strict_merge=False to override the check


This error message is expected! If you want to override it then you can do
```python
t2['photometry'][0]['reference'] = '2021NatAs...5..491H'

t3 = t1.__add__(t2, strict_merge=False)
```

Obviously, this result doesn't makes sense! This has the data from two completely different transients in it. So, be careful using `strict_merge=False`!

### Can then get photometry
This does the conversion for you!!!


```python
db.getPhot?
```
    Get the photometry of the objects matching the arguments. This will do the
    unit conversion for you!

    Args:
        flux_units [astropy.unit.Unit]: Either a valid string to convert
                                        or an astropy.unit.Unit
        date_units [astropy.unit.Unit]: Either a valid string to convert to a date
                                        or an astropy.unit.Unit
        return_type [str]: Either 'astropy' or 'pandas'. If astropy, returns an
                           astropy Table. If pandas, returns a pandas DataFrame.
                           Default is 'astropy'.

        **kwargs : Arguments to pass to Otter.query(). Can be:
                   names [list[str]]: A list of names to get the metadata for
                   coords [SkyCoord]: An astropy SkyCoord object with coordinates to match to
                   radius [float]: The radius in arcseconds for a cone search, default is 0.05"
                   minZ [float]: The minimum redshift to search for
                   maxZ [float]: The maximum redshift to search for
                   refs [list[str]]: A list of ads bibcodes to match to. Will only return
                          metadata for transients that have this as a reference.
                   hasSpec [bool]: if True, only return transients that have spectra.

    Return:
       The photometry for the requested transients that match the arguments.
       Will be an astropy Table sorted by transient default name.


This means you can easily grab photometry in consistent units to plot!
```python
import matplotlib.pyplot as plt
flux_unit = u.ABmag #u.erg/u.s/u.cm**2/u.Hz #'erg/s/cm^2/Hz'
tab = db.getPhot(flux_unit=flux_unit, date_unit='datetime', names=['ASASSN-15oi', 'ASASSN-14li'], return_type='pandas')

tab
```




<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>reference</th>
      <th>date</th>
      <th>filter_key</th>
      <th>computed</th>
      <th>obs_type</th>
      <th>upperlimit</th>
      <th>freq_eff</th>
      <th>freq_units</th>
      <th>human_readable_refs</th>
      <th>converted_flux</th>
      <th>converted_date</th>
      <th>name</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2016ApJ...819L..25A</td>
      <td>57124.871000</td>
      <td>5.0GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>5.0</td>
      <td>GHz</td>
      <td>Alexander et al. (2016)</td>
      <td>15.697417</td>
      <td>2015-04-12 20:54:14.400000</td>
      <td>ASASSN-14li</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2016ApJ...819L..25A</td>
      <td>57190.830000</td>
      <td>5.0GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>5.0</td>
      <td>GHz</td>
      <td>Alexander et al. (2016)</td>
      <td>15.797567</td>
      <td>2015-06-17 19:55:12.000000</td>
      <td>ASASSN-14li</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2016ApJ...819L..25A</td>
      <td>57229.750000</td>
      <td>5.0GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>5.0</td>
      <td>GHz</td>
      <td>Alexander et al. (2016)</td>
      <td>15.915797</td>
      <td>2015-07-26 18:00:00.000000</td>
      <td>ASASSN-14li</td>
    </tr>
    <tr>
      <th>3</th>
      <td>2016ApJ...819L..25A</td>
      <td>57286.514583</td>
      <td>5.0GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>5.0</td>
      <td>GHz</td>
      <td>Alexander et al. (2016)</td>
      <td>16.152094</td>
      <td>2015-09-21 12:20:59.997120</td>
      <td>ASASSN-14li</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2016ApJ...819L..25A</td>
      <td>57362.700000</td>
      <td>5.0GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>5.0</td>
      <td>GHz</td>
      <td>Alexander et al. (2016)</td>
      <td>16.547278</td>
      <td>2015-12-06 16:48:00.000000</td>
      <td>ASASSN-14li</td>
    </tr>
    <tr>
      <th>0</th>
      <td>2021NatAs...5..491H</td>
      <td>57256.200000</td>
      <td>6.1GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>True</td>
      <td>6.1</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>20.103715</td>
      <td>2015-08-22 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2021NatAs...5..491H</td>
      <td>57271.200000</td>
      <td>6.1GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>True</td>
      <td>6.1</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>20.009244</td>
      <td>2015-09-06 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2021NatAs...5..491H</td>
      <td>57338.200000</td>
      <td>6.1GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>True</td>
      <td>6.1</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>19.454622</td>
      <td>2015-11-12 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
    <tr>
      <th>3</th>
      <td>2021NatAs...5..491H</td>
      <td>57430.200000</td>
      <td>4.8GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>4.8</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>16.282787</td>
      <td>2016-02-12 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2021NatAs...5..491H</td>
      <td>57438.200000</td>
      <td>4.8GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>4.8</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>16.515601</td>
      <td>2016-02-20 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
    <tr>
      <th>5</th>
      <td>2021NatAs...5..491H</td>
      <td>57445.200000</td>
      <td>5.5GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>5.5</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>16.537560</td>
      <td>2016-02-27 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
    <tr>
      <th>6</th>
      <td>2021NatAs...5..491H</td>
      <td>57481.200000</td>
      <td>5.5GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>5.5</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>16.610182</td>
      <td>2016-04-03 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
    <tr>
      <th>7</th>
      <td>2021NatAs...5..491H</td>
      <td>57531.200000</td>
      <td>5.0GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>5.0</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>16.571355</td>
      <td>2016-05-23 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
    <tr>
      <th>8</th>
      <td>2021NatAs...5..491H</td>
      <td>57617.200000</td>
      <td>5.0GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>5.0</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>16.924287</td>
      <td>2016-08-17 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
    <tr>
      <th>9</th>
      <td>2021NatAs...5..491H</td>
      <td>57824.200000</td>
      <td>5.0GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>5.0</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>17.854247</td>
      <td>2017-03-12 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
    <tr>
      <th>10</th>
      <td>2021NatAs...5..491H</td>
      <td>58665.200000</td>
      <td>3.0GHz</td>
      <td>False</td>
      <td>radio</td>
      <td>False</td>
      <td>3.0</td>
      <td>GHz</td>
      <td>Horesh et al. (2021)</td>
      <td>14.142275</td>
      <td>2019-07-01 04:48:00.000000</td>
      <td>ASASSN-15oi</td>
    </tr>
  </tbody>
</table>
</div>




```python
fig, ax = plt.subplots()
for key, table in tab.groupby('name'):
    ax.plot(table['converted_date'], table['converted_flux'], label=key, marker='o', linestyle='none')

ax.set_ylabel(f'Flux Density [{flux_unit}]')
ax.set_xlabel('Date')
ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45)
ax.legend();
```


![png](api-test_files/api-test_17_0.png)


### General Queries (shouldn't be used unless you know what you're doing!)

The structure of a more generalized query to OTTER is
```python
# General queries
Otter.query?
```
    Wraps on the super.AQLQuery and queries the OTTER database more intuitively.

    WARNING! This does not do any conversions for you!
    This is how it differs from the `getMeta` method. Users should prefer to use
    `getMeta`, `getPhot`, and `getSpec` independently because it is a better
    workflow and can return the data in an astropy table with everything in the
    same units.

    Args:
        names [list[str]]: A list of names to get the metadata for
        coords [SkyCoord]: An astropy SkyCoord object with coordinates to match to
        radius [float]: The radius in arcseconds for a cone search, default is 0.05"
        minZ [float]: The minimum redshift to search for
        maxZ [float]: The maximum redshift to search for
        refs [list[str]]: A list of ads bibcodes to match to. Will only return
                          metadata for transients that have this as a reference.
        hasPhot [bool]: if True, only returns transients which have photometry.
        hasSpec [bool]: if True, only return transients that have spectra.

    Return:
       Get all of the raw (unconverted!) json data for objects that match the criteria.


An example of this is
```python
res = db.query(names=['ASASSN-15oi', 'AT2020opy'])
print(len(res))
print(type(res[0]))
```

    2
    <class 'otter.transient.Transient'>



```python
res[0].keys()
```




    dict_keys(['_key', '_id', '_rev', 'schema_version', 'name', 'coordinate', 'distance', 'classification', 'reference_alias', 'epoch', 'photometry', 'filter_alias'])



However, Notice how this is simply the raw results!!! This means you need to be very careful with the results you get from these queries because no conversion is done!
```python
res[0]
```




    {'_key': '3776907', '_id': 'tdes/3776907', '_rev': '_hD8NTy2---', 'schema_version': {'value': '0', 'comment': 'Original Dataset'}, 'name': {'default_name': 'ASASSN-15oi', 'alias': [{'value': 'ASASSN-15oi', 'reference': 'ASASSN'}]}, 'coordinate': {'equitorial': [{'ra': '20 39 09.096', 'dec': '-30 45 20.71', 'epoch': 'J2000', 'system': 'ICRS', 'ra_units': 'hourangle', 'dec_units': 'deg', 'reference': ['2021NatAs...5..491H'], 'computed': False, 'default': True, 'uuid': 'a06be641-1601-4737-9a1a-bd25c5dd61e6'}], 'galactic': [{'l': 13.01154485751856, 'b': -35.41877256185317, 'l_units': 'deg', 'b_units': 'deg', 'reference': 'a06be641-1601-4737-9a1a-bd25c5dd61e6', 'computed': True}]}, 'distance': {'redshift': [{'value': '0.0484', 'reference': ['2021NatAs...5..491H'], 'computed': False, 'default': True}]}, 'classification': [{'object_class': 'TDE', 'confidence': 1, 'reference': ['2021NatAs...5..491H'], 'default': True}], 'reference_alias': [{'name': '2021NatAs...5..491H', 'human_readable_name': 'Horesh et al. (2021)'}], 'epoch': {'date_discovery': [{'value': 57248.2, 'date_format': 'MJD', 'reference': ['2021NatAs...5..491H'], 'computed': False}]}, 'photometry': [{'reference': '2021NatAs...5..491H', 'raw': [0.033, 0.036, 0.06, 1.114, 0.899, 0.881, 0.824, 0.854, 0.617, 0.262, 8], 'raw_units': ['mJy', 'mJy', 'mJy', 'mJy', 'mJy', 'mJy', 'mJy', 'mJy', 'mJy', 'mJy', 'mJy'], 'date': [57256.2, 57271.2, 57338.2, 57430.2, 57438.2, 57445.2, 57481.2, 57531.2, 57617.2, 57824.2, 58665.2], 'date_format': ['MJD', 'MJD', 'MJD', 'MJD', 'MJD', 'MJD', 'MJD', 'MJD', 'MJD', 'MJD', 'MJD'], 'filter_key': ['6.1GHz', '6.1GHz', '6.1GHz', '4.8GHz', '4.8GHz', '5.5GHz', '5.5GHz', '5.0GHz', '5.0GHz', '5.0GHz', '3.0GHz'], 'computed': [False, False, False, False, False, False, False, False, False, False, False], 'obs_type': ['radio', 'radio', 'radio', 'radio', 'radio', 'radio', 'radio', 'radio', 'radio', 'radio', 'radio'], 'upperlimit': [True, True, True, False, False, False, False, False, False, False, False]}], 'filter_alias': [{'filter_key': '6.1GHz', 'freq_eff': 6.1, 'freq_units': 'GHz'}, {'filter_key': '4.8GHz', 'freq_eff': 4.8, 'freq_units': 'GHz'}, {'filter_key': '5.5GHz', 'freq_eff': 5.5, 'freq_units': 'GHz'}, {'filter_key': '5.0GHz', 'freq_eff': 5, 'freq_units': 'GHz'}, {'filter_key': '3.0GHz', 'freq_eff': 3, 'freq_units': 'GHz'}]}



### Some helpful methods

To get the Astropy SkyCoord object for a specific transient you can use
```python
skycoord = t1.getSkyCoord()
skycoord
```




    <SkyCoord (ICRS): (ra, dec) in deg
        (309.7879, -30.75575278)>



Or, to get the html code for plotting the photometry for a transient
```python
html = res[0].plotPhotometry()
html[:100] # only show the first part of it
```




    '<div>                        <script type="text/javascript">window.PlotlyConfig = {MathJaxConfig: \'l'



### Uploading/Editing Data
If you have admin access to OTTER you can also upload new data! An example of this is that we first
need to create a new Transient object.

```python
from otter import Otter, Transient
from copy import deepcopy
from collections import Counter
import awkward as ak
import warnings
import numpy as np
import re
from astropy.coordinates import SkyCoord
import json

# generate some test cases
db = Otter()
t1 = db.query(names='2022xkq')[0] #
t2 = deepcopy(t1)
print(t1.keys())

# change t2 for testing
t2['name'] = {'default_name':'2022xkq',
             'alias': [{'value':'foo', 'reference': 'x'},
                      {'value': '2022xkq', 'reference': 'x'}]}
t2['reference_alias'].append({'name': 'x',
   'human_readable_name': 'test, name (year)'}) # add an extra value
del t2['photometry']
t2['for_test'] = {'test': 'bar'} # add a test key that isn't in t1
t2['coordinate/equitorial'][0]['reference'] = 'noah'
t2['filter_alias'].append({'filter_key': 'foo'})
t2['schema_version/value'] = 100
t2['epoch'] = {'date_peak': [{'value': 56983,
    'date_format': 'MJD',
    'reference': ['2016ApJ...819L..25A',
     '2016Sci...351...62V',
     '2016ApJ...832L..10R',
     '2018MNRAS.475.4011B'],
    'computed': False}],

               'date_discovery': [{'value': 56983,
    'date_format': 'MJD',
    'reference': ['2016ApJ...819L..25A',
     '2016Sci...351...62V',
     '2016ApJ...832L..10R',
     '2018MNRAS.475.4011B'],
    'computed': False}],


              'date_discovery': [{'value': 56984,
    'date_format': 'MJD',
    'reference': ['2016ApJ...819L..25A',
     '2016Sci...351...62V',
     '2016ApJ...832L..10R',
     '2018MNRAS.475.4011B'],
    'computed': False}]
              }

t2['distance'] = {
    "redshift": [
      {
        "value": "0.0207",
        "reference": [
          "Noah"
        ],
        "computed": False
      },
        {
        "value": "0.02",
        "reference": [
          "Noah"
        ],
        "computed": False
      }
    ],

    "dispersion_measure": [
      {
        "value": "0.0206",
        "reference": [
          "Noah"
        ],
        "computed": False
      }
    ]
  }

t2['classification'] = [{'object_class':'SN',
                        'confidence': 1,
                         'reference': 'Noah'
                        }]

t2['photometry'] = {'phot_0': {'telescope': 'Noahs Telescope',
                               'reference': 'Noah',
                               'flux': [{'filter': 'z',
                                 'telescope': 'Noahs Telescope',
                                 'upperlimit': True,
                                 'date': 59864.4914116667,
                                 'date_format': 'MJD',
                                 'raw': 20.01,
                                 'raw_units': 'mag(AB)',
                                 'filter_key': 'NoahsTelescope.z',
                                 'obs_type': 'uvoir'}]},
                    'phot_1': {'telescope': 'CAHA',
                               'reference': 'Noah',
                               'flux': [{'filter': 'H',
                                 'telescope': 'CAHA',
                                 'upperlimit': False,
                                 'date': 59898.12077,
                                 'date_format': 'MJD',
                                 'raw': 14.87048,
                                 'raw_err': 0.0187,
                                 'raw_units': 'mag(AB)',
                                 'filter_key': 'CAHA.H',
                                 'obs_type': 'uvoir'}]}

                               }
```

Now we need to connect to OTTER with admin access and upload the new transient object.
```python
db = Otter(username='admin@otter', password='insecure')
#db.upload(t2)
```


## Repo Organization
| Directory | Contents |
|------------|------------|
| `py/otter` | A pip installable API for interfacing with the ArangoDB database|
