===========
specialsoss
===========


.. image:: https://img.shields.io/pypi/v/specialsoss.svg
        :target: https://pypi.python.org/pypi/specialsoss

.. image:: https://img.shields.io/travis/hover2pi/specialsoss.svg
        :target: https://travis-ci.org/hover2pi/specialsoss

.. image:: https://readthedocs.org/projects/specialsoss/badge/?version=latest
        :target: https://specialsoss.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://img.shields.io/coveralls/github/hover2pi/specialsoss.svg
        :target: https://coveralls.io/github/hover2pi/specialsoss

.. image:: https://pyup.io/repos/github/hover2pi/specialsoss/shield.svg
     :target: https://pyup.io/repos/github/hover2pi/specialsoss/
     :alt: Updates


SPECtral Image AnaLysis for SOSS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Authors: Joe Filippazzo, Kevin Stevenson

This pure Python 3.5+ package performs optimal spectral extraction routines
for the Single Object Slitless Spectroscopy (SOSS) mode of the
Near-Infrared Imager and Slitless Spectrograph (NIRISS) instrument
onboard the James Webb Space Telescope (JWST).

Additional resources:

- `Full documentation <https://specialsoss.readthedocs.io/en/latest/>`_
- `Jupyter notebook <https://github.com/spacetelescope/specialsoss/blob/master/notebooks/specialsoss_demo.ipynb>`_
- `Build history <https://travis-ci.org/hover2pi/specialsoss>`_


Extracting Spectra from SOSS Observations
-----------------------------------------

The headers in JWST data products provide almost all the information
needed to perform the spectral extraction, making the path to your data
the only required input. To extract time series 1D spectra, simply do

.. code:: python

   # Imports
   import numpy as np
   from specialsoss import SossObs
   from pkg_resources import resource_filename

   # Run the extraction
   data = resource_filename('specialsoss', 'files/soss_example.fits')
   obs = SossObs(data)

That’s it! Now we can take a look at the extracted time-series spectra:

.. code:: python

   obs.plot()

.. figure:: specialsoss/files/extracted_spectra.png
   :alt: SOSS Extraction
