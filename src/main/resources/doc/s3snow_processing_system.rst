.. role:: underline
    :class: underline

.. index:: SNAP S3-SNOW Processors

.. _s3snow_processing_system:

============================================
The SNAP S3-SNOW Processors
============================================

Overview
========

The key goal of the S3-SNOW project regarding software development, production and dissemination was to
implement the proposed algorithms for several  key  snow  parameters from  Sentinel-3 OLCI data in free and easily
accessible open source toolboxes, notably and foremost ESA’s SNAP toolbox.
During the implementation process, all processing software was distributed within the consortium in frequent cycles
for the purpose of a comprehensive validation from an appropriate OLCI database containing products covering a variety of
selected snow-covered areas. All SNAP
S3-SNOW processors are available as SNAP plugins and can be run within SNAP on any supported platform
(Windows, Linus, MacOS). The procedure for installation and operation is described in this chapter.

.. index:: Theoretical Background Summary

Theoretical Background
======================

The motivation and theoretical background for the retrieval of key snow parameters is summarized in the S3-SNOW project
proposal [`1 <intro.html#References>`_].
The underlying algorithms are described in detail in the corresponding project ATBD  [`2 <intro.html#References>`_].

.. index:: Processing Environment

Processing Environment
======================

As said, the S3-SNOW processors are available as SNAP plugins and can be run within SNAP on any supported platform
(Windows, Linus, MacOS).
Section :ref:`s3snow_installation` describes in more detail how to install the plugins in SNAP.

.. index:: Processing Components

Processing Components
=====================

The SNAP S3-SNOW processing software consists of the following components and auxiliary datasets:

- *snap-core* module
- *snap-gpf* module
- *s3tbx-idepix-olci* plugin (specific version for S3-SNOW)
- *s3tbx-olci-o2corr* plugin
- *s3tbx-snow* plugin
- *snap-slope* plugin
- GIMP Digital Elevation Model for Greenland
- lookup tables for OLCI O2 correction


These components are described in more detail in the following subsections.

The Sentinel Application Platform (SNAP)
----------------------------------------

A common architecture for all Sentinel Toolboxes has been jointly developed by Brockmann Consult, Array Systems
Computing and C-S called the Sentinel Application Platform (SNAP).

The SNAP architecture is ideal for Earth Observation processing and analysis due to various technological
innovations as well as approved concepts from the BEAM toolbox. Most of the software components listed above make
use of various SNAP core capabilities.

A good starting point for much more detailed information is the SNAP homepage [`4 <intro.html#References>`_], and also
the comprehensive help documentation integrated in the SNAP desktop application.

The SNAP Graph Processing Framework
-----------------------------------

One of the key components in SNAP is the Graph
Processing Framework (GPF) for creating user-defined processing chains. All provided S3-SNOW processors make use of this
framework.

Within SNAP, the term data processor refers to a software module which creates an output product from one or more
input products configured by a set of processing parameters.
The GPF framework was originally developed for the BEAM toolbox, the precursor of SNAP.
Since the early days of BEAM, a number of data processors have been developed; some of them are standard modules while others
are contributed by 3rd parties. All of these data processors have been developed using a dedicated processing
framework which was already part of the first version of BEAM.

Based on the experience collected within a number of projects, the SNAP authors have developed what is now the
SNAP Graph Processing Framework.
The GPF provides all the features inherited from BEAM, but adds a number of new ones for developers and
reduces the amount of source code to write while drastically improving its readability and maintainability.

Much more detailed information on the SNAP GPF is provided by
the specific GPF help documentation integrated in the SNAP desktop application.

The IdePix OLCI Pixel Classification Processor
----------------------------------------------

IdePix (Identification of Pixels) is a pixel classification tool which has been developed by BC originally for BEAM
and has been used for a variety of projects. It was transferred to SNAP and is continuously being further
developed.

Among the supported sensors is OLCI, which made IdePix the most appropriate candidate for cloud and snow identification in
the S3-SNOW project.

Originally, IdePix has been developed as an internal component of the SNAP Sentinel-3 toolbox. To increase flexibility,
the sub-processors for the various sensors were recently extracted to make them available as separate plugins.
One of these plugins is the IdePix Sentinel-3 OLCI processor.
The processor described here is a special version of this plugin, being adapted for the specific needs for a pixel
classification within S3-SNOW. This allows to more easily provide special user options which are ultimately not
needed in other projects than S3-Snow, and in return leave out other options which are not relevant for S3-Snow.

The IdePix classification algorithm for Sentinel-3 OLCI is based on a neural network approach. A common neural net
is used for both land and water pixels. As input for the neural net, the square roots of the OLCI TOA reflectances
(obtained from an internal radiance-to-reflectance conversion) in all 21 bands are used. As output, the neural net
finally provides per pixel one of the properties 'cloud sure', 'cloud ambiguous', 'cloud'
(which means sure OR ambiguous), or 'snow/ice'.

Although the IdePix classification for OLCI has been tested and successively improved
within various activities, some limitations and weaknesses in cloud detection (most of them well
known from other existing cloud masking approaches) could not be solved to 100%. Among these is the distiction of
cloud and snow/ice, which is very important for the usage for S3-SNOW, and which has shown to be often rather poor.
Therefore, a new approach to detect clouds over snow/ice has been introduced in the IdePix OLCI version for S3-SNOW
which makes use of the O2 correction algorithm provided by R.Preusker (Spectral Earth, Berlin), and which has been
implemented in the OLCI O2 Correction Processor (see next section). As additional output, a binary band 'cloud_over_snow'
is provided.

The pixel classification with IdePix is an optional processing step in S3-SNOW (although recommended in most cases),
applied on the same OLCI L1b products which are being considered for the snow properties retrieval.

The OLCI O2 Correction Processor
--------------------------------

The OLCI O2 Correction Processor provides a 'harmonisation' of O2 bands, which means a modification of the effective
transmittances in O2A bands 13, 14 and 15 to their values which would be measured at their mean wavelengths and with
nominal bandwidth. The corresponding algorithm was provided by R.Preusker (Spectral Earth, Berlin) and is described
in detail in [`2 <intro.html#References>`_]. Among various outputs, the processor provides the rectified and desmiled
transmission for OLCI band 13 (761.25nm) which is used by the IdePix classification for the retrieval of clouds
over snow (previous section).

The OLCI Snow Properties Processor
----------------------------------

The Snow Properties processor is the key component for the processing in S3-SNOW. The processor provides the
implementation
of the algorithms for the various snow properties of interest. These algorithms are also described
in detail in [`2 <intro.html#References>`_].

As input, the processor requires an OLCI L1b product (original or being Rayleigh corrected in a preprocessing step).
Optionally, an IdePix pixel classification product can be provided as additional input. The output is an amount of
snow properties of interest, defined by the user via processing parameters. This is described in detail in
(...s3snow_usage.rst...).


The SNAP Slope Processor
------------------------

The Slope Processor provides pixelwise terrain slope and aspect angle from an arbitrary input product containing
a band with terrain height (i.e. a DEM product). For S3-SNOW this processor is provided as utility tool, as slope
and aspect are often useful information for the validation of snow properties.


The GIMP Digital Elevation Model for Greenland
----------------------------------------------

todo


Lookup Tables
-------------

Various lookup tables are used for the OLCI O2 correction, which in return is part of the IdePix OLCI
pixel classification, all described in more detail in
[`2 <intro.html#References>`_]. These lookup table are not provided separately, bus as an internal part of the
OLCI O2 correction processor plugin.

.. index:: Processing Flow

Processing Flow
===============

IdePix OLCI Pixel Classification Processor
------------------------------------------

todo

OLCI O2 Correction Processor
----------------------------

todo

S3-SNOW Snow Properties Processor
---------------------------------

The overall processing flow of the SNAP TCWV processor is shown in :numref:`tcwv_chain`.

.. _tcwv_chain:
.. figure::  pix/tcwv_chain.png
   :align:   center
   :scale: 80 %

   Processing flow of the SNAP TCWV processor.

As mentioned, L1b products from MERIS or MODIS are used as input. These products are pre-processed with the IdePix
pixel classification module. Idepix provides a classification flag and the reflectance bands (converted from radiances
in case of MERIS) needed for the TCWV retrieval. Further optional input (per pixel) are prior values for temperature,
pressure, wind speed, and an initial TCWV guess. Ideally, these priors are taken from an external data source to provide
values of good quality. For the S3-SNOW TCWV processing on Calvalus, these data were taken from ERA-Interim
[`12 <intro.html#References>`_] products
which were interpolated and collocated onto the initial L1b/IdePix product grid. If no priors are provided, the
processor will use reasonable constant values, but this is not recommended for good TCWV retrievals.

The IdePix products (optionally including the prior bands) are the input for the TCWV processing step, which
provides the final TCWV products (TCWV + flag band).

SNAP Slope Processor
--------------------

The overall processing flow of the SNAP CTP processor is shown in :numref:`ctp_chain`.

.. _ctp_chain:
.. figure::  pix/ctp_chain.png
    :align:   center
    :scale: 80 %

    Processing flow of the SNAP CTP processor.








