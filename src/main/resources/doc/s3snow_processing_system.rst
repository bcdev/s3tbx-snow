.. role:: underline
    :class: underline

.. index:: SNAP S3-SNOW Processors

.. _s3snow_processing_system:

===========================
The SNAP S3-SNOW Processors
===========================

Overview
========

The key goal of the S3-SNOW project and its successor SICE regarding software development,
production and dissemination was to
implement the proposed algorithms for several  key  snow  parameters from  Sentinel-3 OLCI data in free and easily
accessible open source toolboxes, notably and foremost ESA’s SNAP toolbox.
During the implementation process, all processing software was distributed within the consortium in frequent cycles
for the purpose of a comprehensive validation from an appropriate OLCI database containing products covering a variety of
selected snow-covered areas. All SNAP
S3-SNOW and SICE processors are available as SNAP plugins and can be run within SNAP on any supported platform
(Windows, Linux, MacOS). The :underline:`procedure` for installation and operation is described in this chapter.

.. index:: Theoretical Background Summary

Theoretical Background
======================

The motivation and theoretical background for the retrieval of key snow parameters is summarized in the S3-SNOW project
proposal [`1 <intro.html#References>`_].
The underlying algorithms are described in detail in the corresponding project ATBD  [`2 <intro.html#References>`_].

.. index:: Processing Environment

Processing Environment
======================

As said, the S3-SNOW and SICE processors are available as SNAP plugins and can be run within SNAP on any supported platform
(Windows, Linux, MacOS).
The chapter :doc:`s3snow_installation` describes in more detail how to install the plugins in SNAP.

.. index:: Processing Components

Processing Components
=====================

The SNAP S3-SNOW processing software consists of the following components and auxiliary datasets:

- *SNAP Sentinel-3 toolbox* (current version including latest updates is 7.0.2)
- *snap-slope* processor  (comes with the SNAP Sentinel-3 toolbox)
- *s3tbx-olci-o2corr* processor (comes with the SNAP Sentinel-3 toolbox)
- lookup tables for OLCI O2 harmonisation (come with the SNAP Sentinel-3 toolbox)
- *s3tbx-snow* plugin (current version is 3.0)
- *idepix-core* plugin (current version is 7.0.1)
- *idepix-olci* plugin (current version is 7.0.1)
- GIMP Digital Elevation Model for Greenland

These components are described in more detail in the following subsections.
Note that, compared to the previous version of this SUM, the software package looks a bit different due to
the further SNAP evolution during 2019 towards current version 7.0.2.
I.e., the O2 harmonisation processor (described in more detail later) is now
an internal part of the Sentinel-3 toolbox, whereas the IdePix pixel classification modules are now
provided as separate plugins.

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

The OLCI Snow Properties Processor
----------------------------------

The Snow Properties Processor (SPP) is the key component for the processing in S3-SNOW. The processor provides the
implementation
of the algorithms for the various snow properties of interest. These algorithms are also described
in detail in [`2 <intro.html#References>`_].

As input, the processor requires an OLCI L1b product (original or being Rayleigh corrected in a preprocessing step).
Optionally, an IdePix pixel classification product (see below) can be provided as additional input. The output is a set of
snow properties of interest, defined by the user via processing parameters. This is described in detail in
the chapter :doc:`s3snow_usage`.

The OLCI SICE Snow Properties Processor
----------------------------------

The SICE Snow Properties Processor (SICE SPP) is the most recent processor provided for the retrieval of snow properties.
As it contains various improvements compared to the SPP, this processor is the recommended one for most users.
However, the SPP is still a useful alternative for experienced users as it contains many user options to change
specific algoritnm parameters as well as to generate additional bands in the final snow product.
These underlying algorithms are described in detail in the latest version of [`2 <intro.html#References>`_].

As input, the SICE processor requires both an OLCI L1b product AND a corresponding Rayleigh corrected product from
a preprocessing step.
As for the SPP, an IdePix pixel classification product (see below) can be optionally provided as additional input.
The output is again a set of snow properties of interest, described in detail in the chapter :doc:`s3snow_usage`.

The IdePix OLCI Pixel Classification Processor
----------------------------------------------

IdePix (Identification of Pixels) is a pixel classification tool which has been developed by BC originally for BEAM
and has been used for a variety of projects. It was transferred to SNAP and is continuously being further
developed.

Among the supported sensors is OLCI, which made IdePix the most appropriate candidate for cloud and snow identification in
the S3-SNOW and SICE projects.

Originally, IdePix has been developed as an internal component of the SNAP Sentinel-3 toolbox. To increase flexibility,
the sub-processors for the various sensors were recently extracted to make them available as separate plugins.
One of these plugins is the IdePix Sentinel-3 OLCI processor which can now be used in its standard version
as it has been further improved during 2019 and provides now all the needs for S3-SNOW and SICE,
i.e. the distiction of cloud and snow/ice which now works reasonably well.
(It is no longer necessary to use a 'special version' of Idepix OLCI, as described in previous SUM versions.)

The IdePix classification algorithm for Sentinel-3 OLCI is based on a neural network approach. A common neural net
is used for both land and water pixels. As input for the neural net, the square roots of the OLCI TOA reflectances
(obtained from an internal radiance-to-reflectance conversion) at all 21 wavelengths are used. As output, the neural net
finally provides per pixel one of the properties 'cloud sure', 'cloud ambiguous', 'cloud'
(which means sure OR ambiguous), or 'snow/ice'.

The pixel classification with IdePix is an optional processing step in S3-SNOW as well as in SICE
(although recommended in most cases),
applied on the same OLCI L1b products which are being considered for the snow properties retrieval.

The OLCI O2 Harmonisation Processor
--------------------------------

The OLCI O2 harmonisation Processor provides a 'harmonisation' of O2 wavebands, which means a modification of the effective
transmittances in O2A wavebands 13, 14 and 15 to their values which would be measured at their mean wavelengths and with
nominal bandwidth. The corresponding algorithm was provided by R.Preusker (Spectral Earth, Berlin) and is described
in detail in [`2 <intro.html#References>`_]. Among various outputs, the processor provides the rectified and desmiled
transmission for OLCI waveband 13 (761.25nm) which is used by the IdePix classification for the detection of clouds
over snow (previous subsection).

This processor has now become a part of the current Sentinel-3 toolbox, therefor it is no longer needed to install
it from a separate plugin.

The SNAP Slope Processor
------------------------

The Slope Processor provides pixelwise terrain slope and aspect angle from an arbitrary input product containing
a band with terrain height (i.e. a DEM product). In addition, the variance of elevation over a 3x3 pixel window is
provided. For S3-SNOW this processor is provided as utility tool, as slope
and aspect are often useful information for the validation of snow properties.


The GIMP Digital Elevation Model for Greenland
----------------------------------------------

A Digital Elevation Model for Greenland has been generated within the GIMP project. This product has been post-processed
by BC and is provided in GeoTIFF format with a resolution of ~90m. As only layer in this product, the DEM altitude
given in metres is provided. The altitude is e.g. used as input by the OLCI O2 Harmonisation Processor.
The GIMP DEM product is illustrated in :numref:`gimp_dem`.

.. _gimp_dem:
.. figure::  pix/gimp_dem.png
   :align:   center
   :scale: 80 %
    
   Illustration of the GIMP DEM for Greenland.

Using the SNAP Slope Processor, this product can be used as input to derive the corresponding slope and aspect.


Lookup Tables
-------------

Various lookup tables are used for the OLCI O2 harmonisation, which in return is part of the IdePix OLCI
pixel classification, all described in more detail in
[`2 <intro.html#References>`_]. These lookup tables are not provided separately, but as an internal part of the
OLCI O2 Harmonisation processor.

.. index:: Processing Flow

Processing Flow
===============

The overall processing flow and the interaction of the S3-SNOW components are illustrated in :numref:`processing_flow`.

.. _processing_flow:
.. figure::  pix/processing_flow_2.png
   :align:   center
   :scale: 80 %

   Processing flow of the S3-SNOW processors. See text for details.

The same is illustrated for SICE in :numref:`processing_flow_sice`. The main difference to S3-SNOW is that the
Rayleigh corrected product is needed as mandatory input for the SICE Snow properties processor, thus it needs
to be generated in a pre-processing step.

.. _processing_flow_sice:
.. figure::  pix/processing_flow_sice.png
   :align:   center
   :scale: 80 %

   Processing flow for SICE. See text for details.

The colour and arrow schemes in the diagrams have the following meaning:

- **red** : The standard processing flow for snow properties retrieval. The red boxes indicate the mandatory input
  products and processing modules: An OLCI L1b radiances product is used as input product for the SPP.
  If not provided as pre-processed product, BRRs are computed from an internal call of the SNAP Rayleigh Correction
  Processor, which in return are used for the retrieval of the various snow properties. In opposite to SPP, SICE needs
  the BRR product as mandatory input from pre-processing.
- **orange** : Alternative processing flow in SPP for snow properties retrieval:
  An OLCI BRR product is used as input product
  for the SPP. This BRR product has been computed independently in a preprocessing step, directly
  using the Rayleigh Correction Processor.
- **green** : Optional processing, i.e. cloud classification: An OLCI L1b radiances product is used as input product
  for the IdePix Pixel Classification Processor. The IdePix output product can then be used as optional second input
  product for the SPP or SICE. Internally, IdePix calls the O2 Harmonisation Processor to obtain the
  O2 waveband transmissions being used to generate the improved cloud classification band 'cloud_over_snow'. An optional
  DEM product can be used as input for the O2 Harmonisation Processor. If no DEM is specified by the user, the altitude band
  from the Olci L1b product is used.
- **grey** : Additional processing options, not directly used in the snow properties retrieval. I.e., O2 harmonisation
  and slope/aspect computation, as outlined above.
- **solid arrows** : indicate input/output to/from a processing module
- **dashed arrows** : indicate internal calls of one processing module into another









