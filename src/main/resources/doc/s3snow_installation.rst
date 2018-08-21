.. index:: S3-SNOW Processor Installation
.. _s3snow_installation:

=====================
Software Installation
=====================

.. BC

Overview
========

This chapter describes the overall S3-SNOW software installation procedure.

.. index:: Usage Requirements

Usage Requirements
==================

General Requirements
--------------------

In general, all S3-SNOW processors require SNAP in the latest snapshot version (v7.0.0-SNAPSHOT), which comes with
the S3-SNOW software delivery for all platforms (Windows, Linux, MacOS). [#f1]_

.. rubric:: Footnotes

.. [#f1] This snapshot version can be replaced by the next official release (v7.0.0) when available. This will be announced
         on the SNAP website.


Operating System
----------------

The S3-SNOW software can be run on any operating system which is suppported by SNAP (Windows, Linux, MacOS).

Hardware Requirements
---------------------

The S3-SNOW processing system is a complex piece of software. Although the algorithms for the snow properties
retrieval are mostly relatively simple, the effort for data input/output is fairly high in case of full resolution
OLCI products.
Therefore, up-to-date and  sufficiently powerful and
dimensioned hardware is strongly recommended for reliable and convenient processing.

.. index:: Contents of Software Bundle

Contents of the Processing Software Bundle
==========================================

The S3-SNOW processing software consists of the following components:

- *snap-core* module: part of the SNAP installation package
- *snap-gpf* module:  part of the SNAP installation package
- *s3tbx-idepix-olci*: nbm plugin file
- *s3tbx-olci-o2corr*: nbm plugin file
- *s3tbx-snow*: nbm plugin file
- *snap-slope*: nbm plugin file
- GIMP Digital Elevation Model for Greenland: GeoTIFF auxiliary file

See listing of files in :numref:`sw_bundle`.

.. _sw_bundle:
.. figure::  pix/sw_bundle.png
   :align:   center
   :scale: 80 %

   Contents of the S3-SNOW software bundle.


The processor versions were frequently incremented during the iterative development process.
E.g., the current version of the Snow Properties Processor is v2.0.6.

.. index:: How to get the Software

How to get the Software
=======================

The S3-SNOW processing software bundle can be obtained from the S3-SNOW ftp site hosted at BC with the
following configuration:

- FTP, Port 21
- ftp.brockmann-consult.de
- username: s3snow
- password: $3Sn0W@bc
- subdirectory: software


.. index:: Installation Process

Installation Steps
==================

Installation of the SNAP Software
---------------------------------

Download SNAP (Unix version) from the SNAP web site [`4 <intro.html#References>`_] and follow the
information and instructions for installation given there.

Installation of the S3-SNOW Processor modules
---------------------------------------------

Once SNAP has been installed, the installation of all NBM plugin files needs to be done from the 'Plugins' toolwindow
in the SNAP Desktop application. This is illustrated in the figure sequence :numref:`plugins_in_tools_menu` to
:numref:`add_plugins_confirm_restart`.

.. _plugins_in_tools_menu:
.. figure::  pix/plugins_in_tools_menu.png
   :align:   center
   :scale: 80 %

   The SNAP menu entry for installation of plugins.

.. _add_plugins:
.. figure::  pix/add_plugins.png
   :align:   center
   :scale: 80 %

   Selection of plugins to be installed.

.. _add_plugins_confirm:
.. figure::  pix/add_plugins_confirm.png
   :align:   center
   :scale: 60 %

   Confirmation of selected plugins (step 1 of 4).

.. _add_plugins_confirm_restart:
.. figure::  pix/add_plugins_confirm_restart.png
   :align:   center
   :scale: 60 %

   Final confirmation for restart after selection of plugins.


After restart of SNAP, the installed processors will be available from their dedicated menu entries. This will be
shown in more detail in the next chapter.

