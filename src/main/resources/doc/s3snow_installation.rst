.. index:: S3-SNOW Processor Installation
.. _s3snow_installation:

================================
Processing Software Installation
================================

.. BC

Overview
========

This chapter describes the overall software installation procedure (processing modules and auxiliary data sets) as well
as the system requirements (hardware and software).

.. index:: Usage Requirements

Usage Requirements
==================

General Requirements
--------------------

In general, the S3-SNOW processors require:

- a 1.8 version of the Java Development Kit (JDK)
- Python v2.7
- SNAP latest release (currently v5.0.0), including IdePix


Operating System
----------------

The software has been developed and tested on Virtual machines based on Linux Ubuntu, which is also used on the
Calvalus system. Also, the required FORTRAN shared libraries included in the software bundle were pre-compiled in a
Linux environment. Therefore, the S3-SNOW software can be run on Linux systems only.

Hardware Requirements
---------------------

The S3-SNOW TCWV and CTP Processing System is a complex piece of software which is based on
massive numerical operations and lookup table access. Therefore the system requires sufficiently powerful and sufficiently
dimensioned hardware for reliable processing. The recommended key parameters for the
hardware are:

- Multi-kernel CPU (8 or more), > 3 GHz
- RAM 8GB or more
- sufficient disk space according to sizes of products


.. index:: Contents of Software Bundle

Contents of the Processing Software Bundle
==========================================

The SNAP TCWV and CTP processing software bundle contains the following components in two separate jar files:

- *snap-s3snow* plug-in (Python files, lookup tables, FORTRAN shared libraries, GPF configuration files)
- *snap-s3snow-io* plug-in

The current processor version is v1.2, therefore we have the two files:

- snap-s3snow-1.2.jar
- snap-s3snow-io-1.2.jar

.. index:: How to get the Software

How to get the Software
=======================

The SNAP TCWV and CTP processing software bundle can be obtained from the S3-SNOW ftp site hosted at BC with the
following configuration:

- SFTP, Port 22
- ftp.brockmann-consult.de
- username: s3snow
- password: 7t86.8K9i7z
- subdirectory: s3snow_processor


.. index:: Installation Process

Installation Steps
==================

Installation of the SNAP Software
---------------------------------

Download SNAP (Unix version) from the SNAP web site [`4 <intro.html#References>`_] and follow the
information and instructions for installation given there.

Installation of the Python Software
-----------------------------------

Download Python v2.7 from the Python web site [`16 <intro.html#References>`_] and follow the
information and instructions for installation given there.

Python Configuration
--------------------

Once downloaded and installed, Python needs to be configured to use the SNAP-Python (snappy) interface.
Instructions for this step are given in detail in the SNAP Wiki [`17 <intro.html#References>`_].

Installation of the S3-SNOW Processor modules
---------------------------------------------

The SNAP TCWV and CTP processor modules need to be installed as follows:

- download the *snap-s3snow-1.2.jar* and *snap-s3snow-io-1.2.jar* into an arbitrary directory
- copy the file *snap-s3snow-io-1.2.jar* to $SNAP_INSTALL_DIR/snap/modules
- unpack the *snap-s3snow-1.2.jar* into an arbitrary *snap-s3snow* directory, e.g. */home/snap-s3snow*
- now all required resources should be in another subdirectory */home/snap-s3snow/resources_bundle*
- To link the snap-s3snow directory to SNAP, edit the file $SNAP_INSTALL_DIR/etc/snap.properties
- In this file, at the end of the file, add the line:
        snap.pythonExtraPaths = *snap-s3snow directory*
  e.g.
        snap.pythonExtraPaths = /home/snap-s3snow

