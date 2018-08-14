.. _intro:

============
Introduction
============

Project background
==================

The SEOM S3 ‘advanced Clouds, Aerosols and WAter vapour products for Sentinel-3/OLCI’ CAWA project aims to the
development and improvement of advanced atmospheric retrieval algorithms for the Envisat/MERIS and Sentinel-3/OLCI
mission. A sensor comprehensive and consistent 1D-Var water vapour algorithm has been developed and applied to the MERIS,
MODIS and first available OLCI measurements. An innovative and consistent cloud top pressure 1D-Var procedure was defined
for MERIS and all three OLCI O2 A-band channels, which will significantly improve the retrieval accuracy. The
challenging and innovative GRASP algorithm for the retrieval of aerosols and surface properties has already shown
its advantage in comparison to conventional aerosol retrieval methods. All three algorithms will be further improved,
applied to the complete MERIS dataset, to a four months MODIS global time series and six months of OLCI data. We expect
to create improved consistent datasets of water vapour, cloud properties, namely cloud top pressure, and aerosol and
surface pressure. The intention of the CAWA team is to establish new and improved procedures to estimate atmospheric
properties, which also improve the retrieval of land and ocean properties.


Purpose and Scope
=================

This document is the User Manual for the SNAP TCWV and CTP processors written in `Python <http://www.python.org>`_ and
`Java <http://www.oracle.com/java>`_ which have been developed in the frame of the CAWA
project. Its purpose is to describe in detail how to obtain, install and operate these processors. Also, a
comprehensive overview of all related data products (input as well as intermediate and final products) is provided.

The explicit structure of the document is as follows:

* Chapter 1 is this introduction.
* `Chapter 2 <cawa_processing_system.html>`_ gives an overview of the SNAP CAWA TCWV and CTP processing system.
* `Chapter 3 <cawa_products.html>`_ describes all relevant SNAP CAWA products.
* `Chapter 4 <cawa_installation.html>`_ explains how to get and install the processing software.
* `Chapter 5 <cawa_usage.html>`_ explains how to run the processing software.
* `The Annex <annex.html>`_ contains various annexes.

References
==========

1.  ADVANCED CLOUDS, AEROSOLS AND WATER VAPOUR PRODUCTS FOR SENTINEL-3/OLCI: Technical, Management and
    Financial Proposal. Issue 1.0, 28.03.2014.

2.  Retrieval for Total Coulumn Water Vapor from MERIS/OLCI and MODIS for Land- and Ocean Surfaces.
    CAWA TCWV ATBD,
    available at: https://earth.esa.int/web/sppa/activities/cawa/projects-documents

3.  Retrieval of Cloud Top Pressure from MERIS and  OLCI O2 A-Band Measurements. CAWA CTP ATBD,
    available at: https://earth.esa.int/web/sppa/activities/cawa/projects-documents

4.  The Sentinel Application Platform (SNAP) Web Site,
    available at: http://step.esa.int/main/toolboxes/snap/

5.  Configure Python to use the SNAP-Python (snappy) interface,
    available at: https://senbox.atlassian.net/wiki/display/SNAP/Configure+Python+to+use+the+SNAP-Python+%28snappy%29+interface

6.  CoastColour Project Web Site,
    available at: http://www.coastcolour.org

7.  OceanColour Project Web Site,
    available at: http://www.esa-oceancolour-cci.org

8.  Bourg, L. (2009): MERIS Level 2 Detailed Processing Model. ACRI-ST, Document No. PO-TN-MEL-GS-0006, 15 July 2009.

9.  GlobAlbedo Project Web Site,
    available at: http://globalbedo.org

10. LandCover Project Web Site,
    available at: http://www.esa-landcover-cci.org

11. GlobAlbedo ATBD 'Pixel Classification'. Version 4.1, 26 June 2012.

12. ERA-Interim global atmospheric reanalysis dataset,
    available at: http://www.ecmwf.int/en/research/climate-reanalysis/era-interim

13. European Space Agency: Meris Product Handbook, Issue 3.0, 1 August 2011.

14. MODIS Level 1B Product User’s Guide. For Level 1B Version 6.1.0 (Terra) and Version 6.1.1 (Aqua).
    MODIS Characterization Support Team, Document PUB-01-U-0202- REV C, February 27, 2009.

15. Climate Data Operators (CDO) Web Site,
    available at: https://code.zmaw.de/projects/cdo

16. The Python Download Web Site,
    available at: https://www.python.org/downloads/

17. SNAP Wiki: Configure Python to use the SNAP-Python (snappy) interface,
    available at: https://senbox.atlassian.net/wiki/display/SNAP/Configure+Python+to+use+the+SNAP-Python+%28snappy%29+interface

Acronyms and Abbreviations
==========================

=======================  =============================================================================================
**Acronym**              **Definition**
=======================  =============================================================================================
ATBD                     Algorithm Theoretical Basis Document
-----------------------  ---------------------------------------------------------------------------------------------
BC                       Brockmann Consult
-----------------------  ---------------------------------------------------------------------------------------------
BEAM                     Basic ERS & Envisat (A)ATSR and Meris Toolbox
-----------------------  ---------------------------------------------------------------------------------------------
Calvalus                 CAL/VAL and User Services
-----------------------  ---------------------------------------------------------------------------------------------
CAWA                     Advanced Clouds, Aerosols and WAter vapour products
-----------------------  ---------------------------------------------------------------------------------------------
CCI                      Climate Change Initiative
-----------------------  ---------------------------------------------------------------------------------------------
CTP                      Cloud Top Pressure
-----------------------  ---------------------------------------------------------------------------------------------
DFS                      Distributed File System
-----------------------  ---------------------------------------------------------------------------------------------
DUE                      Data User Element
-----------------------  ---------------------------------------------------------------------------------------------
ESA                      European Space Agency
-----------------------  ---------------------------------------------------------------------------------------------
GPF                      Graph Processing Framework
-----------------------  ---------------------------------------------------------------------------------------------
GRASP                    Generalized Retrieval of Aerosol and Surface Properties
-----------------------  ---------------------------------------------------------------------------------------------
JDK                      Java Development Kit
-----------------------  ---------------------------------------------------------------------------------------------
MERIS                    Medium Resolution Imaging Spectrometer
-----------------------  ---------------------------------------------------------------------------------------------
MODIS                    Moderate Resolution Imaging Spectroradiometer
-----------------------  ---------------------------------------------------------------------------------------------
MR                       Map-Reduce
-----------------------  ---------------------------------------------------------------------------------------------
NetCDF                   Network Common Data Form
-----------------------  ---------------------------------------------------------------------------------------------
OLCI                     Ocean and Land Colour Instrument
-----------------------  ---------------------------------------------------------------------------------------------
RAM                      Random Access Memory
-----------------------  ---------------------------------------------------------------------------------------------
SEOM                     Scientific Exploitation of Operational Missions
-----------------------  ---------------------------------------------------------------------------------------------
SNAP                     Sentinel Application Platform
-----------------------  ---------------------------------------------------------------------------------------------
SNAPPY                   SNAP Python
-----------------------  ---------------------------------------------------------------------------------------------
TCWV                     Total Column of Water Vapour
-----------------------  ---------------------------------------------------------------------------------------------
TOA                      Top Of Atmosphere
-----------------------  ---------------------------------------------------------------------------------------------
UCAR                     University Corporation for Atmospheric Research
=======================  =============================================================================================

