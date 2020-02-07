.. _intro:

.. raw:: latex

   \listoffigures
   \listoftables

============
Introduction
============

Project background
==================

The SEOM Sentinel-3 for Science, Land Study 1: SNOW (referred as 'S3-SNOW' in the following) is “to develop, implement  
and  validate  algorithms  
for deriving  several  key  snow  parameters from  Sentinel  3  optical  satellite  data,  appropriate  for  addressing  
ESA’s  Cryosphere  challenge  C3 : "Seasonal  snow,  lake/river  ice  and  land  ice,  their  effect  on  the
climate  system,  water  resources, energy  and  carbon  cycles:  the  representa<on  of  terrestrial  cryosphere  
in  land  surface,  atmosphere  and climate  models”. See [1] for details.

As a successor of S3-SNOW, a new one-year project named SICE (Sentinel-3 snow and ice products) was established in 2019,
with the particular goal to determine a dry/wet snow and clean/polluted bare ice spectral and broadband optical albedo
1 km daily product for land ice (glaciers, ice caps, ice sheet). Within SICE, both algorithms and corresponding
data processors and software were updated and further improved. The relevant changes of the processors and software
will be described in this updated version of the SUM.

Purpose and Scope
=================

This document is the User Manual for the SNAP processors written in 
`Java <http://www.oracle.com/java>`_ which have been developed in the frame of the S3-SNOW
project and its successor SICE. Its purpose is to describe in detail how to obtain, install and operate these processors.
Also, a comprehensive overview of the related data products (input as well as intermediate and final products) is provided.

The explicit structure of the document is as follows:

* Chapter 1 is this introduction.
* `Chapter 2 <s3snow_processing_system.html>`_ gives an overview of the SNAP S3-SNOW and SICE processors.
* `Chapter 3 <s3snow_products.html>`_ describes all relevant S3-SNOW and SICE products.
* `Chapter 4 <s3snow_installation.html>`_ explains how to get and install the processing software.
* `Chapter 5 <s3snow_usage.html>`_ explains how to run the processing software.

References
==========

1.  Scientific Exploitation of Operational Missions (SEOM) Sentinel-3 for Science, Land Study 1: SNOW:
    Technical, Management/Administrative, Implementation, Financial and Contractual Proposal.
    Issue 1, Revision 1, 16 March 2016.

2.  Kokhanovsky, A., Box, J.E., Lamare, M., Dumont, M., Picard, G., Danne, O., and C. Brockmann:
    Algorithm Theoretical Basis Document: Snow Properties Retrieval from Sentinel-3. Version 2.2, 8 November 2018.

3.  Sentinel-3 OLCI User Guide,
    available at: https://sentinel.esa.int/web/sentinel/user-guides/sentinel-3-olci

4.  S3-SNOW Project Web Site,
    available at: http://seom.esa.int/page_project032.php

5.  The Sentinel Application Platform (SNAP) Web Site,
    available at: http://step.esa.int/main/toolboxes/snap/

6.  CoastColour Project Web Site,
    available at: http://www.coastcolour.org

7.  OceanColour Project Web Site,
    available at: http://www.esa-oceancolour-cci.org

8.  Bourg, L. (2009): MERIS Level 2 Detailed Processing Model. ACRI-ST, Document No. PO-TN-MEL-GS-0006, 15 July 2009.

9.  Sentinel-3A Product Notice: OLCI Level-2 Ocean Colour,
    available at: https://www.eumetsat.int/website/wcm/idc/idcplg?IdcService=GET_FILE&dDocName=PDF_S3A_PN_OLCI_L2_REP&RevisionSelectionMethod=LatestReleased&Rendition=Web

10. Pre-operational Sentinel-3 snow and ice products (SICE) web site.
    available at: http://snow.geus.dk

11. Kokhanovsky, A., Lamare, M., Di Mauro, B., Picard, G., Arnaud, L., Dumont, M., Tuzet, F., Brockmann, C., and J.E. Box:
    On the reflectance spectroscopy of snow. The Cryosphere, 12, 2371-2382, 20 July 2018


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
BOA                      Bottom-of-atmosphere
-----------------------  ---------------------------------------------------------------------------------------------
BOAR                     Bottom-of-atmosphere reflectance
-----------------------  ---------------------------------------------------------------------------------------------
BRR                      Bottom-of-atmosphere Rayleigh Reflectance
-----------------------  ---------------------------------------------------------------------------------------------
CCI                      Climate Change Initiative
-----------------------  ---------------------------------------------------------------------------------------------
CLI                      Command Line Interface
-----------------------  ---------------------------------------------------------------------------------------------
DEM                      Digital Elevation Model
-----------------------  ---------------------------------------------------------------------------------------------
ESA                      European Space Agency
-----------------------  ---------------------------------------------------------------------------------------------
GIMP                     Greenland Ice Mapping Project
-----------------------  ---------------------------------------------------------------------------------------------
GPF                      Graph Processing Framework
-----------------------  ---------------------------------------------------------------------------------------------
GUI                      Graphical User Interface
-----------------------  ---------------------------------------------------------------------------------------------
IdePix                   Identification of Pixels
-----------------------  ---------------------------------------------------------------------------------------------
JDK                      Java Development Kit
-----------------------  ---------------------------------------------------------------------------------------------
MERIS                    Medium Resolution Imaging Spectrometer
-----------------------  ---------------------------------------------------------------------------------------------
NBM                      NetBeans Module
-----------------------  ---------------------------------------------------------------------------------------------
NDSI                     Normalised Difference Snow Index
-----------------------  ---------------------------------------------------------------------------------------------
NetCDF                   Network Common Data Form
-----------------------  ---------------------------------------------------------------------------------------------
OLCI                     Ocean and Land Colour Instrument
-----------------------  ---------------------------------------------------------------------------------------------
PPA                      Probability of Photon Absorption
-----------------------  ---------------------------------------------------------------------------------------------
SEOM                     Scientific Exploitation of Operational Missions
-----------------------  ---------------------------------------------------------------------------------------------
SICE                     Sentinel-3 Snow and Ice
-----------------------  ---------------------------------------------------------------------------------------------
SNAP                     Sentinel Application Platform
-----------------------  ---------------------------------------------------------------------------------------------
SPP                      Snow Properties Processor
-----------------------  ---------------------------------------------------------------------------------------------
SUM                      Software User Manual
-----------------------  ---------------------------------------------------------------------------------------------
TOA                      Top Of Atmosphere
=======================  =============================================================================================

