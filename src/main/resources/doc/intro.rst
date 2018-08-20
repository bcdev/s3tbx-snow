.. _intro:

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

Purpose and Scope
=================

This document is the User Manual for the SNAP processors written in 
`Java <http://www.oracle.com/java>`_ which have been developed in the frame of the S3-SNOW
project. Its purpose is to describe in detail how to obtain, install and operate these processors. Also, a
comprehensive overview of the related data products (input as well as intermediate and final products) is provided.

The explicit structure of the document is as follows:

* Chapter 1 is this introduction.
* `Chapter 2 <s3snow_processing_system.html>`_ gives an overview of the SNAP S3-SNOW processors.
* `Chapter 3 <s3snow_products.html>`_ describes all relevant SNAP S3-SNOW products.
* `Chapter 4 <s3snow_installation.html>`_ explains how to get and install the processing software.
* `Chapter 5 <s3snow_usage.html>`_ explains how to run the processing software.
* `The Annex <annex.html>`_ contains various annexes.

References
==========

1.  Scientific Exploitation of Operational Missions (SEOM) Sentinel-3 for Science, Land Study 1: SNOW:
    Technical, Management/Administrative, Implementation, Financial and Contractual Proposal.
    Issue 1, Revision 1, 16 March 2016.

2.  Kokhanovsky, A., Box, J.E., Lamare, M., Dumont, M., Picard, G., Danne, O., and C. Brockmann:
    Algorithm Theoretical Basis Document: Snow Properties Retrieval from Sentinel-3. Version 3.0, 30 April 2018.

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
BRR                      Bottom-of-Rayleigh Reflectance
-----------------------  ---------------------------------------------------------------------------------------------
CCI                      Climate Change Initiative
-----------------------  ---------------------------------------------------------------------------------------------
DEM                      Digital Elevation Model
-----------------------  ---------------------------------------------------------------------------------------------
ESA                      European Space Agency
-----------------------  ---------------------------------------------------------------------------------------------
GIMP                     Greenland Ice Mapping Project
-----------------------  ---------------------------------------------------------------------------------------------
GPF                      Graph Processing Framework
-----------------------  ---------------------------------------------------------------------------------------------
IdePix                   Identification of Pixels
-----------------------  ---------------------------------------------------------------------------------------------
JDK                      Java Development Kit
-----------------------  ---------------------------------------------------------------------------------------------
MERIS                    Medium Resolution Imaging Spectrometer
-----------------------  ---------------------------------------------------------------------------------------------
nbm                      Netbeans module
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
SNAP                     Sentinel Application Platform
-----------------------  ---------------------------------------------------------------------------------------------
TOA                      Top Of Atmosphere
=======================  =============================================================================================

