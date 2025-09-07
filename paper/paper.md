---
title: A Python API for OTTER
authors:
  - name: Noah Franz
    orcid: 0000-0003-4537-3575
    affiliation: 1
  - name: Kate D Alexander
    orcid: 0000-0002-8297-2473
	affiliation: 1
  - name: Sebastian Gomez
    orcid: 0000-0001-6395-6702
	affiliation: 2
affiliations:
  - name: Department of Astronomy and Steward Observatory, University of Arizona, 933 North Cherry Avenue, Tucson, AZ 85721-0065, USA
    index: 1
  - name: Department of Astronomy, The University of Texas at Austin, 2515 Speedway, Stop C1400, Austin, TX 78712, USA
    index: 2

date: 06 September 2025
bibliography: paper.bib

aas-doi: update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal
---

# Summary
The Open mulTiwavelength Transient Event Repository (or "OTTER") is a new catalog of
published transient metadata and photometry. Here we present a Python wrapper on the
RESTful application programming interface (API) built-in to the OTTER backend database
(hereafter the "OTTER API"). Since the OTTER backend is built on the document database
ArangoDB, using the RESTful API directly requires learning the Arango Query Language
(AQL). For python-savvy astronomers this may be a roadblock to programmatic access to
OTTER. Therefore, the primary goal of the OTTER API is to make programmatic access easy
and fast for astronomers. In addition to wrapping the RESTful API, the OTTER Python API
also provides additional methods for 1) Converting the stored photometry in standard
units; 2) Helper methods for querying additional astronomy database servicesl and 3)
Methods for quickly plotting the photometry stored in OTTER.

# Statement of need
Transient astrophysical events provide a unique view of high energy phenomena evolving
on human timescales. Examples include supernovae, gamma ray bursts, tidal disruption
events, and many more weird ones that do not fit into these typical classifications
[@1966ApJ...143..626C; @1995ApJS..101..181W; @smartt_progenitors_2009;
@1966ApJ...143..626C; @maoz_observational_2014; @kouveliotou_identification_1993;
@norris_frequency_1984; @2003ApJ...591..288H; @eichler_nucleosynthesis_1989;
@narayan_gamma-ray_1992; @hills_possible_1975; @rees_tidal_1988].
Understanding the many astrophysical transients observed is not easy
and typically requires detailed multiwavelength observations and analyses of the
population of events [@costa_discovery_1997; @olivares_e_multiwavelength_2015;
@gezari_x-ray_2017; @pasham_discovery_2018; @eftekhari_radio_2018;
@eftekhari_late-time_2021; @laskar_first_2022; @laskar_radio_2023;
@margutti_luminous_2023; @gomez_type_2024; @christy_peculiar_2024;
@hajela_eight_2024; @guolo_systematic_2024; @masterson_new_2024;
@alexander_multi-wavelength_2025]. However, it can take years to gather the multiwavelength datasets
for a population of transients from the literature.

Furthermore, transient classification is a non-trivial process that typically requires
detailed spectroscopic and multiwavelength follow-up observations [e.g., @arcavi_continuum_2014; @charalampopoulos_detailed_2022]. However, with the
advent of Rubin Observatory's Legacy Survey of Space and Time [@ivezic_lsst_2019] the number of known
transients will increase by at least an order of magnitude [@van_velzen_optical_2011;@bricman_prospects_2020] necessitating additional
methods for classifying transients, such as machine learning [@villar_supernova_2019;
@villar_superraenn_2020;@villar_deep-learning_2021;@gomez_fleet_2020;@gomez_identifying_2023;
@stein_texttttdescore_2023; @de_soto_superphot_2024; @boesky_splash_2025]. However, machine learning
classifiers require large training datasets that can, similar to above, take years to
gather and clean into a coherent training set. This makes services for archiving and
cataloging transient metadata and photometry to maximize the scientific output of
the Rubin Observatory time domain survey.

For both of these reasons we created the Open mulTiwavelength Transient Event Repository
[OTTER; @franz_otter_2025], a scalable catalog of transient events specifically designed to store the
nuances of the multiwavelength photometry available on these events. OTTER is a
successor to the Open Astronomy Catalogs, which have not been maintained since 2022, but
with a focus on multiwavelength datasets [@guillochon_open_2017; @auchettl_new_2017].
To store the various nuances of transient event photometry at multiple wavelengths (e.g., the model used to reduce and extract
a flux from an X-ray observation), we chose to use a flexible document database as our
backend, called ArangoDB. This also enables us to store both default and non-default
values or measurements of different metadata properties associated with the transient.

One of our primary goals of OTTER is ease of access to the dataset, including a way to
programmatically access it to make the curation of large transient samples easier.
ArangoDB has built in RESTful API endpoints for progrommatic access to the
data. However, the API endpoints expect queries in the syntax of the "Arango Query
Language" (AQL). This would provide a barrier to use of the detailed dataset available
for busy astronomers who do not have the time to learn yet another query language (on
top of the existing SQL and ADQL).

To help overcome this barrier, we present an object-oriented Python API for access to
the OTTER dataset. This API acts as a "wrapper" on the RESTful API with many additions
that make it even easier to access and analyze the dataset. Some of the more useful
components include the following.
* In OTTER we store photometry in a "raw" form, or essentially as close to the actual
  publication value as possible. This makes the data more reproducible since it does
  not rely on complex processing scripts that may have various bugs. However, this also
  means that the data is not stored in consistent units (but the unit of the photometry
  point is stored). In the OTTER API we automatically convert the photometry for the
  user into a consistent, requested set of units. This is specifically done in the
  `Otter.get_phot` and `Transient.clean_photometry` methods. Behind the scenes this
  conversion is done using the python package synphot [@stsci_development_team_synphot_2018].
* The raw data from some astronomical observations is made public and reduced by
  multiple people. Depending on the differences in the reduction methodology this may
  produce slightly different flux measurements. If this is the case, we store both flux
  measurements in the OTTER database. However, to help users de-duplicate these
  datasets we provide a method for finding and choosing only one of the multiple
  reductions. This is done in the `Transient.deduplicate_photometry` method.
* Sometimes users just want to quickly view the photometry for a specific transient
  event as either a light curve (flux vs. time) or a spectral energy distribution (flux
  vs. wavelength, frequency, or energy; i.e., an "SED"). In the `plotter` module of the OTTER API we
  provide numerous methods for quickly and automatically plotting the photometry. The
  most powerful of these methods (albeit the least flexible) is the `query_quick_view`
  method which passes the users query to the OTTER API and returns both a light curve
  and SED of the associated photometry [@hunter_matplotlib_2007; @plotly].
* Identifying the host galaxy of a transient event is not necessarily trivial and it is
  important to store the information about the host galaxy (name and coordinates) in
  databases like OTTER. However, there are numerous existing astronomical databases that store galaxy properties and
  we do not want to duplicate their efforts. We therefore provide methods as part of
  the `Host` object to query other public services for e.g., host photometry or spectra.
  These other services include Simbad [@2000A&AS..143....9W], ATLAS [@ATLAS], ZTF
  [@ZTF], iPTF [@iPTF], ASAS-SN [@ASASSN,@2017PASP..129j4502K,@2023arXiv230403791H],
  Vizier [@10.26093/cds/vizier, @vizier2000], WISE
  [@WISE,@NEOWISE,@NEOWISE_Reactivation,@2020MNRAS.493.2271H],
  FIRST [@1997ApJ...475..479W], NVSS [@1998AJ....115.1693C], HEASARC, and Sparcl
  [@juneau_sparcl_2024] --- most of which are available through astroquery [@ginsburg_astroquery_2019].
* As mentioned previously, we allow for the storage of different measurements
  associated with the same value. For example, we allow multiple redshift measurements
  for a single transient. The OTTER API will automatically choose a default value if
  multiple values are present for a single property.

# Conclusions and other considerations

Finally, it should be noted that we also developed a front end web application for
access to the database. This web application uses the API described here for interfacing
with the database. In addition to simply providing a GUI for interaction with the dataset,
the web application adds other useful functionality like uploading new data after a user
publish it. The upload functionality is solely in the web application to allow for
additional vetting of the uploaded datasets.

We hope that OTTER, in its entirety, will be useful infrastructure tool for
time domain science moving forward. Even more, we hope that the OTTER API described here will
make access to that dataset easier for users, without the afformentioned roadblocks. We welcome
GitHub issues with comments and feedback (or even pull requests!) from the community on our
GitHub repository.

# Citations
