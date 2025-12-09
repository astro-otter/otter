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
published transient data. Here we present a thick Python wrapper on the
REpresentational State Transfer (REST) application programming interface (API) built-in to the OTTER backend database
(the "OTTER API"). Since the OTTER backend is built on the document database
ArangoDB, using the REST API directly requires learning the Arango Query Language
(AQL). Since AQL has a niche user base, OTTER users unfamiliar with it may face a roadblock to programmatic access of
OTTER. To overcome this barrier, we created the OTTER Python API to make programmatic access easy
and fast. In addition to wrapping the REST API, the OTTER Python API
also provides additional methods for 1) Converting the stored photometry to standard
units; 2) Helper methods for querying additional astronomy database services; and 3)
Methods for quickly plotting the photometry stored in OTTER.

# Statement of need
Transient astrophysical events provide a unique high energy laboratory that evolves
on human timescales. Examples include supernovae, gamma ray bursts, tidal disruption
events, and many other exotic transients
[@{1966ApJ...143..626C}; @{1995ApJS..101..181W}; @smartt_progenitors_2009;
 @maoz_observational_2014; @kouveliotou_identification_1993;
@norris_frequency_1984; @{2003ApJ...591..288H}; @eichler_nucleosynthesis_1989;
@narayan_gamma-ray_1992; @hills_possible_1975; @rees_tidal_1988].
Developing an understanding of the many astrophysical transients observed is not easy
and typically requires detailed multiwavelength observations and analyses of the
population of events [@costa_discovery_1997; @olivares_e_multiwavelength_2015;
@gezari_x-ray_2017; @pasham_discovery_2018; @eftekhari_radio_2018;
@eftekhari_late-time_2021; @laskar_first_2022; @laskar_radio_2023;
@margutti_luminous_2023; @gomez_type_2024; @christy_peculiar_2024;
@hajela_eight_2024; @guolo_systematic_2024; @masterson_new_2024;
@alexander_multi-wavelength_2025]. However, it can take years to gather the multiwavelength datasets
from the literature, necessitating publicly available transient event data archives.

Furthermore, transient classification is a non-trivial process that typically requires
detailed spectroscopic and/or multiwavelength follow-up observations [e.g., @arcavi_continuum_2014; @charalampopoulos_detailed_2022],
a method that is only feasible with the current transient discovery rates. With the
advent of Rubin Observatory's Legacy Survey of Space and Time [@ivezic_lsst_2019] the number of known
transients will increase by at least an order of magnitude [@van_velzen_optical_2011;@bricman_prospects_2020].
Therefore, additional methods for classifying transients, such as machine learning, are required [@villar_supernova_2019;
@villar_superraenn_2020;@villar_deep-learning_2021;@gomez_fleet_2020;@gomez_identifying_2023;
@stein_texttttdescore_2023; @de_soto_superphot_2024; @boesky_splash_2025]. However, machine learning
classifiers require large training datasets that can be laborious to
curate. This further motivates archival services for cataloging transient metadata and photometry, and
will be necessary to maximize the scientific output of the Rubin time domain survey.

For both of these reasons we created the Open mulTiwavelength Transient Event Repository
[OTTER, @franz_otter_2025], a scalable catalog of transient event metadata and photometry. OTTER is a
successor to the Open Astronomy Catalogs [OAC, @auchettl_new_2017; @guillochon_open_2017] [^1],
but designed and optimized for multiwavelength datasets.
To store the various nuances of multiwavelength photometry (e.g., the model used to reduce and extract
a flux from an X-ray observation), we chose to use a flexible document database management system as our
backend: ArangoDB. The nested structure of the document database files also provides an intuitive way to
store multiple values of a single measurement when different sources disagree.

One of our primary goals of OTTER is ease of access to the dataset, including a way to
programmatically access it to make the curation of large transient samples easier.
ArangoDB has a built-in REST API for programmatic access to the
data. However, the API endpoints expect queries in the syntax of the "Arango Query
Language" (AQL). Learning a new query language creates a barrier for user programmatic
access to the indispensable dataset available in the OTTER catalog.

To help overcome this barrier, we present a Python API for access to
the OTTER dataset. This API acts as a thick wrapper on the AQL-based API, with many additions
that make it easier to access and analyze the dataset. Some of these features include:

* In OTTER we store photometry as close to the actual
  published value as possible to make the data more reproducible. However, this also
  means that the data is not stored in consistent units (but the unit of the photometry
  point is stored). In the OTTER API we automatically convert the photometry into the user-requested units.
  Specifically, the conversion is done in the
  `Otter.get_phot` and `Transient.clean_photometry` methods which use
  `astropy` [@astropy_collaboration_astropy_2013; @astropy_collaboration_astropy_2018;
  @astropy_collaboration_astropy_2022] and `synphot` [@stsci_development_team_synphot_2018].
* The same raw data from an astronomical observation may be reduced[^2] by
  multiple, distinct, teams. Depending on the differences in the reduction methodology this may
  produce different flux measurements. If this is the case, we store both flux
  measurements in the OTTER database to allow the user to choose their preferred reduction.
  However, to help users de-duplicate these datasets while curating large samples, we provide
  an (optional) automated algorithm for finding duplicates and choosing only one of the multiple
  reductions. This is done in the `Transient.deduplicate_photometry` method.
* Sometimes users want to quickly view the photometry for a specific transient
  event as either a light curve (flux as a function of time) or a spectral energy distribution (flux
  as a function of wavelength, frequency, or energy; i.e., an "SED"). In the `plotter` module of the OTTER API we
  provide numerous methods for quickly and automatically plotting the photometry [@hunter_matplotlib_2007; @plotly].
* Identifying the host galaxy of a transient event can be difficult and it is
  important to store identifying information (e.g., name and coordinates) for a host galaxy, if it is known.
  However, there are numerous existing astronomical databases that store galaxy properties and
  we do not want to duplicate their efforts. We therefore provide methods as part of
  the `Host` object to query other public services for host photometry or spectra.
  These other services include Simbad [@{2000A&AS..143....9W}], ATLAS [@ATLAS], ZTF
  [@ZTF], iPTF [@iPTF], ASAS-SN [@ASASSN;@{2017PASP..129j4502K};@2023arXiv230403791H],
  Vizier [@vizier2000], WISE [@WISE;@NEOWISE;@NEOWISE_Reactivation;@2020MNRAS.493.2271H],
  FIRST [@{1997ApJ...475..479W}], NVSS [@{1998AJ....115.1693C}], HEASARC, and Sparcl
  [@juneau_sparcl_2024] --- most of which are queried using the astropy-affiliated `astroquery` package [@ginsburg_astroquery_2019].
* Users may want to compare new observations stored locally with the publicly available data in OTTER.
  As part of the OTTER API we make this very easy as long as their data is stored in a well-documented
  CSV file format (see the OTTER web application upload form or the example jupyter notebook titled
  "Interfacing with Private Data"). When the data is stored like this
  a user is able to use the `Otter.from_csvs` method to construct an `Otter` object that will pass
  their queries to both the public OTTER dataset and the one locally stored and return all relevant
  information in a consistent format.
* We allow for the storage of different measurements (e.g., redshift, discovery date, etc.)
  associated with the same property of the transient. The OTTER API will automatically choose a default value if
  multiple measurements are present for a single property.

[^1]: The OAC was an indispensable resource but has not been maintained since 2022, further necessitating a successor like OTTER.
[^2]: By "reduced" we mean that the proper calibrations are applied and a flux, flux density, or magnitude is extracted from the raw data.

# Software Impact and Conclusions

Moving forward, OTTER, in its entirety, will be a useful infrastructure tool for
time domain science. Even more, the OTTER API described here will
make access to that dataset easier for users by lowering the API learning curve.
Evidence of this is the impact of the Open Astronomy Catalogs, which has $>500$ citations
[@guillochon_open_2017] and is still used today, despite being deprecated.

There are already multiple astronomers using the software for their research, spanning
from undergraduate students to faculty. There are currently two papers citing
OTTER [@alexander_multi-wavelength_2025; @{2025arXiv250914317C}] and at least another in preparation
(Farley et al., in prep.). Additionally, our immediate research
groups has already used OTTER for writing successful telescope observing proposals. We presented
this work at the Kavli Institute for Theoretical Physics: Towards a Physical Understanding of
Tidal Disruption Events session, Astronomical Data Analysis Software and Systems: 2025 Monsoon
Workshop, and the 2025 X-ray Quasi-Periodic Eruptions and Repeating Nuclear
Transients Conference. It was positively received at all of these conferences and we
accrued $\sim 15$ beta testers who provided invaluable feedback. We welcome GitHub
issues with comments and feedback (or even pull requests!) from the community on
our [GitHub repository](https://github.com/astro-otter/otter).

# Acknowledgements
N.F. acknowledges support from the National Science Foundation Graduate Research Fellowship
Program under Grant No. DGE-2137419. KDA acknowledges support provided by the NSF through
award SOSPA9-007 from the NRAO and award AST-2307668. KDA gratefully acknowledges support
from the Alfred P. Sloan Foundation. This research was supported in part by grant NSF PHY-2309135
to the Kavli Institute for Theoretical Physics (KITP).

# References
