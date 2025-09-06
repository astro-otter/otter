---
title: 'OTTER API: A Python Package for API Access to the OTTER Catalog'
authors:
  - name: Noah Franz
    orcid: 0000-0003-4537-3575
    equal-contrib: true
    affiliation: 1
affiliations:
  - name: Department of Astronomy and Steward Observatory, University of Arizona, 933 North Cherry Avenue, Tucson, AZ 85721-0065, USA
    index: 1
date: 29 August 2025
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
events, and many more weird ones! However, understanding the many astrophysical transients observed
is not that easy and typically requires detailed multiwavelength observations and
analyses of the events.

# Citations
