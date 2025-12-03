OTTER Data Descriptions
=======================

Photometry Table Column Definitions
-----------------------------------
Generally, all columns that start with "converted" will be in the requested unit
and should be used. These are the columns that are returned with the `Otter.get_phot`
method. Here are those column descriptions:

==================== ==================
Column Name          Column Description
==================== ==================
name                 The default name of the transient.
converted_flux       The flux/flux density/magnitude in the requested units.
converted_flux_err   The error corresponding to the `converted_flux` column, in the requested units.
converted_date       The date the flux was measured, in the requested format.
converted_wave       The effective (central) wavelength of the filter.
converted_freq       The effective (central) frequency of the filter.
converted_flux_unit  The units of the `converted_flux` column.
converted_date_unit  The units of the `converted_date` column.
converted_wave_unit  The units of the `converted_wave` column.
converted_freq_unit  The units of the `converted_freq` column.
filter_name          The name of the filter (or band) used for the measurement of the flux value.
obs_type             Either "uvoir", "xray", or "radio"; describes the wavelength regime where the flux was measured.
upperlimit           True if the `converted_flux` column is the limiting flux rather than a flux measurement, False if the `converted_flux` is a detection.
reference            The ADS bibcodes(s) for this flux measurement. This can either be a string or list of strings.
human_readable_refs  The "human readable" version of the ADS bibcodes in the `reference` column, useful for labelling e.g., plots.
telescope            The telescope used to make this measurement if we have it stored, otherwise NaN.
==================== ==================

The other columns that are returned when calling either `Otter.get_phot` with the `raw=True`
option or from `Transient.clean_photometry` are described below:

==================== ==================
Column Name          Column Description
==================== ==================
freq_units           The frequency effective _raw_ units for this filter (which may differ from the `converted_freq_unit` column!).
val_k                The value of the k-correction applied to the data, if `corr_k = True`. Otherwise this is NaN.
corr_host            True if the photometry has been host subtracted, False if we know that it has not been, NaN if it is unclear.
wave_eff             The wavelength effective, in _raw_ `wave_units` units, for this filter (which may differ from the `converted_wave_unit` column!). This may change in length from the `converted_wave` column since we only store _either_ the `wave_eff` or `freq_eff` for a given filter.
_flux_err            The error on the flux units in some _unconverted_ format. The is prepended with an underscore because the user should almost always use the `converted_flux` column instead!
date                 The observation date in a _raw_ format. The user should default to the `converted_date` column.
wave_min             The minimum wavelength for this observation. Typically only used for X-ray observations.
raw_err_detail       A dictionary with other details on the raw uncertainty on the flux measurement. This dictionary stores things like asymmetric errors and a breakdown of the statistical and systematic errors that are in the `raw_err`.
wave_units           The units on the _raw_ (or "stored") wavelength for this filter. This can be NaN when wave_eff is NaN. The user should default to the `converted_wave` column.
corr_av              True if the photometry has been corrected for milky way extinction, False if we know that it has not been, NaN if it is unclear.
freq_eff             The frequency effective, in _raw_ `freq_units` units, for this filter (which may differ from the `converted_freq_unit` column!). This may change in length from the `converted_wave` column since we only store _either_ the `wave_eff` or `freq_eff` for a given filter.
xray_model           A dictionary with details about the xray model used to convert from the raw counts to flux value.
raw_err              The error on the `raw` column.
raw_units            The units on the `raw` column.
corr_s               True if the photometry has been s-corrected, False if we know that it has not been, NaN if it is unclear.
computed             True if we computed the `raw` or `value` columns from other information presented in the corresponding `reference`.
_flux_units          The units on the _unconverted_ flux. The is prepended with an underscore because the user should almost always use the `converted_flux` column instead!
val_hostav           The value of the host AV applied to the data, if `corr_hostav = True`. Otherwise this is NaN.
val_host             The value of the host subtraction applied to the data, if `corr_host = True`. Otherwise this is NaN.
val_s                The value of the s-correction applied to the data, if `corr_s = True`. Otherwise this is NaN.
value_err            The error on the `value` column, in the same units as the `value` column.
instrument           The instrument used to make this flux measurement, in some cases this is repetitive with the `telescope` column.
wave_max             The maximum wavelength for this observation. Typically only used for X-ray observations.
val_av               The value of the Milky Way AV applied to the data, if `corr_av = True`. Otherwise this is NaN.
value                The flux value after additional corrections. For example, this is used for X-ray photometry when we store both the raw count rate (in the `raw` column) and the flux value (in the `flux` column) after conversion.
corr_k               True if the photometry has been k-corrected, False if we know that it has not been, NaN if it is unclear.
corr_hostav          True if the photometry has been host extinction corrected, False if we know that it has not been, NaN if it is unclear.
raw                  The raw flux value, can be in units of counts or countrate, unlike the converted flux values or the `value` column. The user should default to using the `converted_flux` column.
filter_key           A unique identifier for this filter. This should be ignored and the `filter_name` column should be used instead.
value_units          The units on the `value` column.
date_format          The format of the _raw_ and _unconverted_ date. The user should default to the `converted_date` column.
_flux                The flux value in some _unconverted_ format. The is prepended with an underscore because the user should almost always use the `converted_flux` column instead!
==================== ==================

OTTER Data Schema
-----------------

.. autopydantic_model:: otter.schema.OtterSchema

.. autopydantic_model:: otter.schema.VersionSchema

.. autopydantic_model:: otter.schema.NameSchema

.. autopydantic_model:: otter.schema.CoordinateSchema

.. autopydantic_model:: otter.schema.ClassificationSchema

.. autopydantic_model:: otter.schema.ReferenceSchema

.. autopydantic_model:: otter.schema.DateSchema

.. autopydantic_model:: otter.schema.PhotometrySchema

.. autopydantic_model:: otter.schema.FilterSchema

.. autopydantic_model:: otter.schema.HostSchema

.. autopydantic_model:: otter.schema._AliasSchema

.. autopydantic_model:: otter.schema._XrayModelSchema

.. autopydantic_model:: otter.schema._ErrDetailSchema
