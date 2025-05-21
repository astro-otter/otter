"""
Pydantic Schema Model of our JSON schema
"""

from pydantic import BaseModel, model_validator, field_validator, ValidationError
from typing import Optional, Union, List


class VersionSchema(BaseModel):
    value: Union[str, int] = None
    comment: str = None


class _AliasSchema(BaseModel):
    value: str
    reference: Union[str, List[str]]


class _XrayModelSchema(BaseModel):
    # the following two lines are needed to prevent annoying warnings
    model_config: dict = {}
    model_config["protected_namespaces"] = ()

    # required keywords
    model_name: str
    param_names: List[str]
    param_values: List[Union[float, int, str]]
    param_units: List[Union[str, None]]
    min_energy: Union[float, int, str]
    max_energy: Union[float, int, str]
    energy_units: str

    # optional keywords
    param_value_err_upper: Optional[List[Union[float, int, str]]] = None
    param_value_err_lower: Optional[List[Union[float, int, str]]] = None
    param_upperlimit: Optional[List[Union[float, int, str]]] = None
    param_descriptions: Optional[List[str]] = None
    model_reference: Optional[Union[str, List[str]]] = None


class _ErrDetailSchema(BaseModel):
    # all optional keywords!
    upper: Optional[List[Union[float, int, str]]] = None
    lower: Optional[List[Union[float, int, str]]] = None
    systematic: Optional[List[Union[float, int, str]]] = None
    statistical: Optional[List[Union[float, int, str]]] = None
    iss: Optional[List[Union[float, int, str]]] = None


class NameSchema(BaseModel):
    default_name: str
    alias: list[_AliasSchema]


class CoordinateSchema(BaseModel):
    reference: Union[List[str], str]
    ra: Union[str, float] = None
    dec: Union[str, float] = None
    l: Union[str, float] = None  # noqa: E741
    b: Union[str, float] = None
    lon: Union[str, float] = None
    lat: Union[str, float] = None
    ra_units: str = None
    dec_units: str = None
    l_units: str = None
    b_units: str = None
    lon_units: str = None
    lat_units: str = None
    ra_error: Union[str, float] = None
    dec_error: Union[str, float] = None
    l_error: Union[str, float] = None
    b_error: Union[str, float] = None
    lon_error: Union[str, float] = None
    lat_error: Union[str, float] = None
    epoch: str = None
    frame: str = "J2000"
    coord_type: str = None
    computed: bool = False
    default: bool = False

    @model_validator(mode="after")
    def _has_coordinate(self):
        uses_ra_dec = self.ra is not None and self.dec is not None
        uses_galactic = self.l is not None and self.b is not None
        uses_lon_lat = self.lon is not None and self.lat is not None

        if uses_ra_dec:
            if self.ra_units is None:
                raise ValidationError("ra_units must be provided for RA!")
            if self.dec_units is None:
                raise ValidationError("dec_units must be provided for Dec!")

        elif uses_galactic:
            if self.l_units is None:
                raise ValidationError("l_units must be provided for RA!")
            if self.b_units is None:
                raise ValidationError("b_units must be provided for Dec!")

        elif uses_lon_lat:
            if self.lon_units is None:
                raise ValidationError("lon_units must be provided for RA!")
            if self.lat_units is None:
                raise ValidationError("lat_units must be provided for Dec!")

        else:
            ValidationError("Must have RA/Dec, l/b, and/or lon/lat!")

        return self


class DistanceSchema(BaseModel):
    value: Union[str, float, int]
    unit: str = None
    reference: Union[str, List[str]]
    distance_type: str
    error: Union[str, float, int] = None
    cosmology: str = None
    computed: bool = False
    uuid: str = None
    default: bool = False

    @model_validator(mode="after")
    def _has_units(self):
        if self.distance_type != "redshift" and self.unit is None:
            raise ValidationError("Need units if the distance_type is not redshift!")

        return self


class ClassificationSchema(BaseModel):
    object_class: str
    confidence: float
    reference: Union[str, List[str]]
    default: bool = False
    class_type: str = None


class ClassificationDictSchema(BaseModel):
    spec_classed: Optional[int] = None
    unambiguous: Optional[bool] = None
    value: list[ClassificationSchema]


class ReferenceSchema(BaseModel):
    name: str
    human_readable_name: str


class DateSchema(BaseModel):
    value: Union[str, int, float]
    date_format: str
    date_type: str
    reference: Union[str, List[str]]
    computed: bool = None


class PhotometrySchema(BaseModel):
    reference: Union[List[str], str]
    raw: list[Union[float, int]]
    raw_err: Optional[List[float]] = []
    raw_units: Union[str, List[str]]
    value: Optional[list[Union[float, int]]] = None
    value_err: Optional[list[Union[float, int]]] = None
    value_units: Optional[Union[str, List[str]]] = None
    epoch_zeropoint: Optional[Union[float, str, int]] = None
    epoch_redshift: Optional[Union[float, int]] = None
    filter: Optional[Union[str, List[str]]] = None
    filter_key: Union[str, List[str]]
    obs_type: Union[str, List[str]]
    telescope_area: Optional[Union[float, List[float]]] = None
    date: Union[str, float, List[Union[str, float]]]
    date_format: Union[str, List[str]]
    date_err: Optional[Union[str, float, List[Union[str, float]]]] = None
    ignore: Optional[Union[bool, List[bool]]] = None
    upperlimit: Optional[Union[bool, List[bool]]] = None
    sigma: Optional[Union[str, float, List[Union[str, float]]]] = None
    sky: Optional[Union[str, float, List[Union[str, float]]]] = None
    telescope: Optional[Union[str, List[str]]] = None
    instrument: Optional[Union[str, List[str]]] = None
    phot_type: Optional[Union[str, List[str]]] = None
    exptime: Optional[Union[str, int, float, List[Union[str, int, float]]]] = None
    aperture: Optional[Union[str, int, float, List[Union[str, int, float]]]] = None
    observer: Optional[Union[str, List[str]]] = None
    reducer: Optional[Union[str, List[str]]] = None
    pipeline: Optional[Union[str, List[str]]] = None
    corr_k: Optional[Union[bool, str, List[Union[bool, str]]]] = None
    corr_s: Optional[Union[bool, str, List[Union[bool, str]]]] = None
    corr_av: Optional[Union[bool, str, List[Union[bool, str]]]] = None
    corr_host: Optional[Union[bool, str, List[Union[bool, str]]]] = None
    corr_hostav: Optional[Union[bool, str, List[Union[bool, str]]]] = None
    val_k: Optional[Union[float, int, str, List[Union[float, int, str]]]] = None
    val_s: Optional[Union[float, int, str, List[Union[float, int, str]]]] = None
    val_av: Optional[Union[float, int, str, List[Union[float, int, str]]]] = None
    val_host: Optional[Union[float, int, str, List[Union[float, int, str]]]] = None
    val_hostav: Optional[Union[float, int, str, List[Union[float, int, str]]]] = None
    xray_model: Optional[Union[List[_XrayModelSchema], List[None]]] = None
    raw_err_detail: Optional[_ErrDetailSchema] = None
    value_err_detail: Optional[_ErrDetailSchema] = None

    @field_validator(
        "raw_units",
        "raw_err",
        "filter_key",
        "obs_type",
        "date_format",
        "upperlimit",
        "date",
        "telescope",
    )
    @classmethod
    def ensure_list(cls, v):
        if not isinstance(v, list):
            return [v]
        return v

    @model_validator(mode="after")
    def _ensure_xray_model(self):
        """
        This will eventually ensure the xray_model key is used if obs_type="xray"

        It will be commented out until we get the data setup correctly
        """
        # if self.obs_type == "xray" and self.xray_model is None:
        #     raise ValidationError(
        #         "Need an xray_model for this xray data!"
        #     )

        return self


class FilterSchema(BaseModel):
    filter_key: str
    filter_name: str
    wave_eff: Union[str, float, int] = None
    wave_min: Union[str, float, int] = None
    wave_max: Union[str, float, int] = None
    freq_eff: Union[str, float, int] = None
    freq_min: Union[str, float, int] = None
    freq_max: Union[str, float, int] = None
    zp: Union[str, float, int] = None
    wave_units: Union[str, float, int] = None
    freq_units: Union[str, float, int] = None
    zp_units: Union[str, float, int] = None
    zp_system: Union[str, float, int] = None


class HostSchema(BaseModel):
    reference: Union[str, List[str]]
    host_ra: Optional[Union[str, float]] = None
    host_dec: Optional[Union[str, float]] = None
    host_ra_units: Optional[str] = None
    host_dec_units: Optional[str] = None
    host_z: Optional[Union[str, int, float]] = None
    host_type: Optional[str] = None
    host_name: Optional[str] = None

    @model_validator(mode="after")
    def _has_coordinate_or_name(self):
        has_coordinate = self.host_ra is not None and self.host_dec is not None
        has_name = self.host_name is not None

        # if it has the RA/Dec keys, make sure it also has ra_unit, dec_unit keys
        if has_coordinate:
            if self.host_ra_units is None:
                raise ValidationError("Need RA unit if coordinates are provided!")
            if self.host_dec_units is None:
                raise ValidationError("Need Dec unit if coordinates are provided!")

        # we need either the coordinate or name to identify this object
        # Both are okay too (more info is always better)
        if not has_coordinate and not has_name:
            raise ValidationError(
                "Need to provide a Host name and/or host coordinates!"
            )

        # Make sure that if one of RA/Dec is given then both are given
        if (self.host_ra is None and self.host_dec is not None) or (
            self.host_ra is not None and self.host_dec is None
        ):
            raise ValidationError(
                "Please provide RA AND Dec, not just one or the other!"
            )

        return self


class OtterSchema(BaseModel):
    schema_version: Optional[VersionSchema] = None
    name: NameSchema
    coordinate: list[CoordinateSchema]
    distance: Optional[list[DistanceSchema]] = None
    classification: Optional[ClassificationDictSchema] = None
    reference_alias: list[ReferenceSchema]
    date_reference: Optional[list[DateSchema]] = None
    photometry: Optional[list[PhotometrySchema]] = None
    filter_alias: Optional[list[FilterSchema]] = None
    host: Optional[list[HostSchema]] = None

    @model_validator(mode="after")
    def _verify_filter_alias(self):
        if self.photometry is not None and self.filter_alias is None:
            raise ValidationError("filter_alias is needed if photometry is given!")
