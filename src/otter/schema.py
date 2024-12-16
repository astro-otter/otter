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


class NameSchema(BaseModel):
    default_name: str
    alias: list[_AliasSchema]


class CoordinateSchema(BaseModel):
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
    reference: Union[List[str], str]
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


class ClassificationSchema(BaseModel):
    object_class: str
    confidence: float
    reference: Union[str, List[str]]
    default: bool = False
    class_type: str = None


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
    date: Union[str, List[str], float, List[float]]
    date_format: Union[str, List[str]]
    date_err: Optional[Union[str, List[str], float, List[float]]] = None
    ignore: Optional[Union[bool, List[bool]]] = None
    upperlimit: Optional[Union[bool, List[bool]]] = None
    sigma: Optional[Union[str, List[str], float, List[float]]] = None
    sky: Optional[Union[str, List[str], float, List[float]]] = None
    telescope: Optional[Union[str, List[str]]] = None
    instrument: Optional[Union[str, List[str]]] = None
    phot_type: Optional[Union[str, List[str]]] = None
    exptime: Optional[Union[str, List[str], int, List[int], float, List[float]]] = None
    aperture: Optional[Union[str, List[str], int, List[int], float, List[float]]] = None
    observer: Optional[Union[str, List[str]]] = None
    reducer: Optional[Union[str, List[str]]] = None
    pipeline: Optional[Union[str, List[str]]] = None
    corr_k: Optional[Union[bool, List[bool]]] = None
    corr_av: Optional[Union[bool, List[bool]]] = None
    corr_host: Optional[Union[bool, List[bool]]] = None
    corr_hostav: Optional[Union[bool, List[bool]]] = None
    val_k: Optional[Union[float, List[float], int, List[int]]] = None
    val_s: Optional[Union[float, List[float], int, List[int]]] = None
    val_av: Optional[Union[float, List[float], int, List[int]]] = None
    val_host: Optional[Union[float, List[float], int, List[int]]] = None
    val_hostav: Optional[Union[float, List[float], int, List[int]]] = None

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


class FilterSchema(BaseModel):
    filter_key: str
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


class OtterSchema(BaseModel):
    schema_version: VersionSchema
    name: NameSchema
    coordinate: list[CoordinateSchema]
    distance: list[DistanceSchema]
    classification: list[ClassificationSchema]
    reference_alias: list[ReferenceSchema]
    date_reference: list[DateSchema]
    photometry: list[PhotometrySchema]
    filter_alias: list[FilterSchema]
