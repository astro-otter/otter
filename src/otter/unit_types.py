"""
Some convenient unit type classes for conversions
"""

import warnings
import numpy as np

from astropy.units import Quantity, Unit, UnitBase, StructuredUnit, LogUnit
from astropy.units.format import Base
import astropy.units as u
import astropy.constants as const
from collections.abc import Iterable

from .util import XRAY_AREAS


class _CustomQuantity(Quantity):
    """
    Customizing the astropy Quantity class to allow it to also handle Magnitude units

    This is potentially dangerous so we just have to be careful to not do any complex
    math after conversion
    """

    def _set_unit(self, unit):
        """Set the unit.

        This is used anywhere the unit is set or modified, i.e., in the
        initializer, in ``__imul__`` and ``__itruediv__`` for in-place
        multiplication and division by another unit, as well as in
        ``__array_finalize__`` for wrapping up views.  For Quantity, it just
        sets the unit, but subclasses can override it to check that, e.g.,
        a unit is consistent.

        **Modified from astropy source code**
        """
        if not isinstance(unit, UnitBase):
            if isinstance(self._unit, StructuredUnit) or isinstance(
                unit, StructuredUnit
            ):
                unit = StructuredUnit(unit, self.dtype)
            else:
                # Trying to go through a string ensures that, e.g., Magnitudes with
                # dimensionless physical unit become Quantity with units of mag.
                unit = Unit(str(unit), parse_strict="silent")
                if not isinstance(unit, (UnitBase, StructuredUnit, LogUnit)):
                    raise UnitTypeError(
                        f"{self.__class__.__name__} instances require normal units, "
                        f"not {unit.__class__} instances."
                    )

        self._unit = unit


class Flux(_CustomQuantity):
    def __init__(self, quantity):
        """
        This is a Flux Unit Type and should have units of Energy / Distance^2 / Time
        """
        super(Flux, self).__init__()

        self._name = "energy flux"

    @staticmethod
    def isflux(quantity):
        """
        Tests if this is a flux
        """
        _name = "energy flux"

        if not isinstance(quantity, list):
            quantity = [quantity]

        return all(_name in q.unit.physical_type for q in quantity)

    def tofluxdensity(
        self,
        wave_eff: Quantity = None,
        freq_eff: Quantity = None,
        out_units: Unit = "mag(AB)",
        **kwargs,
    ):
        """
        Converts to flux density by dividing by the given effective wave/freq

        Args:
            freq_eff [astropy.units.Quantity]: An astropy Quantity with the effective frequency
            wave_eff [astropy.units.Quantity]: An astropy Quantity with the effective wavelength

        Returns:
            A FluxDensity Quantity Object with the adjusted data
        """

        # check that the input unit is a flux density
        if not FluxDensity.isfluxdensity(1 * Unit(out_units)):
            raise ValueError(
                "The output units that were given are not a flux density!!"
            )

        # perform the conversion to out_units
        if freq_eff is not None and wave_eff is None:
            if not isinstance(freq_eff, np.ndarray):
                freq_eff = np.array(freq_eff)

            out = (self / freq_eff).to(
                out_units, equivalencies=u.spectral_density(1 * freq_eff.unit)
            )

        elif freq_eff is None and wave_eff is not None:
            if not isinstance(wave_eff, np.ndarray):
                wave_eff = np.array(wave_eff)

            out = (self / wave_eff).to(
                out_units, equivalencies=u.spectral_density(1 * wave_eff.unit)
            )

        else:
            raise ValueError(
                "Either freq_eff or wave_eff must be provided (NOT both!)"
            )

        # double check this is a fluxdensity now
        if not FluxDensity.isfluxdensity(out):
            raise ValueError(
                "Something went wrong! This is not a flux density still!"
            )

        return out

    def toflux(self, out_units=u.erg / u.cm**2 / u.s, **kwargs):
        """
        Just a wrapper that returns itself. This will make my other code cleaner
        """
        return self.to(out_units)

    def tocountrate(
        self,
        freq_eff: Quantity = None,
        wave_eff: Quantity = None,
        out_units=1 / u.s,
        **kwargs,
    ):
        """
        Convert Flux to a count rate
        """
        raise ValueError(
            "Converting a flux to a count rate is currently not "
            + "supported!"
        )


class FluxDensity(_CustomQuantity):
    def __init__(self, quantity):
        """
        This is a Flux Unit Type and should have units of a Energy / Distance^2 / Time / Frequency

        This also handles AB/ST Magnitudes depending on if it is in frequency or
        wavelength space
        """
        super(FluxDensity, self).__init__()

        self._name_wave = "spectral flux density wav"
        self._name_freq = "spectral flux density"

    @staticmethod
    def isfluxdensity(quantity):
        """
        Tests if this is a flux density
        """
        _name_wave = "spectral flux density wav"
        _name_freq = "spectral flux density"

        if not isinstance(quantity, list):
            quantity = [quantity]

        test1 = all(_name_freq in q.unit.physical_type for q in quantity)
        test2 = all(_name_wave in q.unit.physical_type for q in quantity)
        return test1 or test2

    def toflux(
        self,
        freq_eff: Quantity = None,
        wave_eff: Quantity = None,
        out_units: Unit = u.erg / u.cm**2 / u.s,
        **kwargs,
    ) -> Flux:
        """
        Converts to flux by multiplying by the fiven effective wave/freq

        Args:
            freq_eff [astropy.units.Quantity]: An astropy Quantity with the effective frequency
            wave_eff [astropy.units.Quantity]: An astropy Quantity with the effective wavelength
            out_units [astropy.units.Unit]: An astropy Unit with units of flux

        Returns:
            A Flux Quantity Object with the adjusted data
        """

        # check that the input unit is a flux
        if not Flux.isflux(1 * Unit(out_units)):
            raise ValueError(
                "The output units that were given are not a flux!!"
            )

        # perform the conversion
        if freq_eff is not None and wave_eff is None:
            if not isinstance(freq_eff, np.ndarray):
                freq_eff = np.array(freq_eff)

            if "wav" in self.unit.physical_type:
                eff = const.c / freq_eff  # convert to wavelengths
            else:
                eff = freq_eff
            out = (self * eff).to(
                out_units, equivalencies=u.spectral_density(1 * freq_eff.unit)
            )

        elif freq_eff is None and wave_eff is not None:
            if not isinstance(wave_eff, np.ndarray):
                wave_eff = np.array(wave_eff)

            if "wav" not in self.unit.physical_type:
                eff = const.c / wave_eff  # convert to frequency
            else:
                eff = wave_eff
            out = (self * eff).to(
                out_units, equivalencies=u.spectral_density(1 * wave_eff.unit)
            )

        else:
            raise ValueError(
                "Either freq_eff or wave_eff must be provided (NOT both!)"
            )

        # double check this is a flux now
        if not Flux.isflux(out):
            raise ValueError("Something went wrong! This is not a flux still!")

        return Flux(out)

    def tofluxdensity(
        self,
        out_units="mag(AB)",
        freq_eff: Quantity = None,
        wave_eff: Quantity = None,
        **kwargs,
    ):
        """
        Just a wrapper that returns itself. This will make my other code cleaner
        """

        wav = self._name_wave
        nowav = self._name_freq

        # check if out_units is in frequency or wavelength space
        if (
            wav in Unit(out_units).physical_type
            and nowav in self.unit.physical_type
        ):
            # There is a discrepancy
            # we have to put self in wavelength space
            if wave_eff is not None:
                # this is ideal, nothing needs to be done
                pass
            elif freq_eff is not None:
                # we can work with this
                wave_eff = const.c / freq_eff
            else:
                raise ValueError(
                    "We need a wave_eff or freq_eff to make this conversion!"
                )

            eq = u.spectral_density(wave_eff)

        elif (
            nowav in Unit(out_units).physical_type
            and wav in self.unit.physical_type
        ):
            if freq_eff is not None:
                pass
            elif wave_eff is not None:
                freq_eff = const.c / wave_eff
            else:
                raise ValueError("We need a wave_eff to make this conversion!")
            eq = u.spectral_density(freq_eff)

        else:
            eq = None

        return self.to(out_units, equivalencies=eq)

    def tocountrate(
        self,
        out_units=1 / u.s,
        freq_eff: Quantity = None,
        wave_eff: Quantity = None,
        **kwargs,
    ):
        """
        Convert the flux density to a count rate
        """
        raise ValueError(
            "Converting flux density to a count rate is currently not "
            + "supported!"
        )


class CountRate(_CustomQuantity):
    """
    Represents the counts of electrons in an x-ray detector
    """

    def __init__(self, quantity):
        """
        This is a Count Rate Unit Type and should have units of a 1 / Time
        """
        super(CountRate, self).__init__()

        self._name = "frequency"  # the name astropy gave to this

    @staticmethod
    def iscountrate(quantity):
        """
        Checks if the input is actually a count rate
        """
        _name = "frequency"

        if not isinstance(quantity, list):
            quantity = [quantity]

        return all(_name in q.unit.physical_type for q in quantity)

    def toflux(
        self,
        wave_eff: Quantity = None,
        freq_eff: Quantity = None,
        out_units=u.erg / u.s / u.cm**2,
        telescope: str = None,
        **kwargs,
    ):
        """
        Converts the count rate to a flux
        """
        warnings.warn(
            "Converting count rate to a flux is still under development! "
            + "This is just an order of magnitude estimate!"
        )

        if telescope is None:
            raise ValueError(
                "Can not convert count rate to flux without a telescope!"
            )

        # check that the input unit is a flux
        if not Flux.isflux(1 * Unit(out_units)):
            raise ValueError(
                "The output units that were given are not a flux!!"
            )

        # perform the conversion
        if freq_eff is not None and wave_eff is None:
            if not isinstance(freq_eff, np.ndarray):
                freq_eff = np.array(freq_eff)

            out_erg_s = self * (freq_eff.to(u.erg, equivalencies=u.spectral()))
            out_erg_s_cm2 = out_erg_s / XRAY_AREAS[telescope.lower()]
            out = out_erg_s_cm2.to(
                out_units, equivalencies=u.spectral_density(1 * freq_eff.unit)
            )

        elif freq_eff is None and wave_eff is not None:
            if not isinstance(wave_eff, np.ndarray):
                wave_eff = np.array(wave_eff)

            out_erg_s = self * (wave_eff.to(u.erg, equivalencies=u.spectral()))
            out_erg_s_cm2 = out_erg_s / XRAY_AREAS[telescope.lower()]
            out = out_erg_s_cm2.to(
                out_units, equivalencies=u.spectral_density(1 * wave_eff.unit)
            )

        else:
            raise ValueError(
                "Either freq_eff or wave_eff must be provided (NOT both!)"
            )

        # double check this is a flux now
        if not Flux.isflux(out):
            raise ValueError("Something went wrong! This is not a flux still!")

        return Flux(out)

    def tofluxdensity(
        self,
        wave_eff: Quantity = None,
        freq_eff: Quantity = None,
        out_units=u.erg / u.s / u.cm**2 / u.Hz,
        telescope: str = None,
        **kwargs,
    ):
        """
        Converts the count rate to a flux density
        """

        warnings.warn(
            "Converting count rate to a flux density is still under development! "
            + "This is just an order of magnitude estimate!"
        )

        if telescope is None:
            raise ValueError(
                "Can not convert count rate to flux density without a telescope!"
            )

        # check that the input unit is a flux
        if not FluxDensity.isfluxdensity(1 * Unit(out_units)):
            raise ValueError(
                "The output units that were given are not a flux!!"
            )

        # perform the conversion
        if freq_eff is not None and wave_eff is None:
            if not isinstance(freq_eff, np.ndarray):
                freq_eff = np.array(freq_eff)

            out_erg_s = self * (freq_eff.to(u.erg, equivalencies=u.spectral()))
            out_erg_s_cm2 = out_erg_s / XRAY_AREAS[telescope.lower()]
            out_erg_s_cm2_freq = out_erg_s_cm2 / freq_eff
            out = out_erg_s_cm2_freq.to(
                out_units, equivalencies=u.spectral_density(1 * freq_eff.unit)
            )

        elif freq_eff is None and wave_eff is not None:
            if not isinstance(wave_eff, np.ndarray):
                wave_eff = np.array(wave_eff)

            out_erg_s = self * (wave_eff.to(u.erg, equivalencies=u.spectral()))
            out_erg_s_cm2 = out_erg_s / XRAY_AREAS[telescope.lower()]
            out_erg_s_cm2_wav = out_erg_s_cm2 / wave_eff
            out = out_erg_s_cm2_wav.to(
                out_units, equivalencies=u.spectral_density(1 * wave_eff.unit)
            )

        else:
            raise ValueError(
                "Either freq_eff or wave_eff must be provided (NOT both!)"
            )

        # double check this is a flux now
        if not FluxDensity.isfluxdensity(out):
            raise ValueError(
                "Something went wrong! This is not a flux density still!"
            )

        try:
            return FluxDensity(out)
        except:
            import pdb

            pdb.set_trace()

    def tocountrate(self, out_units=1 / u.s, **kwargs):
        """
        Just a wrapper to return the count rate in out_units
        """
        return self.to(out_units)


def get_type(quantity):
    """
    Finds the type of the input quantity and returns the quantity as that class
    """

    if FluxDensity.isfluxdensity(quantity):
        # we have to be a little careful cause some of these are mags
        if str(quantity.unit) == "mag(AB)":
            print("fouund magab")
            return FluxDensity(quantity.to(u.erg / u.cm**2 / u.s / u.Hz))
        elif str(quantity.unit) == "mag(ST)":
            return FluxDensity(quantity.to(u.erg / u.cm**2 / u.s / u.nm))
        else:
            return FluxDensity(quantity)
    elif Flux.isflux(quantity):
        return Flux(quantity)
    elif CountRate.iscountrate(quantity):
        return CountRate(quantity)
    else:
        raise ValueError(f"{quantity} does not have valid units!")
