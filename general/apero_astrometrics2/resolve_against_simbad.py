#!/usr/bin/env python3
"""
Resolve astrometrics YAML files against SIMBAD / Gaia / VizieR TAP services.

For each YAML in the astrometrics/ directory:
  - Resolve the object name in SIMBAD and fill SIMBAD_NAME
  - Fill any null photometry (G, GBP, GRP, J, H, Ks, W1, W2, W3, W4)
  - Fill parallax, proper motion, RV, Teff, SpT if missing
  - Write back the YAML, touching only null values (existing data is preserved)

Keep the YAML I/O and the SIMBAD resolution code clearly separated so that
resolve_from_name() can be called directly in future for new targets.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# TAP endpoints (must stay near the top of the file)
# ---------------------------------------------------------------------------
SIMBAD_TAP = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"
VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
# ---------------------------------------------------------------------------

import argparse
import csv
import io
import json
import math
import os
import sys
import time
from pathlib import Path
from typing import Any
from urllib.parse import quote_plus
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import yaml


OBLIQUITY_DEG = 23.4392911
EARTH_ORBITAL_SPEED_KMS = 29.79
TELLURIC_THRESHOLD_KMS = 5.0
NON_LEAP_MONTH_LENGTHS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
NON_LEAP_MONTH_NAMES = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
EQUATORIAL_TO_GALACTIC_MATRIX = [
    [-0.0548755604, -0.8734370902, -0.4838350155],
    [0.4941094279, -0.44482963, 0.7469822445],
    [-0.867666149, -0.1980763734, 0.4559837762],
]


# ---------------------------------------------------------------------------
# Helpers – value normalisation
# ---------------------------------------------------------------------------

def _nv(value: Any) -> Any:
    """Return None for empty / null-like values, otherwise the stripped string."""
    if value is None:
        return None
    s = str(value).strip()
    if s == "" or s.lower() in {"none", "null", "nan"}:
        return None
    return s


def _pf(value: Any) -> float | None:
    s = _nv(value)
    if s is None:
        return None
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------------------
# Low-level TAP helpers (synchronous, plain urllib)
# ---------------------------------------------------------------------------

def _tap_get(
    endpoint: str,
    adql: str,
    fmt: str = "json",
    timeout: int = 30,
    request_key: str = "request",
    lang_key: str = "lang",
    format_key: str = "format",
    query_key: str = "query",
) -> bytes | None:
    params = (
        f"{request_key}=doQuery"
        f"&{lang_key}=adql"
        f"&{format_key}={fmt}"
        f"&{query_key}={quote_plus(adql)}"
    )
    url = f"{endpoint}?{params}"
    req = Request(url, headers={"User-Agent": "apero-astrometrics/1.0"})
    try:
        with urlopen(req, timeout=timeout) as resp:
            return resp.read()
    except HTTPError as exc:
        detail = ""
        try:
            body = exc.read().decode("utf-8", errors="replace")
            if body:
                detail = f" | {body[:300].replace(chr(10), ' ').replace(chr(13), ' ')}"
        except Exception:
            detail = ""
        print(f"  [TAP WARNING] {endpoint} → {exc}{detail}", file=sys.stderr)
        return None
    except (URLError, HTTPError, OSError) as exc:
        print(f"  [TAP WARNING] {endpoint} → {exc}", file=sys.stderr)
        return None


def _simbad_json(adql: str, timeout: int = 30) -> dict | None:
    raw = _tap_get(SIMBAD_TAP, adql, fmt="json", timeout=timeout)
    if raw is None:
        return None
    try:
        return json.loads(raw)
    except Exception:
        return None


def _vizier_json(adql: str, timeout: int = 30) -> dict | None:
    # VizieR TAP can be strict on parameter names; use uppercase keys.
    raw = _tap_get(
        VIZIER_TAP,
        adql,
        fmt="json",
        timeout=timeout,
        request_key="REQUEST",
        lang_key="LANG",
        format_key="FORMAT",
        query_key="QUERY",
    )
    if raw is None:
        return None
    try:
        return json.loads(raw)
    except Exception:
        return None


def _gaia_csv(adql: str, timeout: int = 30) -> list[dict] | None:
    """Query Gaia TAP in CSV format; return a list of dicts."""
    raw = _tap_get(GAIA_TAP, adql, fmt="csv", timeout=timeout)
    if raw is None:
        return None
    try:
        text = raw.decode("utf-8", errors="replace")
        reader = csv.DictReader(io.StringIO(text))
        return list(reader)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def _sep_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    cos_sep = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(r1 - r2)
    cos_sep = max(-1.0, min(1.0, cos_sep))
    return math.degrees(math.acos(cos_sep)) * 3600.0


def _propagate(ra_deg: float, dec_deg: float,
               pmra: float | None, pmdec: float | None,
               dt_yr: float) -> tuple[float, float]:
    if pmra is None or pmdec is None:
        return ra_deg, dec_deg
    cos_dec = math.cos(math.radians(dec_deg))
    if abs(cos_dec) < 1e-8:
        return ra_deg, dec_deg
    new_ra = (ra_deg + (pmra * dt_yr) / (3.6e6 * cos_dec)) % 360.0
    new_dec = max(-90.0, min(90.0, dec_deg + (pmdec * dt_yr) / 3.6e6))
    return new_ra, new_dec


def _jd_to_jyear(epoch_jd: float | None) -> float | None:
    if epoch_jd is None:
        return None
    try:
        return float(Time(epoch_jd, format="jd").jyear)
    except Exception:
        return None


def _propagate_to_epoch(
    ra_deg: float | None,
    dec_deg: float | None,
    pmra: float | None,
    pmdec: float | None,
    source_epoch_jyear: float | None,
    target_epoch_jyear: float,
    plx_mas: float | None = None,
    rv_kms: float | None = None,
) -> tuple[float | None, float | None]:
    """Propagate ICRS coordinates between epochs using astropy, with safe fallback."""
    if ra_deg is None or dec_deg is None:
        return ra_deg, dec_deg
    if pmra is None or pmdec is None or source_epoch_jyear is None:
        return ra_deg, dec_deg

    try:
        kwargs: dict[str, Any] = {
            "pm_ra_cosdec": pmra * u.mas / u.yr,
            "pm_dec": pmdec * u.mas / u.yr,
            "obstime": Time(source_epoch_jyear, format="jyear"),
        }
        if plx_mas is not None and plx_mas > 0:
            kwargs["distance"] = (1000.0 / plx_mas) * u.pc
        if rv_kms is not None:
            kwargs["radial_velocity"] = rv_kms * u.km / u.s

        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs", **kwargs)
        moved = coord.apply_space_motion(new_obstime=Time(target_epoch_jyear, format="jyear"))
        return float(moved.ra.deg), float(moved.dec.deg)
    except Exception:
        # Fallback to simple PM-only propagation if astropy propagation cannot be applied.
        return _propagate(ra_deg, dec_deg, pmra, pmdec, target_epoch_jyear - source_epoch_jyear)


def _set_if_none(d: dict, key: str, value: Any) -> None:
    """Set d[key] = value only if d[key] is currently None."""
    if d.get(key) is None and value is not None:
        d[key] = value


def _format_ra_hms(ra_deg: float | None) -> str | None:
    if ra_deg is None:
        return None
    total_seconds = (ra_deg % 360.0) * 240.0
    hours = int(total_seconds // 3600)
    minutes = int((total_seconds % 3600) // 60)
    seconds = total_seconds - 3600 * hours - 60 * minutes
    if seconds >= 59.995:
        seconds = 0.0
        minutes += 1
    if minutes >= 60:
        minutes = 0
        hours = (hours + 1) % 24
    return f"{hours:02d}:{minutes:02d}:{seconds:06.3f}"


def _format_dec_dms(dec_deg: float | None) -> str | None:
    if dec_deg is None:
        return None
    sign = "+" if dec_deg >= 0 else "-"
    total_arcsec = abs(dec_deg) * 3600.0
    degrees = int(total_arcsec // 3600)
    minutes = int((total_arcsec % 3600) // 60)
    seconds = total_arcsec - 3600 * degrees - 60 * minutes
    if seconds >= 59.995:
        seconds = 0.0
        minutes += 1
    if minutes >= 60:
        minutes = 0
        degrees += 1
    return f"{sign}{degrees:02d}:{minutes:02d}:{seconds:05.2f}"


def _galactic_from_radec(ra_deg: float | None, dec_deg: float | None) -> tuple[float | None, float | None]:
    if ra_deg is None or dec_deg is None:
        return None, None
    ra_r = math.radians(ra_deg)
    dec_r = math.radians(dec_deg)
    eq_vec = [
        math.cos(dec_r) * math.cos(ra_r),
        math.cos(dec_r) * math.sin(ra_r),
        math.sin(dec_r),
    ]
    gal_vec = [
        sum(EQUATORIAL_TO_GALACTIC_MATRIX[i][j] * eq_vec[j] for j in range(3))
        for i in range(3)
    ]
    gal_lon = math.degrees(math.atan2(gal_vec[1], gal_vec[0])) % 360.0
    gal_lat = math.degrees(math.asin(max(-1.0, min(1.0, gal_vec[2]))))
    return gal_lon, gal_lat


def _ecliptic_from_radec(ra_deg: float | None, dec_deg: float | None) -> tuple[float | None, float | None]:
    if ra_deg is None or dec_deg is None:
        return None, None
    ra_r = math.radians(ra_deg)
    dec_r = math.radians(dec_deg)
    eps_r = math.radians(OBLIQUITY_DEG)
    sin_beta = math.sin(dec_r) * math.cos(eps_r) - math.cos(dec_r) * math.sin(eps_r) * math.sin(ra_r)
    beta_r = math.asin(max(-1.0, min(1.0, sin_beta)))
    y = math.sin(ra_r) * math.cos(eps_r) + math.tan(dec_r) * math.sin(eps_r)
    x = math.cos(ra_r)
    lambda_r = math.atan2(y, x)
    return math.degrees(lambda_r) % 360.0, math.degrees(beta_r)


def _doy_label(day_index: int) -> str:
    day_number = day_index + 1
    remaining = day_number
    for month_name, month_length in zip(NON_LEAP_MONTH_NAMES, NON_LEAP_MONTH_LENGTHS):
        if remaining <= month_length:
            return f"{month_name} {remaining:02d}"
        remaining -= month_length
    return "Dec 31"


def _telluric_windows(ra_deg: float | None, dec_deg: float | None, rv_kms: float | None) -> str | None:
    ecl_lon_deg, ecl_lat_deg = _ecliptic_from_radec(ra_deg, dec_deg)
    if ecl_lon_deg is None or ecl_lat_deg is None or rv_kms is None:
        return None

    flagged = []
    for day_index in range(365):
        sun_lon_deg = ((day_index + 1) - 80.0) * 360.0 / 365.0
        vbary = EARTH_ORBITAL_SPEED_KMS * math.cos(math.radians(ecl_lat_deg)) * math.sin(
            math.radians(sun_lon_deg - ecl_lon_deg)
        )
        flagged.append(abs(rv_kms + vbary) < TELLURIC_THRESHOLD_KMS)

    ranges: list[tuple[int, int]] = []
    start = None
    for day_index, is_flagged in enumerate(flagged):
        if is_flagged and start is None:
            start = day_index
        elif not is_flagged and start is not None:
            ranges.append((start, day_index - 1))
            start = None
    if start is not None:
        ranges.append((start, 364))

    if len(ranges) > 1 and flagged[0] and flagged[-1]:
        first_start, first_end = ranges[0]
        last_start, _last_end = ranges[-1]
        ranges = [(last_start, first_end)] + ranges[1:-1]

    if not ranges:
        return "always > 5 km/s"

    parts = []
    for r0, r1 in ranges:
        if r0 == r1:
            parts.append(_doy_label(r0))
        else:
            parts.append(f"{_doy_label(r0)} to {_doy_label(r1)}")
    return "; ".join(parts)


def _absolute_mag(apparent_mag: float | None, parallax_mas: float | None) -> float | None:
    if apparent_mag is None or parallax_mas is None or parallax_mas <= 0:
        return None
    distance_pc = 1000.0 / parallax_mas
    return apparent_mag - 5.0 * math.log10(distance_pc) + 5.0


def _teff_from_gaia_colors(gbp: float | None, grp: float | None, jmag: float | None, hmag: float | None) -> tuple[float | None, float | None]:
    """Return (TEFF_GAIA_JH, TEFF_GAIA) using Mann+2015 relations used in starometer."""
    if gbp is None or grp is None:
        return None, None
    col = gbp - grp
    if col < 1.5 or col > 4.5:
        return None, None

    teff_gaia_jh = None
    if jmag is not None and hmag is not None:
        jh = jmag - hmag
        a, b, c, d, e, f, g = 3.172, -2.475, 1.082, -0.2231, 0.01738, 0.08776, 0.04355
        teff_gaia_jh = round(3500.0 * (a + b * col + c * col**2 + d * col**3 + e * col**4 + f * jh + g * jh**2), 1)

    a2, b2, c2, d2, e2 = 3.245, -2.4309, 1.043, -0.2127, 0.01649
    teff_gaia = round(3500.0 * (a2 + b2 * col + c2 * col**2 + d2 * col**3 + e2 * col**4), 1)
    return teff_gaia_jh, teff_gaia


def derive_fields(yaml_data: dict) -> None:
    """Populate derived fields; only fill values currently set to None."""
    ra = _pf(yaml_data.get("RA", {}).get("value"))
    dec = _pf(yaml_data.get("DEC", {}).get("value"))
    plx = _pf(yaml_data.get("PLX", {}).get("value"))
    pmra = _pf(yaml_data.get("PMRA", {}).get("value"))
    pmde = _pf(yaml_data.get("PMDE", {}).get("value"))
    rv = _pf(yaml_data.get("RV", {}).get("value"))

    gmag = _pf(yaml_data.get("G_MAG", {}).get("value"))
    gbp = _pf(yaml_data.get("GBP_MAG", {}).get("value"))
    grp = _pf(yaml_data.get("GRP_MAG", {}).get("value"))
    jmag = _pf(yaml_data.get("J_MAG", {}).get("value"))
    hmag = _pf(yaml_data.get("H_MAG", {}).get("value"))
    kmag = _pf(yaml_data.get("KS_MAG", {}).get("value"))
    w1 = _pf(yaml_data.get("W1_MAG", {}).get("value"))
    w2 = _pf(yaml_data.get("W2_MAG", {}).get("value"))

    # If EPOCH is provided and is not J2000, propagate stored RA/DEC back to J2000.
    # This is needed when coordinates are stored at a Gaia epoch (e.g. 2016.0).
    # J2000 = JD 2451545.0
    epoch_jd = _pf(yaml_data.get("EPOCH"))
    epoch_jyear = _jd_to_jyear(epoch_jd)
    ra_coord = ra
    dec_coord = dec
    ra_source = _nv(yaml_data.get("RA", {}).get("source"))
    is_gaia_epoch_coord = ra_source is not None and "GAIA" in ra_source.upper()
    if (ra is not None and dec is not None and pmra is not None and pmde is not None
            and is_gaia_epoch_coord and epoch_jyear is not None and abs(epoch_jd - 2451545.0) > 1.0):
        ra_coord, dec_coord = _propagate_to_epoch(
            ra_deg=ra,
            dec_deg=dec,
            pmra=pmra,
            pmdec=pmde,
            source_epoch_jyear=epoch_jyear,
            target_epoch_jyear=2000.0,
            plx_mas=plx,
            rv_kms=rv,
        )

    if yaml_data.get("RA_J2000_DEG") is None:
        yaml_data["RA_J2000_DEG"] = ra_coord
    if yaml_data.get("DEC_J2000_DEG") is None:
        yaml_data["DEC_J2000_DEG"] = dec_coord

    # Coordinate formats
    if yaml_data.get("RA_HMS") is None:
        yaml_data["RA_HMS"] = _format_ra_hms(ra_coord)
    if yaml_data.get("DEC_DMS") is None:
        yaml_data["DEC_DMS"] = _format_dec_dms(dec_coord)

    gal_l, gal_b = _galactic_from_radec(ra_coord, dec_coord)
    if yaml_data.get("GALACTIC_LON") is None:
        yaml_data["GALACTIC_LON"] = gal_l
    if yaml_data.get("GALACTIC_LAT") is None:
        yaml_data["GALACTIC_LAT"] = gal_b

    ecl_l, ecl_b = _ecliptic_from_radec(ra_coord, dec_coord)
    if yaml_data.get("ECLIPTIC_LON") is None:
        yaml_data["ECLIPTIC_LON"] = ecl_l
    if yaml_data.get("ECLIPTIC_LAT") is None:
        yaml_data["ECLIPTIC_LAT"] = ecl_b

    # Telluric RV limits and windows (vsys + vbary)
    if rv is not None and ecl_b is not None:
        vbary_amp = EARTH_ORBITAL_SPEED_KMS * math.cos(math.radians(ecl_b))
        if yaml_data.get("TELLURIC_VSYS_PLUS_VBARY_MIN") is None:
            yaml_data["TELLURIC_VSYS_PLUS_VBARY_MIN"] = rv - vbary_amp
        if yaml_data.get("TELLURIC_VSYS_PLUS_VBARY_MAX") is None:
            yaml_data["TELLURIC_VSYS_PLUS_VBARY_MAX"] = rv + vbary_amp
        if yaml_data.get("TELLURIC_LIMIT_WINDOWS") is None:
            yaml_data["TELLURIC_LIMIT_WINDOWS"] = _telluric_windows(ra_coord, dec_coord, rv)

    # Kinematics
    if plx is not None and plx > 0 and pmra is not None and pmde is not None:
        mu_tot = math.sqrt(pmra**2 + pmde**2)
        d_pc = 1000.0 / plx
        v_sky = 4.74047 * d_pc * mu_tot / 1000.0
        if yaml_data.get("V_SKY") is None:
            yaml_data["V_SKY"] = v_sky
        if rv is not None and yaml_data.get("V3D") is None:
            yaml_data["V3D"] = math.sqrt(v_sky**2 + rv**2)

    if all(v is not None for v in (ra_coord, dec_coord, plx, pmra, pmde, rv)) and plx > 0:
        ra_r = math.radians(ra_coord)
        dec_r = math.radians(dec_coord)
        d_pc = 1000.0 / plx
        k = 4.74047
        cos_ra = math.cos(ra_r)
        sin_ra = math.sin(ra_r)
        cos_dec = math.cos(dec_r)
        sin_dec = math.sin(dec_r)
        a_matrix = [
            [-sin_ra, -cos_ra * sin_dec, cos_ra * cos_dec],
            [cos_ra, -sin_ra * sin_dec, sin_ra * cos_dec],
            [0.0, cos_dec, sin_dec],
        ]
        velocity_components = [
            k * d_pc * pmra / 1000.0,
            k * d_pc * pmde / 1000.0,
            rv,
        ]
        v_eq = [sum(a_matrix[i][j] * velocity_components[j] for j in range(3)) for i in range(3)]
        u = sum(EQUATORIAL_TO_GALACTIC_MATRIX[0][j] * v_eq[j] for j in range(3))
        v = sum(EQUATORIAL_TO_GALACTIC_MATRIX[1][j] * v_eq[j] for j in range(3))
        w = sum(EQUATORIAL_TO_GALACTIC_MATRIX[2][j] * v_eq[j] for j in range(3))
        if yaml_data.get("U") is None:
            yaml_data["U"] = u
        if yaml_data.get("V") is None:
            yaml_data["V"] = v
        if yaml_data.get("W") is None:
            yaml_data["W"] = w

    # Absolute magnitudes
    if yaml_data.get("AMAG_G") is None:
        yaml_data["AMAG_G"] = _absolute_mag(gmag, plx)
    if yaml_data.get("AMAG_KS") is None:
        yaml_data["AMAG_KS"] = _absolute_mag(kmag, plx)

    # Gaia-color Teff relations
    teff_gaia_jh, teff_gaia = _teff_from_gaia_colors(gbp, grp, jmag, hmag)
    if yaml_data.get("TEFF_GAIA_JH") is None:
        yaml_data["TEFF_GAIA_JH"] = teff_gaia_jh
    if yaml_data.get("TEFF_GAIA") is None:
        yaml_data["TEFF_GAIA"] = teff_gaia

    # Duque-Arribas photometric [Fe/H]
    m_ks = _absolute_mag(kmag, plx)
    m_g = _absolute_mag(gmag, plx)
    spt = _nv(yaml_data.get("SPT", {}).get("value"))
    is_m_candidate = (m_g is not None and 7.5 <= m_g <= 16.5) or (spt is not None and spt.upper().startswith("M"))

    if is_m_candidate and all(v is not None for v in (m_ks, gbp, grp, w1, w2)):
        x = w1 - w2
        denom = 0.618 + 0.960 * x
        if abs(denom) > 1e-9:
            feh = ((gbp - grp) - 0.596 - 2.336 * x - 0.498 * (x ** 2) - 0.254 * m_ks) / denom
            yaml_data["FE_H"] = feh

    if yaml_data.get("FE_H") is None:
        feh_gaia = _pf(yaml_data.get("GAIA_MH_GSPPHOT"))
        if feh_gaia is not None:
            yaml_data["FE_H"] = feh_gaia

    # Mann+2015 / Delfosse+2000 radius+mass, then logg and luminosity
    if m_ks is not None and is_m_candidate:
        if yaml_data.get("R_STAR_MKS") is None:
            r_star = 1.9515 - 0.3520 * m_ks + 0.01680 * (m_ks ** 2)
            if r_star > 0:
                yaml_data["R_STAR_MKS"] = r_star

        if yaml_data.get("R_STAR_MKS_FEH") is None and yaml_data.get("FE_H") is not None:
            feh = _pf(yaml_data.get("FE_H"))
            if feh is not None:
                r_star_feh = 1.9305 - 0.3466 * m_ks + 0.01647 * (m_ks ** 2) + 0.04458 * feh
                if r_star_feh > 0:
                    yaml_data["R_STAR_MKS_FEH"] = r_star_feh

        if yaml_data.get("MASS_STAR_MANN15") is None:
            m_star = 0.5858 + 0.3872 * m_ks - 0.1217 * (m_ks ** 2) + 0.0106 * (m_ks ** 3) - 2.7262e-4 * (m_ks ** 4)
            if m_star > 0:
                yaml_data["MASS_STAR_MANN15"] = m_star

        if yaml_data.get("MASS_STAR_DELFOSSE00") is None and 4.5 <= m_ks <= 9.5:
            log_m_del = 1e-3 * (1.8 + 6.12 * m_ks + 13.205 * m_ks**2 - 6.2315 * m_ks**3 + 0.37529 * m_ks**4)
            yaml_data["MASS_STAR_DELFOSSE00"] = 10.0 ** log_m_del

    if not is_m_candidate:
        feh_gaia = _pf(yaml_data.get("GAIA_MH_GSPPHOT"))
        yaml_data["FE_H"] = feh_gaia
        for key in ("R_STAR_MKS", "R_STAR_MKS_FEH", "MASS_STAR_MANN15", "MASS_STAR_DELFOSSE00"):
            if yaml_data.get(key) is not None:
                yaml_data[key] = None

    mass_for_logg = _pf(yaml_data.get("MASS_STAR_MANN15"))
    radius_for_logg = _pf(yaml_data.get("R_STAR_MKS"))
    if not is_m_candidate:
        if mass_for_logg is None:
            mass_for_logg = _pf(yaml_data.get("GAIA_MASS_FLAME"))
        if radius_for_logg is None:
            radius_for_logg = _pf(yaml_data.get("GAIA_RADIUS_FLAME"))
    else:
        if radius_for_logg is None:
            radius_for_logg = _pf(yaml_data.get("R_STAR_MKS_FEH"))
    if yaml_data.get("LOG_G") is None and mass_for_logg is not None and radius_for_logg is not None and mass_for_logg > 0 and radius_for_logg > 0:
        yaml_data["LOG_G"] = 4.438 + math.log10(mass_for_logg) - 2.0 * math.log10(radius_for_logg)

    teff_for_l = _pf(yaml_data.get("TEFF_GAIA_JH"))
    if teff_for_l is None:
        teff_for_l = _pf(yaml_data.get("TEFF_GAIA"))
    if teff_for_l is None:
        teff_for_l = _pf(yaml_data.get("TEFF", {}).get("value"))
    if is_m_candidate:
        radius_for_l = _pf(yaml_data.get("R_STAR_MKS"))
        if radius_for_l is None:
            radius_for_l = radius_for_logg
        if radius_for_l is not None and teff_for_l is not None and teff_for_l > 0:
            yaml_data["L_STAR"] = (radius_for_l ** 2) * (teff_for_l / 5778.0) ** 4
    else:
        lum_gaia = _pf(yaml_data.get("GAIA_LUM_FLAME"))
        if lum_gaia is not None:
            yaml_data["L_STAR"] = lum_gaia

    if not is_m_candidate:
        mass_gaia = _pf(yaml_data.get("GAIA_MASS_FLAME"))
        radius_gaia = _pf(yaml_data.get("GAIA_RADIUS_FLAME"))
        lum_gaia = _pf(yaml_data.get("GAIA_LUM_FLAME"))
        if mass_gaia is not None and radius_gaia is not None and mass_gaia > 0 and radius_gaia > 0:
            yaml_data["LOG_G"] = 4.438 + math.log10(mass_gaia) - 2.0 * math.log10(radius_gaia)
        if lum_gaia is not None:
            yaml_data["L_STAR"] = lum_gaia


# ---------------------------------------------------------------------------
# SIMBAD / Gaia / VizieR resolution  (I/O-free, call directly for new targets)
# ---------------------------------------------------------------------------

def _escape_simbad(name: str) -> str:
    return name.replace("'", "''")


def _extract_wisea_designation(identifier: Any) -> str | None:
    sval = _nv(identifier)
    if sval is None:
        return None
    if sval.startswith("WISEA "):
        return sval.replace("WISEA ", "", 1).strip()
    return None


def _fetch_wise_by_designation(designation: str) -> tuple[Any, ...] | None:
    safe = designation.replace("'", "''")
    adql = f'''
    SELECT TOP 1 "AllWISE", W1mag, e_W1mag, W2mag, e_W2mag, W3mag, e_W3mag, W4mag, e_W4mag, RAJ2000, DEJ2000
    FROM "II/328/allwise"
    WHERE "AllWISE" = '{safe}'
    '''
    payload = _vizier_json(adql)
    if payload is None:
        return None
    rows = payload.get("data", [])
    if not rows:
        return None
    row = rows[0]
    if len(row) < 11:
        return None
    return tuple(row)


def resolve_from_name(name: str) -> dict | None:
    """
    Query SIMBAD TAP for *name* and return a normalised result dict.

    All keys may be None.  Keys:
      simbad_main_id, ra_deg, dec_deg, parallax_mas, pmra_masyr, pmdec_masyr,
      rv_kms, sp_type, otype_txt,
      G_mag, J_mag, H_mag, K_mag,
      GBP_mag, GBP_err, GRP_mag, GRP_err,
      W1_mag, W1_err, W2_mag, W2_err, W3_mag, W3_err, W4_mag, W4_err,
      gaia_source_id,
      gaia_teff_gspphot, gaia_logg_gspphot, gaia_mh_gspphot,
      gaia_radius_flame, gaia_lum_flame, gaia_mass_flame,
      teff_simbad,
    """
    safe = _escape_simbad(name)

    # ------------------------------------------------------------------
    # 1. Core SIMBAD query
    # ------------------------------------------------------------------
    adql = f"""
    SELECT TOP 1
      b.main_id,
      b.ra,
      b.dec,
      b.plx_value,
      b.pmra,
      b.pmdec,
      b.rvz_radvel,
      b.sp_type,
      b.otype_txt,
      f.G,
      f.J,
      f.H,
      f.K
    FROM ident i
    JOIN basic b ON i.oidref = b.oid
    LEFT JOIN allfluxes f ON b.oid = f.oidref
    WHERE i.id = '{safe}'
    """
    payload = _simbad_json(adql)
    if payload is None:
        return None
    rows = payload.get("data", [])
    if not rows:
        return None

    row = rows[0]
    keys = ["simbad_main_id", "ra_deg", "dec_deg", "parallax_mas",
            "pmra_masyr", "pmdec_masyr", "rv_kms",
            "sp_type", "otype_txt",
            "G_mag", "J_mag", "H_mag", "K_mag"]
    result: dict[str, Any] = {k: _nv(v) for k, v in zip(keys, row)}

    # Null-fill all optional keys
    for k in ("GBP_mag", "GBP_err", "GRP_mag", "GRP_err",
               "W1_mag", "W1_err", "W2_mag", "W2_err",
               "W3_mag", "W3_err", "W4_mag", "W4_err",
               "gaia_source_id",
               "gaia_teff_gspphot", "gaia_logg_gspphot", "gaia_mh_gspphot",
               "gaia_radius_flame", "gaia_lum_flame", "gaia_mass_flame",
               "teff_simbad"):
        result.setdefault(k, None)

    # ------------------------------------------------------------------
    # 2. SIMBAD Teff measurements (median of available values)
    # ------------------------------------------------------------------
    adql_teff = f"""
    SELECT TOP 10 m.teff
    FROM mesFe_h m
    JOIN ident i ON m.oidref = i.oidref
    WHERE i.id = '{safe}'
      AND m.teff IS NOT NULL
    """
    teff_p = _simbad_json(adql_teff)
    if teff_p:
        teffs = [_pf(r[0]) for r in teff_p.get("data", []) if _pf(r[0]) is not None]
        if teffs:
            result["teff_simbad"] = sorted(teffs)[len(teffs) // 2]

    # ------------------------------------------------------------------
    # 3. Gaia identifier via SIMBAD cross-match
    # ------------------------------------------------------------------
    adql_idents = f"""
    SELECT i2.id
    FROM ident i
    JOIN ident i2 ON i.oidref = i2.oidref
    WHERE i.id = '{safe}'
      AND (
           i2.id LIKE 'Gaia DR3 %'
        OR i2.id LIKE 'Gaia EDR3 %'
        OR i2.id LIKE 'Gaia DR2 %'
        OR i2.id LIKE 'WISEA %'
      )
    """
    idents_p = _simbad_json(adql_idents)
    gaia_release: str | None = None
    gaia_source_id: str | None = None
    wise_designation: str | None = None
    if idents_p:
        for irow in idents_p.get("data", []):
            ident_val = _nv(irow[0]) if irow else None
            if ident_val is None:
                continue
            for prefix, rel in (("Gaia DR3 ", "dr3"), ("Gaia EDR3 ", "edr3"), ("Gaia DR2 ", "dr2")):
                if ident_val.startswith(prefix):
                    candidate = ident_val[len(prefix):].strip()
                    if candidate.isdigit():
                        gaia_source_id = candidate
                        gaia_release = rel
                        break
            if wise_designation is None:
                wise_designation = _extract_wisea_designation(ident_val)
            if gaia_source_id is not None and wise_designation is not None:
                break

    result["gaia_source_id"] = gaia_source_id

    # ------------------------------------------------------------------
    # 4. Gaia BP/RP + proper motions from VizieR
    # ------------------------------------------------------------------
    ra_deg = _pf(result["ra_deg"])
    dec_deg = _pf(result["dec_deg"])
    simbad_ra_deg = ra_deg
    simbad_dec_deg = dec_deg
    simbad_pmra = _pf(result.get("pmra_masyr"))
    simbad_pmdec = _pf(result.get("pmdec_masyr"))
    simbad_plx = _pf(result.get("parallax_mas"))
    simbad_rv = _pf(result.get("rv_kms"))
    gaia_ref_epoch = None

    if gaia_source_id is not None:
        if gaia_release == "dr2":
            adql_gaia = f"""
            SELECT TOP 1 phot_bp_mean_mag, phot_rp_mean_mag, phot_bp_mean_mag_error, phot_rp_mean_mag_error,
                         ra, dec, pmra, pmdec
            FROM "I/345/gaia2"
            WHERE source_id = {gaia_source_id}
            """
        else:
            table_map = {"dr3": '"I/355/gaiadr3"', "edr3": '"I/350/gaiaedr3"'}
            table = table_map.get(gaia_release or "dr3", '"I/355/gaiadr3"')
            adql_gaia = f"""
            SELECT TOP 1 BPmag, RPmag, e_BPmag, e_RPmag, RA_ICRS, DE_ICRS, pmRA, pmDE
            FROM {table}
            WHERE Source = {gaia_source_id}
            """
        gaia_p = _vizier_json(adql_gaia)
        if gaia_p:
            gdata = gaia_p.get("data", [])
            if gdata:
                gv = gdata[0]
                _set_if_none(result, "GBP_mag", _nv(gv[0]) if len(gv) > 0 else None)
                _set_if_none(result, "GRP_mag", _nv(gv[1]) if len(gv) > 1 else None)
                _set_if_none(result, "GBP_err", _nv(gv[2]) if len(gv) > 2 else None)
                _set_if_none(result, "GRP_err", _nv(gv[3]) if len(gv) > 3 else None)
                # Prefer Gaia astrometric position as the epoch anchor for WISE propagation.
                ra_gaia = _pf(gv[4]) if len(gv) > 4 else None
                dec_gaia = _pf(gv[5]) if len(gv) > 5 else None
                if ra_gaia is not None:
                    ra_deg = ra_gaia
                if dec_gaia is not None:
                    dec_deg = dec_gaia
                if gaia_release in {"dr3", "edr3"}:
                    gaia_ref_epoch = 2016.0
                elif gaia_release == "dr2":
                    gaia_ref_epoch = 2015.5
                _set_if_none(result, "pmra_masyr", _nv(gv[6]) if len(gv) > 6 else None)
                _set_if_none(result, "pmdec_masyr", _nv(gv[7]) if len(gv) > 7 else None)

    # Positional Gaia fallback for BP/RP (in case source ID lookup failed)
    gaia_query_ra = ra_deg
    gaia_query_dec = dec_deg
    if gaia_ref_epoch is None and simbad_ra_deg is not None and simbad_dec_deg is not None:
        gaia_query_ra, gaia_query_dec = _propagate_to_epoch(
            ra_deg=simbad_ra_deg,
            dec_deg=simbad_dec_deg,
            pmra=simbad_pmra,
            pmdec=simbad_pmdec,
            source_epoch_jyear=2000.0,
            target_epoch_jyear=2016.0,
            plx_mas=simbad_plx,
            rv_kms=simbad_rv,
        )

    if (result.get("GBP_mag") is None or result.get("GRP_mag") is None) and gaia_query_ra is not None and gaia_query_dec is not None:
        adql_pos = f"""
        SELECT TOP 10 RA_ICRS, DE_ICRS, BPmag, RPmag, e_BPmag, e_RPmag
        FROM "I/355/gaiadr3"
        WHERE 1 = CONTAINS(
          POINT('ICRS', RA_ICRS, DE_ICRS),
          CIRCLE('ICRS', {gaia_query_ra}, {gaia_query_dec}, 0.0005)
        )
        """
        gaia_pos_p = _vizier_json(adql_pos)
        if gaia_pos_p:
            best_gaia: tuple | None = None
            for gv in gaia_pos_p.get("data", []):
                if len(gv) < 4:
                    continue
                rra, rdec = _pf(gv[0]), _pf(gv[1])
                if rra is None or rdec is None:
                    continue
                dist = _sep_arcsec(gaia_query_ra, gaia_query_dec, rra, rdec)
                if best_gaia is None or dist < best_gaia[0]:
                    best_gaia = (dist, gv)
            if best_gaia:
                gv = best_gaia[1]
                # Positional Gaia fallback uses DR3 table.
                ra_gaia = _pf(gv[0]) if len(gv) > 0 else None
                dec_gaia = _pf(gv[1]) if len(gv) > 1 else None
                if ra_gaia is not None:
                    ra_deg = ra_gaia
                if dec_gaia is not None:
                    dec_deg = dec_gaia
                gaia_ref_epoch = 2016.0
                _set_if_none(result, "GBP_mag", _nv(gv[2]) if len(gv) > 2 else None)
                _set_if_none(result, "GRP_mag", _nv(gv[3]) if len(gv) > 3 else None)
                _set_if_none(result, "GBP_err", _nv(gv[4]) if len(gv) > 4 else None)
                _set_if_none(result, "GRP_err", _nv(gv[5]) if len(gv) > 5 else None)

    # ------------------------------------------------------------------
    # 5. Gaia DR3 astrophysical parameters (GSP-Phot / FLAME)
    # ------------------------------------------------------------------
    if gaia_source_id is not None:
        adql_astro = f"""
        SELECT TOP 1 teff_gspphot, logg_gspphot, mh_gspphot, radius_flame, lum_flame, mass_flame
        FROM gaiadr3.astrophysical_parameters
        WHERE source_id = {gaia_source_id}
        """
        astro_rows = _gaia_csv(adql_astro)
        if astro_rows:
            ar = astro_rows[0]
            _set_if_none(result, "gaia_teff_gspphot", _nv(ar.get("teff_gspphot")))
            _set_if_none(result, "gaia_logg_gspphot", _nv(ar.get("logg_gspphot")))
            _set_if_none(result, "gaia_mh_gspphot",   _nv(ar.get("mh_gspphot")))
            _set_if_none(result, "gaia_radius_flame",  _nv(ar.get("radius_flame")))
            _set_if_none(result, "gaia_lum_flame",     _nv(ar.get("lum_flame")))
            _set_if_none(result, "gaia_mass_flame",    _nv(ar.get("mass_flame")))

    # ------------------------------------------------------------------
    # 6. AllWISE photometry (VizieR positional xmatch, proper-motion-aware)
    # ------------------------------------------------------------------
    pmra_v  = _pf(result.get("pmra_masyr"))
    pmdec_v = _pf(result.get("pmdec_masyr"))

    if wise_designation is not None:
        wise_row = _fetch_wise_by_designation(wise_designation)
        if wise_row is not None:
            _set_if_none(result, "W1_mag", _nv(wise_row[1]))
            _set_if_none(result, "W1_err", _nv(wise_row[2]))
            _set_if_none(result, "W2_mag", _nv(wise_row[3]))
            _set_if_none(result, "W2_err", _nv(wise_row[4]))
            _set_if_none(result, "W3_mag", _nv(wise_row[5]))
            _set_if_none(result, "W3_err", _nv(wise_row[6]))
            _set_if_none(result, "W4_mag", _nv(wise_row[7]))
            _set_if_none(result, "W4_err", _nv(wise_row[8]))

    if (result.get("W1_mag") is None or result.get("W2_mag") is None) and (ra_deg is not None and dec_deg is not None):
        if gaia_ref_epoch is not None and pmra_v is not None and pmdec_v is not None:
            # Propagate Gaia coordinates at their native reference epoch to AllWISE epoch.
            dt_yr = 2010.5 - gaia_ref_epoch
            mu_tot = math.sqrt(pmra_v ** 2 + pmdec_v ** 2)
            radius = max(20.0, min(300.0, abs(dt_yr) * mu_tot / 1000.0 + 12.0))
            ra_wise, dec_wise = _propagate_to_epoch(
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                pmra=pmra_v,
                pmdec=pmdec_v,
                source_epoch_jyear=gaia_ref_epoch,
                target_epoch_jyear=2010.5,
                plx_mas=simbad_plx,
                rv_kms=simbad_rv,
            )
        else:
            # Fallback: propagate SIMBAD J2000 coordinates to AllWISE epoch when PM is available.
            if simbad_ra_deg is not None and simbad_dec_deg is not None:
                ra_wise, dec_wise = _propagate_to_epoch(
                    ra_deg=simbad_ra_deg,
                    dec_deg=simbad_dec_deg,
                    pmra=simbad_pmra,
                    pmdec=simbad_pmdec,
                    source_epoch_jyear=2000.0,
                    target_epoch_jyear=2010.5,
                    plx_mas=simbad_plx,
                    rv_kms=simbad_rv,
                )
                if simbad_pmra is not None and simbad_pmdec is not None:
                    mu_tot = math.sqrt(simbad_pmra ** 2 + simbad_pmdec ** 2)
                    radius = max(20.0, min(300.0, abs(2010.5 - 2000.0) * mu_tot / 1000.0 + 12.0))
                else:
                    radius = 30.0
            else:
                radius = 30.0
                ra_wise, dec_wise = ra_deg, dec_deg

        adql_wise = f"""
        SELECT TOP 10 RAJ2000, DEJ2000, W1mag, e_W1mag, W2mag, e_W2mag, W3mag, e_W3mag, W4mag, e_W4mag
        FROM "II/328/allwise"
        WHERE 1 = CONTAINS(
          POINT('ICRS', RAJ2000, DEJ2000),
          CIRCLE('ICRS', {ra_wise}, {dec_wise}, {radius}/3600.0)
        )
        """
        wise_p = _vizier_json(adql_wise)
        if wise_p:
            best_wise: tuple | None = None
            for wv in wise_p.get("data", []):
                if len(wv) < 6:
                    continue
                wra, wdec = _pf(wv[0]), _pf(wv[1])
                if wra is None or wdec is None:
                    continue
                dist = _sep_arcsec(ra_wise, dec_wise, wra, wdec)
                if best_wise is None or dist < best_wise[0]:
                    best_wise = (dist, wv)
            if best_wise:
                wv = best_wise[1]
                _set_if_none(result, "W1_mag", _nv(wv[2]) if len(wv) > 2 else None)
                _set_if_none(result, "W1_err", _nv(wv[3]) if len(wv) > 3 else None)
                _set_if_none(result, "W2_mag", _nv(wv[4]) if len(wv) > 4 else None)
                _set_if_none(result, "W2_err", _nv(wv[5]) if len(wv) > 5 else None)
                _set_if_none(result, "W3_mag", _nv(wv[6]) if len(wv) > 6 else None)
                _set_if_none(result, "W3_err", _nv(wv[7]) if len(wv) > 7 else None)
                _set_if_none(result, "W4_mag", _nv(wv[8]) if len(wv) > 8 else None)
                _set_if_none(result, "W4_err", _nv(wv[9]) if len(wv) > 9 else None)

    return result


# ---------------------------------------------------------------------------
# YAML I/O  (kept separate so resolve_from_name can be used standalone)
# ---------------------------------------------------------------------------

def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as fh:
        return yaml.safe_load(fh) or {}


def save_yaml(path: Path, payload: dict) -> None:
    text = yaml.safe_dump(payload, sort_keys=False, allow_unicode=False)
    path.write_text(text, encoding="utf-8")


def _full_schema(d: dict) -> dict:
    """Ensure every expected top-level key exists (null if absent)."""
    for top in ("APERO_NAME", "ORIGINAL_NAME", "SIMBAD_NAME", "APERO_CLASS", "EPOCH"):
        d.setdefault(top, None)
    for block, units in [
        ("RA", "deg"), ("DEC", "deg"),
        ("PMRA", "mas/yr"), ("PMDE", "mas/yr"),
        ("PLX", "mas"), ("RV", "km/s"),
        ("TEFF", "K"),
    ]:
        sub = d.setdefault(block, {})
        sub.setdefault("value", None)
        sub.setdefault("source", None)
        sub.setdefault("units", units)
    for block in ("SPT",):
        sub = d.setdefault(block, {})
        sub.setdefault("value", None)
        sub.setdefault("source", None)
    vsini = d.setdefault("VSINI", {})
    for k in ("value", "err", "source"):
        vsini.setdefault(k, None)
    vsini.setdefault("units", "km/s")
    for block in ("G_MAG", "GBP_MAG", "GRP_MAG",
                  "J_MAG", "H_MAG", "KS_MAG",
                  "W1_MAG", "W2_MAG", "W3_MAG", "W4_MAG"):
        sub = d.setdefault(block, {})
        sub.setdefault("value", None)
        sub.setdefault("source", None)
    for k in ("KEYWORDS", "ALIASES"):
        d.setdefault(k, None)
    d.setdefault("NOTES", None)
    # Gaia-derived extra fields
    for k in ("GAIA_SOURCE_ID",
               "GAIA_TEFF_GSPPHOT", "GAIA_LOGG_GSPPHOT", "GAIA_MH_GSPPHOT",
               "GAIA_RADIUS_FLAME", "GAIA_LUM_FLAME", "GAIA_MASS_FLAME"):
        d.setdefault(k, None)
    for k in (
        "RA_HMS", "DEC_DMS", "RA_J2000_DEG", "DEC_J2000_DEG",
        "GALACTIC_LON", "GALACTIC_LAT", "ECLIPTIC_LON", "ECLIPTIC_LAT",
        "TELLURIC_VSYS_PLUS_VBARY_MIN", "TELLURIC_VSYS_PLUS_VBARY_MAX", "TELLURIC_LIMIT_WINDOWS",
        "V_SKY", "V3D", "U", "V", "W",
        "AMAG_G", "AMAG_KS",
        "TEFF_GAIA_JH", "TEFF_GAIA",
        "FE_H", "R_STAR_MKS", "R_STAR_MKS_FEH", "MASS_STAR_MANN15",
        "MASS_STAR_DELFOSSE00", "LOG_G", "L_STAR",
    ):
        d.setdefault(k, None)
    return d


def _merge_scalar(yaml_data: dict, yaml_block: str,
                  simbad: dict, simbad_key: str, source: str) -> None:
    block = yaml_data[yaml_block]
    if block.get("value") is None:
        val = _pf(simbad.get(simbad_key))
        if val is not None:
            block["value"] = val
            block["source"] = source


def _merge_mag(yaml_data: dict, yaml_block: str,
               simbad: dict, simbad_key: str, source: str) -> None:
    block = yaml_data.get(yaml_block)
    if not isinstance(block, dict):
        yaml_data[yaml_block] = {"value": None, "source": None}
        block = yaml_data[yaml_block]
    if block.get("value") is None:
        val = _pf(simbad.get(simbad_key))
        if val is not None:
            block["value"] = val
            block["source"] = source


def update_yaml_from_simbad(yaml_data: dict, simbad: dict) -> dict:
    """
    Merge a resolve_from_name() result into yaml_data, filling only None values.
    Modifies yaml_data in-place and returns it.
    """
    _full_schema(yaml_data)

    # SIMBAD_NAME
    if yaml_data.get("SIMBAD_NAME") is None:
        yaml_data["SIMBAD_NAME"] = simbad.get("simbad_main_id")

    # Astrometry
    _merge_scalar(yaml_data, "RA",   simbad, "ra_deg",       "SIMBAD")
    _merge_scalar(yaml_data, "DEC",  simbad, "dec_deg",       "SIMBAD")
    _merge_scalar(yaml_data, "PLX",  simbad, "parallax_mas",  "SIMBAD")
    _merge_scalar(yaml_data, "PMRA", simbad, "pmra_masyr",    "SIMBAD")
    _merge_scalar(yaml_data, "PMDE", simbad, "pmdec_masyr",   "SIMBAD")
    _merge_scalar(yaml_data, "RV",   simbad, "rv_kms",        "SIMBAD")

    # Teff – prefer Gaia GSP-Phot, then SIMBAD mes_teff
    if yaml_data["TEFF"]["value"] is None:
        gaia_teff  = _pf(simbad.get("gaia_teff_gspphot"))
        simbad_teff = _pf(simbad.get("teff_simbad"))
        if gaia_teff is not None:
            yaml_data["TEFF"]["value"]  = gaia_teff
            yaml_data["TEFF"]["source"] = "Gaia DR3 GSP-Phot"
        elif simbad_teff is not None:
            yaml_data["TEFF"]["value"]  = simbad_teff
            yaml_data["TEFF"]["source"] = "SIMBAD mes_teff"

    # SpT
    if yaml_data["SPT"]["value"] is None and _nv(simbad.get("sp_type")) is not None:
        yaml_data["SPT"]["value"]  = simbad["sp_type"]
        yaml_data["SPT"]["source"] = "SIMBAD"

    # Photometry
    _merge_mag(yaml_data, "G_MAG",   simbad, "G_mag",  "SIMBAD/Gaia")
    _merge_mag(yaml_data, "J_MAG",   simbad, "J_mag",  "SIMBAD/2MASS")
    _merge_mag(yaml_data, "H_MAG",   simbad, "H_mag",  "SIMBAD/2MASS")
    _merge_mag(yaml_data, "KS_MAG",  simbad, "K_mag",  "SIMBAD/2MASS")
    _merge_mag(yaml_data, "GBP_MAG", simbad, "GBP_mag","Gaia VizieR")
    _merge_mag(yaml_data, "GRP_MAG", simbad, "GRP_mag","Gaia VizieR")
    _merge_mag(yaml_data, "W1_MAG",  simbad, "W1_mag", "AllWISE")
    _merge_mag(yaml_data, "W2_MAG",  simbad, "W2_mag", "AllWISE")
    _merge_mag(yaml_data, "W3_MAG",  simbad, "W3_mag", "AllWISE")
    _merge_mag(yaml_data, "W4_MAG",  simbad, "W4_mag", "AllWISE")

    # Gaia extra scalar fields
    for yaml_key, src_key in [
        ("GAIA_TEFF_GSPPHOT", "gaia_teff_gspphot"),
        ("GAIA_LOGG_GSPPHOT", "gaia_logg_gspphot"),
        ("GAIA_MH_GSPPHOT",   "gaia_mh_gspphot"),
        ("GAIA_RADIUS_FLAME", "gaia_radius_flame"),
        ("GAIA_LUM_FLAME",    "gaia_lum_flame"),
        ("GAIA_MASS_FLAME",   "gaia_mass_flame"),
    ]:
        if yaml_data.get(yaml_key) is None and _nv(simbad.get(src_key)) is not None:
            v = _pf(simbad[src_key])
            yaml_data[yaml_key] = v if v is not None else simbad[src_key]

    # Keep Gaia source ID as string to avoid scientific notation.
    if yaml_data.get("GAIA_SOURCE_ID") is None and _nv(simbad.get("gaia_source_id")) is not None:
        yaml_data["GAIA_SOURCE_ID"] = str(simbad.get("gaia_source_id"))

    derive_fields(yaml_data)

    return yaml_data


# ---------------------------------------------------------------------------
# Main driver – iterate YAML files
# ---------------------------------------------------------------------------

def _best_search_name(yaml_data: dict) -> str | None:
    for key in ("SIMBAD_NAME", "ORIGINAL_NAME", "APERO_NAME"):
        v = _nv(yaml_data.get(key))
        if v is not None:
            return v
    return None


# Fetchable fields: (yaml_key_or_block, sub_key_or_None)
# If sub_key is None the top-level value is tested directly.
_FETCHABLE_FIELDS: list[tuple[str, str | None]] = [
    # Core astrometry
    ("RA",   "value"),
    ("DEC",  "value"),
    ("PLX",  "value"),
    ("PMRA", "value"),
    ("PMDE", "value"),
    ("RV",   "value"),
    # Stellar parameters from SIMBAD
    ("TEFF", "value"),
    ("SPT",  "value"),
    # Photometry
    ("G_MAG",   "value"),
    ("J_MAG",   "value"),
    ("H_MAG",   "value"),
    ("KS_MAG",  "value"),
    ("GBP_MAG", "value"),
    ("GRP_MAG", "value"),
    ("W1_MAG",  "value"),
    ("W2_MAG",  "value"),
    ("W3_MAG",  "value"),
    ("W4_MAG",  "value"),
    # Gaia extras
    ("GAIA_SOURCE_ID",    None),
    ("GAIA_TEFF_GSPPHOT", None),
    ("GAIA_LOGG_GSPPHOT", None),
    ("GAIA_MH_GSPPHOT",   None),
    ("GAIA_RADIUS_FLAME", None),
    ("GAIA_LUM_FLAME",    None),
    ("GAIA_MASS_FLAME",   None),
]

_DERIVED_FIELDS: tuple[str, ...] = (
    "RA_HMS", "DEC_DMS", "RA_J2000_DEG", "DEC_J2000_DEG",
    "GALACTIC_LON", "GALACTIC_LAT", "ECLIPTIC_LON", "ECLIPTIC_LAT",
    "TELLURIC_VSYS_PLUS_VBARY_MIN", "TELLURIC_VSYS_PLUS_VBARY_MAX", "TELLURIC_LIMIT_WINDOWS",
    "V_SKY", "V3D", "U", "V", "W",
    "AMAG_G", "AMAG_KS",
    "TEFF_GAIA_JH", "TEFF_GAIA",
    "FE_H", "R_STAR_MKS", "R_STAR_MKS_FEH", "MASS_STAR_MANN15",
    "MASS_STAR_DELFOSSE00", "LOG_G", "L_STAR",
)


def _missing_fetchable_fields(yaml_data: dict) -> list[str]:
    """Return list of field names whose value is still None."""
    missing = []
    for block, sub in _FETCHABLE_FIELDS:
        if sub is None:
            val = yaml_data.get(block)
        else:
            val = yaml_data.get(block, {}).get(sub) if isinstance(yaml_data.get(block), dict) else None
        if _nv(val) is None:
            missing.append(block if sub is None else f"{block}.{sub}")
    return missing


def _clear_fetchable_fields(yaml_data: dict) -> None:
    """Clear values in fetchable fields so update_yaml_from_simbad can repopulate them."""
    for block, sub in _FETCHABLE_FIELDS:
        if sub is None:
            yaml_data[block] = None
            continue
        if not isinstance(yaml_data.get(block), dict):
            yaml_data[block] = {"value": None, "source": None}
        else:
            yaml_data[block][sub] = None
            if "source" in yaml_data[block]:
                yaml_data[block]["source"] = None


def _clear_derived_fields(yaml_data: dict) -> None:
    for key in _DERIVED_FIELDS:
        yaml_data[key] = None


def run(astrometrics_dir: Path,
        overwrite_existing: bool = False,
        dry_run: bool = False,
        limit: int | None = None,
        delay_s: float = 0.3,
        single: str | None = None) -> None:
    yaml_files = sorted(astrometrics_dir.glob("*.yaml"))
    if single is not None:
        yaml_files = [f for f in yaml_files
                      if f.stem == single or f.name == single]

    if limit is not None:
        yaml_files = yaml_files[:limit]

    total = len(yaml_files)
    print(f"Processing {total} YAML files in {astrometrics_dir}.")

    resolved = 0
    failed = 0
    skipped = 0

    for idx, yaml_path in enumerate(yaml_files, start=1):
        yaml_data = load_yaml(yaml_path)
        _full_schema(yaml_data)

        name = _best_search_name(yaml_data)
        if name is None:
            print(f"[{idx}/{total}] {yaml_path.name}: no resolvable name – skip.")
            skipped += 1
            continue

        already_resolved = yaml_data.get("SIMBAD_NAME") is not None and not overwrite_existing
        missing = _missing_fetchable_fields(yaml_data)

        if already_resolved and not missing:
            # Everything fetchable is present – only recompute derived quantities.
            _clear_derived_fields(yaml_data)
            derive_fields(yaml_data)
            if not dry_run:
                save_yaml(yaml_path, yaml_data)
            print(f"[{idx}/{total}] {yaml_path.name}: complete – derived fields refreshed.")
            skipped += 1
            continue

        if already_resolved and missing:
            # Has SIMBAD_NAME but some fetchable fields are missing – partial refetch.
            print(
                f"[{idx}/{total}] {yaml_path.name}: partial – missing {len(missing)} field(s): "
                f"{', '.join(missing[:6])}{'…' if len(missing) > 6 else ''} – refetching...",
                end=" ", flush=True,
            )
            simbad = resolve_from_name(name)
            if simbad is not None:
                _clear_derived_fields(yaml_data)
                update_yaml_from_simbad(yaml_data, simbad)
                resolved += 1
                print("OK")
            else:
                # Still run derive_fields with what we have
                _clear_derived_fields(yaml_data)
                derive_fields(yaml_data)
                failed += 1
                print("refetch failed – derived fields refreshed.")
            if not dry_run:
                save_yaml(yaml_path, yaml_data)
            if delay_s > 0 and idx < total:
                time.sleep(delay_s)
            continue

        # Not yet resolved at all – full resolution.
        print(f"[{idx}/{total}] {yaml_path.name}: resolving '{name}'...", end=" ", flush=True)

        simbad = resolve_from_name(name)

        # Alias fallback
        if simbad is None:
            aliases = yaml_data.get("ALIASES") or []
            for alias in (aliases or [])[:5]:
                a = _nv(alias)
                if a and a != name:
                    simbad = resolve_from_name(a)
                    if simbad is not None:
                        print(f"(via alias '{a}')", end=" ", flush=True)
                        break

        if simbad is None:
            print("NOT FOUND.")
            yaml_data["SIMBAD_NAME"] = None
            failed += 1
        else:
            if overwrite_existing:
                _clear_fetchable_fields(yaml_data)
            _clear_derived_fields(yaml_data)
            update_yaml_from_simbad(yaml_data, simbad)
            resolved += 1
            print(f"OK → {yaml_data.get('SIMBAD_NAME') or '?'}")

        if not dry_run:
            save_yaml(yaml_path, yaml_data)

        if delay_s > 0 and idx < total:
            time.sleep(delay_s)

    print(f"\nDone: {resolved} resolved/refetched, {failed} failed, {skipped} already-complete.")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dir",
        default=str(script_dir / "astrometrics"),
        help="Directory containing OBJNAME.yaml files (default: ./astrometrics)",
    )
    parser.add_argument(
        "--overwrite", action="store_true",
        help="Re-query SIMBAD even if SIMBAD_NAME is already set",
    )
    parser.add_argument("--dry-run", action="store_true", help="Do not write YAML files")
    parser.add_argument("--limit", type=int, default=None, help="Process only first N files")
    parser.add_argument(
        "--delay", type=float, default=0.3,
        help="Seconds to sleep between TAP queries (default 0.3)",
    )
    parser.add_argument(
        "--name", default=None,
        help="Process only the YAML whose stem or filename matches this value",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    adir = Path(args.dir).expanduser().resolve()
    if not adir.is_dir():
        print(f"ERROR: {adir} is not a directory.", file=sys.stderr)
        sys.exit(1)
    run(
        astrometrics_dir=adir,
        overwrite_existing=args.overwrite,
        dry_run=args.dry_run,
        limit=args.limit,
        delay_s=args.delay,
        single=args.name,
    )


if __name__ == "__main__":
    try:
        _ = __file__
    except NameError:
        __file__ = os.path.join(os.getcwd(), "resolve_against_simbad.py")
    main()


