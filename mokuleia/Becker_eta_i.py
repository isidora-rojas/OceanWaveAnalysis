import numpy as np
import pandas as pd
import xarray as xr

def compute_eta_i(
    ds_eta_f: xr.Dataset,
    h_i,
    h_i_time,
    *,
    align: str = "nearest",            # "nearest" or "interp"
    tolerance: str | None = "10min",   # only for nearest
    force_positive_depth: bool = True,
    return_new_ds: bool = True,
) -> xr.Dataset | xr.DataArray:
    """
    Wave setup at reef flat based on equation (3) in Becker (2014)
    -----------------------------------------------------------

    eta_i = h_i - h_f + eta_f

    h: water depth (m)
    t: time
    eta_f; setdown at reef face


    This module computs reef setup from an existing dataset with eta_f and h_f. This is
    computed using Becker_eta_f.py to get eta_f from surface elevation spectra. 


    Parameters
    -----------
    ds_eta_f : xr.Dataset
        Must contain variables 'eta_f' and 'h_f' and a 'time' coordinate.
        (This is typically the output of `compute_eta_f`.)
    h_i : array-like | pd.Series | xr.DataArray
        Inner-reef depth time series aligned to ds_eta_f.time.
        If a plain array is given, its length must equal len(ds_eta_f.time).
    force_positive_depth : bool, default True
        If True, take absolute value of h_i and h_f before computing.
    return_new_ds : bool, default True
        If True, return a new Dataset with 'eta_i' added.
        If False, return just the eta_i DataArray.

    Returns
    -------
    xr.Dataset or xr.DataArray
        If return_new_ds=True: dataset with a new variable 'eta_i'.
        Otherwise: the 'eta_i' DataArray.
    """

    # --- checks
    if not isinstance(ds_eta_f, xr.Dataset):
        raise TypeError("ds_eta_f must be an xarray.Dataset.")
    for v in ("eta_f", "h_f"):
        if v not in ds_eta_f:
            raise ValueError(f"ds_eta_f is missing required variable '{v}'.")
    if "time" not in ds_eta_f.coords:
        raise ValueError("ds_eta_f must have a 'time' coordinate.")

    # --- target time from eta_f dataset ---
    target_time = pd.to_datetime(ds_eta_f["time"].values)

    # --- coerce h_i into a Series with its own time index ---
    h_i_arr = np.asarray(h_i, dtype=float)
    ti = pd.to_datetime(h_i_time)
    if getattr(ti, "tz", None) is not None:
        ti = ti.tz_convert("UTC").tz_localize(None)
    h_i_ser = pd.Series(h_i_arr, index=ti)

    # --- Intersect timestamps before aligning ---
    common_times = target_time.intersection(h_i_ser.index)

    # Subselect both to shared timestamps
    h_i_common = h_i_ser.loc[common_times]
    target_common = pd.to_datetime(common_times)

    # Now build DataArray for further xr operations
    h_i_da = xr.DataArray(
        h_i_common.values.astype(float),
        coords={"time": target_common.values.astype("datetime64[ns]")},
        dims=("time",),
        name="h_i"
    )

    # --- align to dataset time using selected common times
    if align == "nearest":
        tol = None if tolerance is None else pd.to_timedelta(tolerance)
        if tol is None:
            h_i_on_face = h_i_da.sel(time=target_common, method="nearest")
        else:
            h_i_on_face = h_i_da.sel(time=target_common, method="nearest", tolerance=tol)
    elif align == "interp":
        h_i_on_face = h_i_da.interp(time=target_common)
    else:
        raise ValueError("align must be 'nearest' or 'interp'.")

    # --- depths from ds (subselect to matching time window) ---
    h_f = ds_eta_f["h_f"].sel(time=target_common)

        # --- enforce positive depth ---
    if force_positive_depth:
        h_i_on_face = np.abs(h_i_on_face)
        h_f = np.abs(h_f)

    h_i_on_face = xr.where(h_i_on_face > 0.0, h_i_on_face, np.nan)
    h_f = xr.where(h_f > 0.0, h_f, np.nan)

    # --- enforce positive depth ---
    if force_positive_depth:
        h_i_on_face = np.abs(h_i_on_face)
        h_f = np.abs(h_f)

    h_i_on_face = xr.where(h_i_on_face > 0.0, h_i_on_face, np.nan)
    h_f = xr.where(h_f > 0.0, h_f, np.nan)

    # --- compute eta_i ---
    eta_i = (ds_eta_f["eta_f"].sel(time=target_common) +
             (h_i_on_face - h_f)).rename("eta_i")

    eta_i.attrs.update({
        "equation": "eta_i = eta_f + (h_i - h_f)",
        "aligned": align,
        "tolerance": tolerance,
    })

    if return_new_ds:
        out = ds_eta_f.sel(time=target_common).copy()
        out["h_i"] = h_i_on_face
        out["eta_i"] = eta_i
        return out
    else:
        return eta_i