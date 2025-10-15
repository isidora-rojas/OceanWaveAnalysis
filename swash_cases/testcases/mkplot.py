# mkplot_py.py
import os, glob, time
import numpy as np
import matplotlib.pyplot as plt

def try_netcdf():
    nc_files = sorted(glob.glob("*.nc"))
    if not nc_files:
        return False
    try:
        import netCDF4 as nc
        ds = nc.Dataset(nc_files[0])
        # Try to guess common var names
        keys = {k.lower(): k for k in ds.variables.keys()}
        # x
        xvar = None
        for k in ("x", "x-coordinate", "xcoord"):
            if k in keys: xvar = keys[k]; break
        # time
        tvar = None
        for k in ("time", "t"):
            if k in keys: tvar = keys[k]; break
        # eta (surface)
        etavar = None
        for k in ds.variables.keys():
            kl = k.lower()
            if "eta" in kl or "zs" in kl or "zeta" in kl:
                etavar = k; break
        if xvar and tvar and etavar:
            x = ds.variables[xvar][:]
            t = ds.variables[tvar][:]
            eta = ds.variables[etavar][:]   # expect (time, x) or (time, y, x)
            print(f"[NetCDF] {nc_files[0]}: eta shape {eta.shape}")
            if eta.ndim == 2 and eta.shape[1] == x.shape[0]:
                plt.figure()
                plt.imshow(eta.T, aspect='auto', origin='lower',
                           extent=[t.min(), t.max(), x.min(), x.max()])
                plt.xlabel("Time (s)"); plt.ylabel("x (m)")
                plt.title("Surface elevation η(t, x)")
                plt.colorbar(label="η (m)")
                plt.show()
            else:
                print("NetCDF present but shape not (time, x). Skipping image.")
        else:
            print("NetCDF found, but expected vars not present:", list(ds.variables.keys()))
        ds.close()
        return True
    except Exception as e:
        print("NetCDF read failed:", e)
        return False

def load_text(fname):
    """Load whitespace-delimited numeric data, skipping commented lines."""
    try:
        return np.loadtxt(fname, comments=('#','%','!','/','*'))
    except Exception:
        # try genfromtxt with auto-skip
        return np.genfromtxt(fname, comments='#')

def recent_numeric_files(hours=6):
    exts = ("*.dat", "*.txt", "*.prt", "*.erf")
    cand = []
    now = time.time()
    for pat in exts:
        for f in glob.glob(pat):
            try:
                if now - os.path.getmtime(f) <= hours*3600:
                    cand.append(f)
            except Exception:
                pass
    return sorted(cand, key=lambda f: os.path.getmtime(f))

def try_text():
    files = recent_numeric_files()
    if not files:
        # fallback: all numeric-looking files
        files = sorted(glob.glob("*.dat") + glob.glob("*.txt") + glob.glob("*.prt") + glob.glob("*.erf"))
    plotted_any = False
    for f in files:
        try:
            arr = load_text(f)
            if not isinstance(arr, np.ndarray) or arr.size == 0:
                continue
            print(f"[Text] {f} shape {arr.shape}")
            # 1D vector?
            if arr.ndim == 1 and arr.size > 1:
                plt.figure()
                plt.plot(arr)
                plt.title(f"{f} (index vs value)")
                plt.xlabel("Index"); plt.ylabel("Value")
                plt.show()
                plotted_any = True
                continue
            # 2D matrix
            if arr.ndim == 2:
                nrow, ncol = arr.shape
                if ncol == 2:
                    # assume col0 = time or x, col1 = value
                    plt.figure()
                    plt.plot(arr[:,0], arr[:,1])
                    plt.title(f"{f} (col2 vs col1)")
                    plt.xlabel("Col 1"); plt.ylabel("Col 2")
                    plt.show()
                    plotted_any = True
                elif ncol > 2:
                    # try: first column is time, remaining are probes/space
                    t = arr[:,0]
                    V = arr[:,1:]
                    # If looks like time vs x snapshots, show image
                    plt.figure()
                    plt.imshow(V.T, aspect='auto', origin='lower',
                               extent=[t.min(), t.max(), 0, V.shape[1]])
                    plt.xlabel("Time (s)"); plt.ylabel("Index (probe/space)")
                    plt.title(f"{f} (matrix: values vs time)")
                    plt.colorbar(label="Value")
                    plt.show()
                    # also plot up to 3 columns as time series
                    for j in range(min(3, V.shape[1])):
                        plt.figure()
                        plt.plot(t, V[:,j])
                        plt.xlabel("Time (s)"); plt.ylabel(f"Series {j+1}")
                        plt.title(f"{f} — series {j+1}")
                        plt.show()
                    plotted_any = True
        except Exception as e:
            print(f"Skipping {f} due to read error:", e)
    return plotted_any

if __name__ == "__main__":
    ok = try_netcdf()
    if not ok:
        ok = try_text()
    if not ok:
        print("No recognizable outputs to plot. Files present:")
        os.system("ls -lt")
