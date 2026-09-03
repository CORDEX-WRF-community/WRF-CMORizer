#!/usr/bin/env python
# coding: utf-8
#
# runctrl.vars.nml generator
# The script generates namelist for CORDEX CMIP6 variables based on the list published on github:
# https://raw.githubusercontent.com/WCRP-CORDEX/data-request-table/main/data-request/dreq_default.csv
#
# To run the script a csv file with variables containing all the metadata is necessary to be placed in the running directory.
# The csv file is CORDEX_CMIP6_variables.csv
#
#
# To run the script in the command line:
#     python generate_vars_namelist.py custom var1 var2 var3  --> runctrl.vars.custom.nml will be created, with info on var1 var2 var3
#     python generate_vars_namelist.py var1  --> runctrl.vars.var1.nml will be created, with info on var1
#     python generate_vars_namelist.py core  --> runctrl.vars.core.nml will be created, with info on core variables
#     python generate_vars_namelist.py trier1  --> runctrl.vars.trier1.nml will be created, with info on trier1 variables
#
#
# Contact: milovacj@unican.es


import csv, sys, urllib.request
from pathlib import Path

URL = "https://raw.githubusercontent.com/WCRP-CORDEX/data-request-table/refs/heads/main/cmor-table/datasets.csv"
URL_FALLBACK = "https://raw.githubusercontent.com/impetus4change/T32-CPRCM/refs/heads/main/data-request-fpsurbrcc.csv"
REQUEST = Path("data-request.csv")
REQUEST_FALLBACK = Path("data-request_backup.csv")
WRF = Path("CORDEX_CMIP6_variables.csv")
FILETYPE = "s"  # s=wrfout; p=wrfpress; x=wrfxtrm

FREQ = {
    "fx":  ("timefx",  "cmfx"),
    "1hr": ("time1hr", "cm1hr"),
    "3hr": ("time3hr", "cm3hr"),
    "6hr": ("time6hlsr", "cm6hr"),
    "day": ("timeDay", "cmDay"),
    "mon": ("timeMon", "cmMon"),
    "sea": ("timeSea", "cmSea"),
}

PRIORITIES = {"CORE", "TIER1", "TIER2"}


def read_csv(path, key, required):
    with open(path, encoding="utf-8-sig", newline="") as f:
        r = csv.DictReader(f)
        r.fieldnames = [x.strip() for x in r.fieldnames]
        missing = set(required) - set(r.fieldnames)
        if missing:
            raise RuntimeError(f"{path}: missing {', '.join(sorted(missing))}")
        data = {}
        for row in r:
            row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
            if row.get(key):
                data.setdefault(row[key], []).append(row)
        return data


def levels(var):
    v = var.lower()
    if any(x in v for x in ("tas", "huss", "hurs")):
        return 2, "-999"
    if any(x in v for x in ("uas", "vas", "sfcwind")):
        return 10, "-999"

    n = "".join(c for c in var if c.isdigit())
    if not any(x in var for x in ("ua", "va", "wa", "ta", "zg", "hus")):
        return "-999", "-999"
    return (n or "-999", "-999") if var.endswith("m") else ("-999", n or "-999")


def get_variables(args, request):
    groups = {x.upper() for x in args if x.upper() in PRIORITIES}
    if groups:
        return [
            v for v, rows in request.items()
            if any(r.get("priority", "").upper() in groups for r in rows)
        ]
    return [v for v in dict.fromkeys(args) if v in request]

def create_metadata(var, request, wrf):
    rows, wrf_rows = request.get(var, []), wrf.get(var, [])
    if not rows or not wrf_rows:
        return None

    height, plevel = levels(var)

    result = {
        "var_wrf": f'"{wrf_rows[0]["WRF variable"]}"',
        "var_cmip": f'"{var}"',
        "standard_name": f'"{rows[0]["standard_name"]}"',
        "long_name": f'"{rows[0]["long_name"]}"',
        "units": f'"{rows[0]["units"]}"',
        "cell_methods": f'"{rows[0]["cell_methods"]}"',
        "cell_measures": f'"{rows[0]["cell_measures"]}"',
        "height": height,
        "plevel": plevel,
        "positive": f'"{rows[0]["positive"]}"',
        "filetype": f'"{FILETYPE}"',
        "var_comm": f'"{rows[0]["comment"]}"',
    }

    # Store frequency-specific metadata
    for row in rows:
        freq = row.get("frequency")

        if freq in FREQ:

            freq_name = {
                "fx": "fx",
                "1hr": "1hr",
                "3hr": "3hr",
                "6hr": "6hr",
                "day": "Day",
                "mon": "Mon",
                "sea": "Sea",
            }[freq]

            result[f"time{freq_name}"] = "T"
            result[f"cm{freq_name}"] = f'"{row["cell_methods"]}"'
            result[f"cms{freq_name}"] = f'"{row["cell_measures"]}"'

    # ------------------------------------------------------------
    # Ensure all hourly frequencies are available
    #
    # If one hourly frequency exists, use its metadata for
    # the other hourly frequencies.
    #
    # If none exists, use the default values.
    # ------------------------------------------------------------

    hourly = ["1hr", "3hr", "6hr"]

    # Find an available hourly frequency
    source_freq = next(
        (
            freq for freq in hourly
            if f"time{freq}" in result
        ),
        None
    )

    if source_freq:
        # Copy the available hourly information
        for freq in hourly:
            result[f"time{freq}"] = result[f"time{source_freq}"]
            result[f"cm{freq}"] = result[f"cm{source_freq}"]
            result[f"cms{freq}"] = result[f"cms{source_freq}"]

    else:
        # No hourly frequency available: use defaults
        for freq in hourly:
            result[f"time{freq}"] = "T"
            result[f"cm{freq}"] = '"area: mean time: point"'
            result[f"cms{freq}"] = '"area: areacella"'

    return result


def namelist(variables):

    base_titles = [
        "var_wrf",
        "var_cmip",
        "standard_name",
        "long_name",
        "units",
        "height",
        "plevel",
        "positive"
    ]

    end_titles = [
        "filetype",
        "var_comm"
    ]

    freq_titles = []

    # ------------------------------------------------------------
    # Fixed fields
    # ------------------------------------------------------------
    if any(
        key in v
        for v in variables
        for key in ("timeFx", "cmFx", "cmsFx",
                    "timefx", "cmfx", "cmsfx")
    ):
        freq_titles.extend([
            "timefx",
            "cmfx",
            "cmsfx",
        ])
            
        for v in variables:
            if "timefx" not in v:
                v["timefx"] = "T"

            if "cmfx" not in v:
                v["cmfx"] = '"area: mean"'

            if "cmsfx" not in v:
                v["cmsfx"] = '"area: areacella"'

    # ------------------------------------------------------------
    # Always include all hourly frequencies
    # ------------------------------------------------------------
    for suffix in ("1hr", "3hr", "6hr"):
        freq_titles.extend([
            f"time{suffix}",
            f"cm{suffix}",
            f"cms{suffix}",
        ])

    # ------------------------------------------------------------
    # Daily, monthly and seasonal frequencies
    # ------------------------------------------------------------
    for suffix in ("Day", "Mon", "Sea"):
        for key in (
            f"time{suffix}",
            f"cm{suffix}",
            f"cms{suffix}",
        ):
            if any(key in v for v in variables):
                freq_titles.append(key)

    titles = base_titles + freq_titles + end_titles

    width = max(
        len(str(v.get(t, "")))
        for t in titles
        for v in variables
    )

    return "&vars\n" + "\n".join(
        f"   {t:<17} = "
        + " ".join(
            f"{str(v.get(t, '')):<{width}} ,"
            for v in variables
        )
        for t in titles
    ) + "\n/\n"

# ============================================================
# RUN
# ============================================================

args = sys.argv[1:]
name = "_".join(args)

if not args:
    raise SystemExit("Usage: python generate_namelist.py <variable> [variable ...]")

requests = [(URL, REQUEST)]

# Download request files
if "URL_FALLBACK" in globals():
    requests.append((URL_FALLBACK, REQUEST_FALLBACK))

for url, file in requests:
    if not Path(file).exists():
        urllib.request.urlretrieve(url, file)


fields = {
    "out_name", "frequency", "units", "long_name",
    "standard_name", "cell_methods", "cell_measures",
    "positive", "comment"
}

request = read_csv(REQUEST, "out_name", fields)

if "URL_FALLBACK" in globals():
    request = {
        **read_csv(REQUEST_FALLBACK, "out_name", fields),
        **request,
    }


wrf = read_csv(
    WRF, "output variable name",
    {"output variable name", "WRF variable"}
)

result = [
    v for x in get_variables(args, request)
    if (v := create_metadata(x, request, wrf))
]

if not result:
    raise RuntimeError("No variables could be created.")

output = Path(f"runctrl.vars.{name}.nml")
output.write_text(namelist(result), encoding="utf-8")

print(f"Created: {output}")
