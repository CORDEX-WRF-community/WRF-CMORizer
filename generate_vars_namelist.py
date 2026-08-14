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

URL = "https://raw.githubusercontent.com/WCRP-CORDEX/data-request-table/main/data-request/dreq_default.csv"
REQUEST = Path("data-request.csv")
WRF = Path("CORDEX_CMIP6_variables.csv")
FILETYPE = "s"  # s=wrfout; p=wrfpress; x=wrfxtrm

FREQ = {
    "fx":  ("timefx",  "cmfx"),
    "1hr": ("time1hr", "cm1hr"),
    "3hr": ("time3hr", "cm3hr"),
    "6hr": ("time6hr", "cm6hr"),
    "day": ("timeDay", "cmDay"),
    "mon": ("timeMon", "cmMon"),
    "sea": ("timeSea", "cmSea"),
}

PRIORITIES = {"CORE", "TIER1", "TIER2"}


def read_csv(path, key, required):
    with open(path, encoding="utf-8-sig", newline="") as f:
        r = csv.DictReader(f)
        r.fieldnames = [x.strip() for x in r.fieldnames]
        if missing := set(required) - set(r.fieldnames):
            raise RuntimeError(f"{path}: missing {', '.join(sorted(missing))}")
        data = {}
        for row in r:
            row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
            if row.get(key):
                data.setdefault(row[key], []).append(row)
        return data


def agg(row):
    text = row.get("cell_methods", "").lower()
    if "time:" not in text:
        return ""
    return {
        "mean": "mean", "point": "point", "maximum": "max",
        "minimum": "min", "max": "max", "min": "min", "sum": "sum"
    }.get(text.split("time:")[-1].strip().split()[0], "")


def positive(var, rows):
    s = " ".join(r.get("standard_name", "").lower() for r in rows)
    if "down" in s or "incoming" in s or var in {"rsdt", "hfso"}:
        return "'down'"
    if "upward" in s or "upwelling" in s or "outgoing" in s:
        return "'up'"
    return "'-999'"


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
        "cordexID": "999",
        "var_wrf": f"'{wrf_rows[0]['WRF variable']}'",
        "var_cmip": f"'{var}'",
        "standard_name": f"'{rows[0]['standard_name']}'",
        "long_name": f"'{rows[0]['long_name']}'",
        "units": f"'{rows[0]['units']}'",
        "height": height,
        "plevel": plevel,
        "positive": positive(var, rows),
        "filetype": f"'{FILETYPE}'",
        "comment": f"'{rows[0]['comment']}'",
    }

    a = {
        r["frequency"]: agg(r)
        for r in rows
        if r.get("frequency") in FREQ
    }

    a["fx"] = "point"
    a["1hr"] = a.get("1hr") or a.get("3hr") or a.get("6hr") or "point"
    a["3hr"] = a.get("3hr") or "point"
    a["6hr"] = a.get("6hr") or "point"
    a["day"] = a.get("day") or "mean"
    a["mon"] = a.get("mon") or "mean"
    a["sea"] = a.get("sea") or "mean"

    for f, (t, c) in FREQ.items():
        result[t], result[c] = "T", f"'{a[f]}'"

    return result


def namelist(variables):
    titles = [
        "cordexID", "var_wrf", "var_cmip", "standard_name",
        "long_name", "units", "height", "plevel", "positive",
        "timefx", "cmfx", "time1hr", "cm1hr", "time3hr", "cm3hr",
        "time6hr", "cm6hr", "timeDay", "cmDay", "timeMon", "cmMon",
        "timeSea", "cmSea", "filetype", "comment"
    ]

    width = max(len(str(v.get(t, ""))) for t in titles for v in variables)

    return "&vars\n" + "\n".join(
        f"   {t:<17} = " +
        " ".join(f"{v.get(t, '')!s:<{width}} ," for v in variables)
        for t in titles
    ) + "\n/\n"


# ============================================================
# RUN
# ============================================================

args = sys.argv[1:]
name = "_".join(args)

if not args:
    raise SystemExit("Usage: python generate_namelist.py <variable> [variable ...]")

urllib.request.urlretrieve(URL, REQUEST)

request = read_csv(
    REQUEST, "out_name",
    {"out_name", "frequency", "units", "long_name",
     "standard_name", "cell_methods", "priority", "comment"}
)

wrf = read_csv(
    WRF, "output variable name",
    {"output variable name", "WRF variable"}
)

variables = get_variables(args, request)
result = [v for x in variables if (v := create_metadata(x, request, wrf))]

if not result:
    raise RuntimeError("No variables could be created.")

output = Path(f"runctrl.vars.{name}.nml")
output.write_text(namelist(result), encoding="utf-8")

print(f"Created: {output}")
