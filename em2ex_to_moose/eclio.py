"""
eclio.py - dependency-light readers for Eclipse-format reservoir model files
as written by Petrel / Eclipse / OPM Flow / tNavigator / IX (Eclipse export).

  * read_unformatted(path)  -> list of (KEYWORD, numpy array) from
                               .EGRID / .INIT / .UNRST / .Xnnnn / .GRID (binary, big-endian)
  * read_grdecl(path)       -> dict of keyword -> array from ASCII GRDECL (Petrel "Eclipse grid" export)
  * Deck(path)              -> light parser for the bits of a .DATA deck we need
                               (EQUIL, DENSITY, PVTW, RTEMP, saturation tables, START/DATES/TSTEP,
                               WELSPECS, COMPDAT, WCONHIST/WCONPROD/WCONINJE/WCONINJH,
                               AQUFETP/AQUCT/AQUANCON, unit system)

Only numpy is required.
"""
from __future__ import annotations

import datetime as _dt
import os
import re
import struct
from collections import OrderedDict, defaultdict

import numpy as np

# ----------------------------------------------------------------------------------------------
# Units
# ----------------------------------------------------------------------------------------------
UNITS = {
    # name: (length->m, pressure->Pa, gas surf rate->sm3/d, liquid surf rate->sm3/d, density->kg/m3)
    "METRIC": dict(length=1.0, pressure=1.0e5, gas_rate=1.0, liq_rate=1.0, density=1.0),
    "FIELD": dict(length=0.3048, pressure=6894.757293168, gas_rate=28.316846592,
                  liq_rate=0.158987294928, density=16.01846337396),
    "LAB": dict(length=0.01, pressure=101325.0, gas_rate=1e-6 * 24.0, liq_rate=1e-6 * 24.0,
                density=1000.0),
}
MD_TO_M2 = 9.869233e-16


# ----------------------------------------------------------------------------------------------
# Binary (unformatted) Eclipse files
# ----------------------------------------------------------------------------------------------
_DTYPES = {
    "INTE": (">i4", 4), "REAL": (">f4", 4), "DOUB": (">f8", 8), "LOGI": (">i4", 4),
    "CHAR": ("S8", 8), "MESS": (None, 0),
}


def _dtype_for(typ: str):
    if typ in _DTYPES:
        return _DTYPES[typ]
    if typ.startswith("C0"):  # C0nn : strings of length nn
        n = int(typ[2:])
        return ("S%d" % n, n)
    raise ValueError("Unknown Eclipse data type %r" % typ)


def _read_record(f):
    head = f.read(4)
    if len(head) < 4:
        return None
    (n,) = struct.unpack(">i", head)
    data = f.read(n)
    tail = f.read(4)
    if len(data) != n or len(tail) != 4 or struct.unpack(">i", tail)[0] != n:
        raise IOError("Corrupt Fortran record in %s" % getattr(f, "name", "?"))
    return data


def read_unformatted(path, stop_after=None):
    """Read an Eclipse unformatted (binary) file. Returns a list of (kw, ndarray).

    stop_after: optional keyword; reading stops after the first occurrence (e.g. 'ENDGRID').
    """
    out = []
    with open(path, "rb") as f:
        while True:
            head = _read_record(f)
            if head is None:
                break
            if len(head) != 16:
                raise IOError("%s: unexpected header length %d (formatted file?)" % (path, len(head)))
            kw = head[:8].decode("ascii", "replace").strip()
            (n,) = struct.unpack(">i", head[8:12])
            typ = head[12:16].decode("ascii", "replace")
            dt, size = _dtype_for(typ)
            chunks, got = [], 0
            while got < n and size > 0:
                rec = _read_record(f)
                if rec is None:
                    raise IOError("%s: truncated data for %s" % (path, kw))
                chunks.append(rec)
                got += len(rec) // size
            if size == 0:
                arr = np.zeros(0)
            else:
                arr = np.frombuffer(b"".join(chunks), dtype=dt, count=n)
                if typ in ("INTE", "LOGI"):
                    arr = arr.astype(np.int32)
                elif typ in ("REAL", "DOUB"):
                    arr = arr.astype(np.float64)
                else:
                    arr = np.array([s.decode("ascii", "replace").strip() for s in arr])
            out.append((kw, arr))
            if stop_after and kw == stop_after:
                break
    return out


def write_unformatted(path, items):
    """Writer used by the test-suite (and handy for debugging). items: list of (kw, array, type)."""
    with open(path, "wb") as f:
        def rec(b):
            f.write(struct.pack(">i", len(b)))
            f.write(b)
            f.write(struct.pack(">i", len(b)))
        for kw, arr, typ in items:
            arr = np.asarray(arr)
            n = len(arr)
            rec(kw.ljust(8).encode() + struct.pack(">i", n) + typ.encode())
            if typ == "CHAR":
                block = 105
                for s in range(0, n, block):
                    rec(b"".join(str(x).ljust(8)[:8].encode() for x in arr[s:s + block]))
            elif typ != "MESS":
                dt = {"INTE": ">i4", "REAL": ">f4", "DOUB": ">f8", "LOGI": ">i4"}[typ]
                for s in range(0, n, 1000):
                    rec(arr[s:s + 1000].astype(dt).tobytes())


def first(items, kw, default=None):
    for k, a in items:
        if k == kw:
            return a
    return default


# ----------------------------------------------------------------------------------------------
# GRDECL (ASCII) reader
# ----------------------------------------------------------------------------------------------
_KW_RE = re.compile(r"^[A-Z][A-Z0-9_+\-]{0,7}$")
_NO_DATA_KWS = {"ECHO", "NOECHO", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY",
                "SCHEDULE", "RUNSPEC", "END", "METRIC", "FIELD", "LAB", "OIL", "WATER", "GAS",
                "DISGAS", "VAPOIL", "NONNC", "NEWTRAN", "OLDTRAN", "UNIFOUT", "UNIFIN",
                "FMTOUT", "FMTIN", "INIT", "GRIDFILE", "NOSIM", "SKIPREST", "MULTOUT", "MULTIN"}


def _strip_comment(line: str) -> str:
    # '--' starts a comment, unless inside quotes
    out, q = [], None
    i = 0
    while i < len(line):
        c = line[i]
        if q:
            out.append(c)
            if c == q:
                q = None
        elif c in "'\"":
            q = c
            out.append(c)
        elif line.startswith("--", i):
            break
        else:
            out.append(c)
        i += 1
    return "".join(out)


def _expand_tokens(tokens, default=np.nan):
    vals = []
    for t in tokens:
        if "*" in t and not t.startswith(("'", '"')):
            n, _, v = t.partition("*")
            n = int(n)
            vals.extend([default if v == "" else v] * n)
        else:
            vals.append(t)
    return vals


def _resolve_path(base_dir, name):
    name = name.strip().strip("'\"").replace("\\", "/")
    cand = name if os.path.isabs(name) else os.path.join(base_dir, name)
    if os.path.exists(cand):
        return cand
    # case-insensitive fallback (files copied off Windows)
    d, b = os.path.split(cand)
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.lower() == b.lower():
                return os.path.join(d, f)
    # last resort: same basename next to the deck
    for f in os.listdir(base_dir):
        if f.lower() == os.path.basename(cand).lower():
            return os.path.join(base_dir, f)
    raise FileNotFoundError("INCLUDE file not found: %s" % name)


_SECTIONS = {"RUNSPEC", "GRID", "EDIT", "PROPS", "REGIONS", "SOLUTION", "SUMMARY", "SCHEDULE"}


def _logical_lines(path, _depth=0, skip_sections=(), _state=None, warnings=None):
    """Yield comment-stripped lines, following INCLUDE statements.

    skip_sections: INCLUDEs met while inside these deck sections are not opened (e.g. the huge
    GRID-section files when only the schedule is needed). Missing includes raise unless a
    `warnings` list is given, in which case they are noted there and skipped.
    """
    if _state is None:
        _state = {"section": None}
    base = os.path.dirname(os.path.abspath(path))
    with open(path, "r", errors="replace") as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        s = _strip_comment(lines[i]).strip()
        i += 1
        if s.upper() in _SECTIONS:
            _state["section"] = s.upper()
        if s.upper() == "INCLUDE":
            # file name may be on the next non-empty line(s), terminated by '/'
            buf = ""
            while "/" not in buf.replace("\\", "") or buf.count("'") % 2 == 1:
                if i >= len(lines):
                    break
                nxt = _strip_comment(lines[i]).strip()
                i += 1
                buf += " " + nxt
                if re.search(r"['\"][^'\"]*['\"]\s*/", buf) or re.search(r"^\s*\S+\s+/", buf):
                    break
            m = re.search(r"['\"]([^'\"]+)['\"]", buf)
            fname = m.group(1) if m else buf.split("/")[0].strip()
            if _state["section"] in skip_sections:
                continue
            try:
                inc = _resolve_path(base, fname)
            except FileNotFoundError as err:
                if warnings is None:
                    raise
                warnings.append(str(err))
                continue
            if _depth > 20:
                raise RecursionError("INCLUDE nesting too deep")
            yield from _logical_lines(inc, _depth + 1, skip_sections, _state, warnings)
            continue
        yield s


def _tokenize(line):
    return re.findall(r"'[^']*'|\"[^\"]*\"|/|[^\s/]+", line)


def read_grdecl(path):
    """Parse a GRDECL file (and its INCLUDEs). Returns OrderedDict kw -> np.ndarray (float or str)."""
    out = OrderedDict()
    kw, buf = None, []
    for line in _logical_lines(path):
        if not line:
            continue
        toks = _tokenize(line)
        if kw is None:
            if len(toks) == 1 and _KW_RE.match(toks[0].upper()):
                kw, buf = toks[0].upper(), []
                if kw in _NO_DATA_KWS:
                    kw = None
            continue
        # inside keyword data
        if len(toks) == 1 and _KW_RE.match(toks[0]) and not re.match(r"^[0-9]", toks[0]) \
                and toks[0] not in ("F", "T"):
            # a new keyword started before '/', previous one had no data
            kw, buf = toks[0].upper(), []
            if kw in _NO_DATA_KWS:
                kw = None
            continue
        if "/" in toks:
            buf.extend(toks[:toks.index("/")])
            vals = _expand_tokens(buf)
            try:
                arr = np.array([float(v) for v in vals])
            except (ValueError, TypeError):
                arr = np.array([str(v).strip("'\"") for v in vals])
            out[kw] = arr
            kw, buf = None, []
        else:
            buf.extend(toks)
    return out


# ----------------------------------------------------------------------------------------------
# .DATA deck parser (only what the converter needs)
# ----------------------------------------------------------------------------------------------
_MULTI_RECORD = {"COMPDAT", "WELSPECS", "WCONHIST", "WCONPROD", "WCONINJE", "WCONINJH", "DATES",
                 "AQUFETP", "AQUCT", "AQUANCON", "COMPDATL", "WELOPEN", "MULTREGT"}
_ONE_RECORD = {"START", "RTEMP", "TSTEP", "TABDIMS", "EQLDIMS", "AQUDIMS", "DIMENS", "WELLDIMS",
               "RTEMPA", "PETGRID", "GRIDFILE", "GDFILE", "SATOPTS", "EHYSTR"}
_TABLES = {"SWFN", "SGFN", "SGWFN", "SWOF", "SGOF", "PVDG", "PVTW", "DENSITY", "EQUIL", "ROCK"}
_MONTHS = {m: i + 1 for i, m in enumerate(
    ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JLY", "AUG", "SEP", "OCT", "NOV", "DEC"])}
_MONTHS["JUL"] = 7


def _parse_value(t):
    if t is None:
        return None
    if isinstance(t, float):
        return t
    s = str(t)
    if s.startswith(("'", '"')):
        return s.strip("'\"")
    try:
        return float(s)
    except ValueError:
        return s


class Deck:
    """Very small Eclipse deck reader. Keeps raw records for the keywords it understands."""

    def __init__(self, path):
        self.path = path
        self.kw = defaultdict(list)          # kw -> list of (records, sched_index)
        self.sequence = []                   # schedule-ordered (kw, records)
        self.units = "METRIC"
        self.phases = set()
        self.warnings = []
        self._parse()

    # ---- low level
    def _parse(self):
        lines = [ln for ln in _logical_lines(self.path, skip_sections=("GRID", "EDIT", "REGIONS"),
                                             warnings=self.warnings)]
        dims = dict(TABDIMS=[1, 1], EQLDIMS=[1])
        i = 0
        n = len(lines)

        def read_records(i, nrec=None, tables=False):
            """Read records ('/'-terminated). nrec=None: until empty record."""
            recs, cur, ntab = [], [], 0
            while i < n:
                ln = lines[i]
                if not ln:
                    i += 1
                    continue
                toks = _tokenize(ln)
                # keyword line appears => stop (defensive)
                if nrec is None and not cur and len(toks) == 1 and _KW_RE.match(toks[0]) \
                        and not toks[0].startswith("'") and toks[0] not in ("/",):
                    break
                i += 1
                for t in toks:
                    if t == "/":
                        recs.append(cur)
                        cur = []
                        if nrec is None and not recs[-1]:
                            recs.pop()
                            return recs, i
                        if nrec is not None and len(recs) >= nrec:
                            return recs, i
                        break  # anything after '/' on a line is a comment
                    cur.append(t)
            if cur:
                recs.append(cur)
            return recs, i

        while i < n:
            ln = lines[i]
            toks = _tokenize(ln) if ln else []
            if len(toks) == 1 and _KW_RE.match(toks[0].upper()):
                kw = toks[0].upper()
                i += 1
                if kw in ("METRIC", "FIELD", "LAB"):
                    self.units = kw
                    continue
                if kw in ("OIL", "WATER", "GAS"):
                    self.phases.add(kw)
                    continue
                if kw in _MULTI_RECORD:
                    recs, i = read_records(i)
                elif kw in _ONE_RECORD:
                    recs, i = read_records(i, 1)
                elif kw in _TABLES:
                    if kw in ("PVTW", "DENSITY", "PVDG", "ROCK"):
                        ntab = dims["TABDIMS"][1]
                    elif kw == "EQUIL":
                        ntab = dims["EQLDIMS"][0]
                    else:
                        ntab = dims["TABDIMS"][0]
                    recs, i = read_records(i, ntab, tables=True)
                else:
                    continue  # unknown keyword: its data lines are skipped by the outer loop
                recs = [[_parse_value(v) for v in _expand_tokens(r, default=None)] for r in recs]
                if kw in ("TABDIMS", "EQLDIMS") and recs:
                    for k, v in enumerate(recs[0][:2]):
                        if v is not None and k < len(dims[kw]):
                            dims[kw][k] = int(v)
                self.kw[kw].append(recs)
                self.sequence.append((kw, recs))
            else:
                i += 1

    # ---- helpers
    def get(self, kw):
        return self.kw.get(kw, [None])[0]

    @property
    def u(self):
        return UNITS[self.units]

    def start_date(self):
        r = self.get("START")
        if not r:
            return None
        d, m, y = r[0][0], r[0][1], r[0][2]
        return _dt.datetime(int(y), _MONTHS[str(m).upper()[:3]], int(d))

    def schedule(self):
        """Walk the schedule. Returns dict with per-well time series and connections.

        times are seconds since START. Rates: mass is computed later (needs densities).
        """
        start = self.start_date() or _dt.datetime(2000, 1, 1)
        t = 0.0
        times = [0.0]
        welspecs, compdat = {}, defaultdict(list)
        events = defaultdict(list)  # well -> list of (t, kind, gas_sm3d, wat_sm3d, status)
        u = self.u
        for kw, recs in self.sequence:
            if kw == "DATES":
                for r in recs:
                    d, m, y = r[0], r[1], r[2]
                    hh = 0
                    if len(r) > 3 and r[3]:
                        try:
                            hh = sum(float(x) * f for x, f in zip(str(r[3]).split(":"), (3600, 60, 1)))
                        except ValueError:
                            hh = 0
                    dt = _dt.datetime(int(y), _MONTHS[str(m).upper()[:3]], int(d)) + _dt.timedelta(seconds=hh)
                    t = (dt - start).total_seconds()
                    times.append(t)
            elif kw == "TSTEP":
                for v in recs[0]:
                    t += float(v) * 86400.0
                    times.append(t)
            elif kw == "WELSPECS":
                for r in recs:
                    r = r + [None] * 6
                    welspecs[r[0]] = dict(group=r[1], i=r[2], j=r[3], ref_depth=r[4], phase=r[5])
            elif kw == "COMPDAT":
                for r in recs:
                    r = r + [None] * 14
                    name = r[0]
                    if "*" in str(name) or "?" in str(name):
                        continue
                    ws = welspecs.get(name, {})
                    i_ = r[1] if r[1] not in (None, 0) else ws.get("i")
                    j_ = r[2] if r[2] not in (None, 0) else ws.get("j")
                    status = (r[5] or "OPEN").upper()
                    k1, k2 = int(r[3]), int(r[4])
                    for k in range(k1, k2 + 1):
                        compdat[name] = [c for c in compdat[name] if (c["i"], c["j"], c["k"]) != (int(i_), int(j_), k)]
                        compdat[name].append(dict(i=int(i_), j=int(j_), k=k, status=status,
                                                  cf=r[7], diam=r[8], kh=r[9], skin=r[10]))
            elif kw in ("WCONHIST", "WCONPROD"):
                for r in recs:
                    r = r + [None] * 12
                    name, status = r[0], str(r[1] or "OPEN").upper()
                    o, w, g = [float(x) if isinstance(x, float) else 0.0 for x in r[3:6]]
                    if kw == "WCONPROD":
                        cmode = str(r[2] or "").upper()
                        # prediction mode: only explicit rate targets can be represented
                        if cmode == "GRAT":
                            w = 0.0
                        elif cmode == "WRAT":
                            g = 0.0
                        else:
                            self.warnings.append("WCONPROD %s control %s at t=%.0f d is not a rate target; "
                                                 "well set to zero rate" % (name, cmode, t / 86400))
                            g = w = 0.0
                    if status != "OPEN":
                        g = w = 0.0
                    events[name].append((t, "PROD", g * u["gas_rate"], w * u["liq_rate"], status, kw))
            elif kw in ("WCONINJE", "WCONINJH"):
                for r in recs:
                    r = r + [None] * 10
                    name, typ, status = r[0], str(r[1] or "").upper(), str(r[2] or "OPEN").upper()
                    if kw == "WCONINJE":
                        rate = r[4] if isinstance(r[4], float) else 0.0
                    else:
                        rate = r[3] if isinstance(r[3], float) else 0.0
                    if status != "OPEN":
                        rate = 0.0
                    conv = u["gas_rate"] if typ.startswith("GA") else u["liq_rate"]
                    g = rate * conv if typ.startswith("GA") else 0.0
                    w = rate * conv if typ.startswith("WA") else 0.0
                    events[name].append((t, "INJ", g, w, status, kw))
        return dict(start=start, times=times, end=max(times), welspecs=welspecs,
                    compdat=dict(compdat), events=dict(events))

    # ---- physical data
    def densities(self):
        """(oil, water, gas) surface densities in kg/m3."""
        r = self.get("DENSITY")
        d = self.u["density"]
        if r and r[0]:
            vals = (r[0] + [None] * 3)[:3]
            defaults = (859.0, 1000.0, 0.8)
            return tuple(float(v) * d if isinstance(v, float) else dflt for v, dflt in zip(vals, defaults))
        return (859.0, 1000.0, 0.8)

    def pvtw(self):
        """PVTW table 1 -> dict(pref [Pa], bw [-], cw [1/Pa], mu [Pa s], cv [1/Pa])."""
        r = self.get("PVTW")
        if not r or not r[0]:
            return None
        x = (r[0] + [None] * 5)[:5]
        if not isinstance(x[0], float):
            return None
        P = self.u["pressure"]
        return dict(pref=x[0] * P, bw=x[1] if isinstance(x[1], float) else 1.0,
                    cw=(x[2] if isinstance(x[2], float) else 4e-5 / 1e5 * P) / P,
                    mu=(x[3] if isinstance(x[3], float) else 0.5) * 1e-3,
                    cv=(x[4] if isinstance(x[4], float) else 0.0) / P)

    def pvdg(self):
        """PVDG table 1 -> (p [Pa], Bg [rm3/sm3], mu [Pa s]) arrays, or None."""
        r = self.get("PVDG")
        if not r or not r[0]:
            return None
        v = np.array([x for x in r[0] if isinstance(x, float)])
        v = v[: len(v) // 3 * 3].reshape(-1, 3)
        bg_scale = self.u["liq_rate"] / self.u["gas_rate"]   # rb/Mscf -> rm3/sm3 (FIELD); 1 for METRIC
        return v[:, 0] * self.u["pressure"], v[:, 1] * bg_scale, v[:, 2] * 1e-3

    def equil(self):
        r = self.get("EQUIL")
        if not r:
            return None
        x = (r[0] + [None] * 9)
        L, P = self.u["length"], self.u["pressure"]
        return dict(datum_depth=float(x[0]) * L, datum_p=float(x[1]) * P,
                    contact1=float(x[2]) * L if isinstance(x[2], float) else None,
                    contact2=float(x[4]) * L if isinstance(x[4], float) else None)

    def temperature(self):
        for k in ("RTEMP", "RTEMPA"):
            r = self.get(k)
            if r and r[0] and isinstance(r[0][0], float):
                v = r[0][0]
                return (v - 32.0) * 5 / 9 + 273.15 if self.units == "FIELD" else v + 273.15
        return None

    def satfunc(self, table=1):
        """Gas-water curves from saturation table number `table` (1-based), if present.

        keys: sw, krw, sg, krg (numpy arrays), src
        """
        def tab(kw, ncol):
            recs = self.get(kw)
            if not recs or len(recs) < table or not recs[table - 1]:
                return None
            v = [x for x in recs[table - 1] if isinstance(x, float)]
            return np.array(v[: len(v) // ncol * ncol]).reshape(-1, ncol)
        out = {}
        P = self.u["pressure"]
        t = tab("SGWFN", 4)
        if t is not None:
            # SGWFN: Sg, krg, krw, Pcgw (Pcgw tabulated against Sg; stored here against Sw)
            out.update(sg=t[:, 0], krg=t[:, 1], sw=1 - t[:, 0], krw=t[:, 2], pc=t[:, 3] * P, src="SGWFN")
        a, b = tab("SWFN", 3), tab("SGFN", 3)
        if a is not None:
            out.update(sw=a[:, 0], krw=a[:, 1], pc=a[:, 2] * P, src="SWFN/SGFN")
        if b is not None:
            out.update(sg=b[:, 0], krg=b[:, 1])
        if not out:
            a, b = tab("SWOF", 4), tab("SGOF", 4)
            if a is not None:
                out.update(sw=a[:, 0], krw=a[:, 1], pc=a[:, 3] * P, src="SWOF/SGOF")
            if b is not None:
                out.update(sg=b[:, 0], krg=b[:, 1])
        if out and "src" in out:
            out["src"] += " table %d" % table
        return out or None

    def hysteresis(self):
        """Saturation-function hysteresis (SATOPTS HYSTER + EHYSTR), or None if it is off.

        kr_model: EHYSTR item 2 (0/1 Carlson, 2/3 Killough for the non-wetting phase; 1/3: imbibition
        curve for the wetting phase; 4+: Killough for both / other options), applies_to: item 5
        (KR, PC or BOTH), pc_curvature: item 1.
        """
        so = self.get("SATOPTS")
        if not so or not any(isinstance(v, str) and v.upper().startswith("HYSTER") for v in so[0]):
            return None
        r = (self.get("EHYSTR") or [[]])[0]

        def item(k, d):
            return r[k] if len(r) > k and r[k] is not None else d
        return dict(kr_model=int(item(1, 0)), applies_to=str(item(4, "BOTH")).upper(),
                    pc_curvature=float(item(0, 0.1)))

    def n_satfunc(self):
        for kw in ("SGWFN", "SWFN", "SWOF"):
            r = self.get(kw)
            if r:
                return len(r)
        return 0

    def equil_all(self):
        """All EQUIL records (one per EQLNUM region)."""
        r = self.get("EQUIL") or []
        L, P = self.u["length"], self.u["pressure"]
        out = []
        for rec in r:
            x = rec + [None] * 9
            if not isinstance(x[0], float):
                continue
            out.append(dict(datum_depth=float(x[0]) * L, datum_p=float(x[1]) * P,
                            contact1=float(x[2]) * L if isinstance(x[2], float) else None))
        return out

    def multregt(self):
        """MULTREGT records -> list of (r1, r2, mult, region_array_letter)."""
        out = []
        for recs in self.kw.get("MULTREGT", []):
            for r in recs:
                r = r + [None] * 6
                if not isinstance(r[0], float) or not isinstance(r[1], float) or not isinstance(r[2], float):
                    continue
                out.append((int(r[0]), int(r[1]), float(r[2]), str(r[5] or "M").upper()[:1]))
        return out

    def aquifers(self):
        """Analytic aquifers: dict id -> dict(type, props..., connections=[...])"""
        L, P = self.u["length"], self.u["pressure"]
        aq = {}
        for recs in self.kw.get("AQUFETP", []):
            for r in recs:
                r = r + [None] * 10
                aq[int(r[0])] = dict(type="FETKOVICH", datum_depth=_f(r[1], L), p0=_f(r[2], P),
                                     volume=_f(r[3], self.u["liq_rate"]), ct=_f(r[4], 1 / P), pi_raw=r[5],
                                     connections=[])
        for recs in self.kw.get("AQUCT", []):
            for r in recs:
                r = r + [None] * 14
                aq[int(r[0])] = dict(type="CARTER_TRACY", datum_depth=_f(r[1], L), p0=_f(r[2], P),
                                     perm_md=r[3], poro=r[4], ct=_f(r[5], 1 / P), r0=_f(r[6], L),
                                     h=_f(r[7], L), theta=r[8], connections=[])
        for recs in self.kw.get("AQUANCON", []):
            for r in recs:
                r = r + [None] * 11
                aid = int(r[0])
                aq.setdefault(aid, dict(type="UNKNOWN", datum_depth=None, p0=None, connections=[]))
                aq[aid]["connections"].append(dict(i1=int(r[1]), i2=int(r[2]), j1=int(r[3]), j2=int(r[4]),
                                                   k1=int(r[5]), k2=int(r[6]), face=str(r[7]).upper()))
        return aq


def _f(v, scale):
    return float(v) * scale if isinstance(v, float) else None
