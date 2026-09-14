r"""ISMIP7 output, part one: what the forward records every year.

The data request (``isschecker/data/ISMIP7_variable_request.csv`` in
``ismip/ISM_SimulationChecker``, bundled here as
``icepack2_tools/ismip7_variable_request.csv``) asks for yearly 2D fields on
the 8 km AIS grid and yearly scalars. Regridding is a serial post-processing step
(``antarctica/scripts/write_ismip7_output.py``); this module is the part that
runs inside the parallel forward, gated by ``ISMIP7_OUTPUT=1``:

* state variables (``ST``) are snapshots at the end of each year, stamped
  1 January of the following year by the writer;
* flux variables (``FL``) are the year's means, accumulated every transport
  advance from the sources REQUESTED of the transport, stamped 1 July. That
  is the forcing SMB as handed to the transport, BEFORE the positivity
  limiter clips a net sink that would draw a cell below ``h_clamp``, plus the
  ocean melt. The withheld part is the run's ``clamp`` budget column and is
  not any ISMIP7 variable, so in the thin front cells where the limiter fires
  the grid budget does not close against ``dlithkdt``;
* the scalars are the integrals of the same fields, written to a CSV as the
  run goes so an early stop loses nothing.

Everything is kept on the model's own mesh in a Firedrake checkpoint,
``<results>/<experiment>_<lc>_ismip7_annual.h5``, one entry per variable and
year (0-based ``idx``, with the ``years`` attribute mapping position to
year), in the model's units (m, m/yr, MPa); the writer converts to the
request's SI units and applies the fill policies. A chained run appends to
the file it finds, so the six links of a 285-year projection leave one
continuous series.

Conventions (from the request and discussions #16, #19, #22):

* ``acabf`` is the forcing surface mass balance REQUESTED of the transport
  (RACMO climatology plus the re-referenced anomaly), not what survived the
  positivity limiter. The apparent-mass-balance reference ``a_ref`` is NOT
  part of it: it cancels the discrete flux
  divergence spike by spike (up to ~1000 m/yr at the Pine Island grounding
  zone) and reported as SMB it would sit two orders of magnitude outside the
  request's range. It is recorded separately as ``acabf_correction`` (m/yr
  ice, not a request variable) so the grid budget can be closed by anyone
  who needs it, and the README states the convention.
* ``libmassbffl`` is the ocean melt on floating cells (negative = loss);
  ``libmassbfgr`` is zero (no grounded basal melt in the model);
  ``lifmassbf`` is zero (no frontal melt distinct from the basal melt).
* ``licalvf`` is the ice removed at the front, negative = loss, booked in the
  cell it was removed from (whole-cell removal, sub-cell shed, retreat
  slivers), the same tallies as the ``calv`` budget column.
* ``ligroundf`` is the flux across the grounding line, booked as a specific
  mass flux into the first FLOATING cell (discussion #22), positive for ice
  leaving grounded ice; the grid sum of ``ligroundf * area`` is the
  grounding-line discharge.
* ``dlithkdt`` is the change in thickness over the year divided by the year.
* the three area fractions are cell indicators here (0 or 1 per DG0 cell);
  the conservative regridding to 8 km turns them into fractions.
"""
import os

import numpy as np
from firedrake import Function, TestFunction, assemble, dS, dx
import firedrake as fd

SECONDS_PER_YEAR = 31556926.0
RHO_I = 917.0

#: request name -> (type, how this model produces it)
VARIABLES_2D = {
    "lithk": "ST", "orog": "ST", "topg": "ST", "base": "ST",
    "sftgif": "ST", "sftgrf": "ST", "sftflf": "ST",
    "xvelmean": "ST", "yvelmean": "ST", "xvelsurf": "ST", "yvelsurf": "ST",
    "xvelbase": "ST", "yvelbase": "ST", "strbasemag": "ST",
    "acabf": "FL", "libmassbfgr": "FL", "libmassbffl": "FL", "dlithkdt": "FL",
    "licalvf": "FL", "ligroundf": "FL", "lifmassbf": "FL",
}
SCALARS = ("lim", "limnsw", "iareagr", "iareafl", "tendacabf", "tendlibmassbfgr",
           "tendlibmassbffl", "tendlicalvf", "tendlifmassbf", "tendligroundf")


class AnnualOutput:
    r"""Accumulates a year of a forward run and writes it at the year end.

    ``begin_step`` / ``commit_step`` bracket one time step so a rewound
    (subcycled) attempt does not double-count: ``_advance`` books into the
    step tallies, and only a completed step is added to the year.

    A chained run resumes into the same files: when ``out_path`` already
    exists it is opened in append mode and its ``years`` attribute seeds the
    written years. The partly accumulated year survives the link boundary
    too: ``state_fields`` / ``state_attrs`` hand it to the run's own
    checkpoint and ``resume`` takes it back, so a job that stops at 2021.4
    goes on accumulating 2021 rather than losing four months of flux or
    relabelling them. The resumed year is taken from that state and must be
    exactly one past the last year written, so no year is ever renamed.
    """

    #: the per-cell year sums carried across a chained resume
    ACCUMULATORS = ("acabf", "acabf_correction", "libmassbffl", "licalvf",
                    "ligroundf")
    #: dataset names of the in-progress year inside the run's checkpoint
    STATE_PREFIX = "ismip7_acc_"
    STATE_THICKNESS = "ismip7_year_start_thickness"
    STATE_YEAR = "ismip7_year"
    STATE_YEAR_TIME = "ismip7_year_time"

    def __init__(self, mesh, Q_dg, V, out_path, scalars_path, first_year, rho_ratio,
                 comm=None, log=None, resume=None):
        self.mesh, self.Q_dg, self.V = mesh, Q_dg, V
        self.out_path, self.scalars_path = out_path, scalars_path
        self.year = int(np.floor(float(first_year) + 1e-6))   # the year being accumulated (its Jan 1 has passed)
        self.rho_ratio = float(rho_ratio)
        self.comm = comm or mesh.comm
        self.log = log or (lambda s: None)
        # owned cells only: dat.data_ro is the owned slice, dof_count counts
        # the halo too (5650 vs 3772 on one of two ranks of the 32 km mesh)
        self.cell_area = assemble(TestFunction(Q_dg) * dx).dat.data_ro.copy()
        n = len(self.cell_area)
        self.year_acc = {k: np.zeros(n) for k in self.ACCUMULATORS}
        self.step_acc = {k: np.zeros(n) for k in self.year_acc}
        self.year_time = 0.0
        self.h_year_start = None
        if resume is not None:
            if int(resume["year"]) != self.year:
                raise ValueError(
                    f"the restart checkpoint was accumulating ISMIP7 year "
                    f"{int(resume['year'])} but its timeline resumes at "
                    f"t={float(first_year)!r}, which falls in year {self.year}; "
                    f"the checkpoint's ISMIP7 state does not belong to it."
                )
            for k in self.ACCUMULATORS:
                self.year_acc[k][:] = resume["acc"][k]
            self.h_year_start = np.asarray(resume["h_year_start"]).copy()
            self.year_time = float(resume["year_time"])
        self._phi = TestFunction(Q_dg)
        self._n = fd.FacetNormal(mesh)
        self._gl_cof = fd.Cofunction(Q_dg.dual())
        # A chained job re-enters run_simulation and rebuilds this object, so
        # an existing annual file is appended to: opening it "w" would
        # truncate every year the earlier links wrote while the scalars CSV
        # (opened "a") kept them, and the writer reads its year list here.
        self._written_years = []
        self._chk_mode = "w"
        if os.path.exists(out_path):
            with fd.CheckpointFile(out_path, "a") as chk:
                years_attr = chk.get_attr("/", "years") if chk.has_attr("/", "years") else ""
            self._written_years = [int(y) for y in str(years_attr).split(",") if y]
            self._chk_mode = "a"
        if self._written_years:
            last = self._written_years[-1]
            if self.year != last + 1:
                raise ValueError(
                    f"{os.path.basename(out_path)} already holds years "
                    f"{self._written_years[0]}-{last}, so the next year to "
                    f"accumulate is {last + 1}, but this run resumes inside "
                    f"year {self.year}: writing it would rename a year of the "
                    f"submitted series. Resume from a checkpoint inside "
                    f"{last + 1}, or move the annual file aside."
                )
            self.log(f"  ISMIP7 output: appending to {os.path.basename(out_path)} "
                     f"({len(self._written_years)} years through {last}; "
                     f"{self.year_time:.2f} yr of {self.year} carried over)")
        if self.comm.rank == 0:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            new = not os.path.exists(scalars_path)
            self._csv = open(scalars_path, "a")
            if new:
                self._csv.write("year," + ",".join(SCALARS) + "\n"); self._csv.flush()
        else:
            self._csv = None

    # ---- per-step bookkeeping -------------------------------------------
    def begin_step(self):
        for a in self.step_acc.values():
            a[:] = 0.0
        self.step_time = 0.0

    def book_advance(self, dt, accum, ocean_melt, a_ref, h_dg, u, grounded_cells):
        r"""Called by ``_advance`` after the transport solve, BEFORE removal:
        books the sources REQUESTED of this advance (the forcing SMB before
        the positivity limiter, and the melt) and the grounding-line flux
        with the velocity the transport used."""
        smb = assemble(accum * self._phi * dx).dat.data_ro / self.cell_area        # m/yr, cell mean
        melt = assemble(ocean_melt * self._phi * dx).dat.data_ro / self.cell_area
        self.step_acc["acabf"] += smb * dt
        if a_ref is not None:
            corr = assemble(a_ref * self._phi * dx).dat.data_ro / self.cell_area
            self.step_acc["acabf_correction"] += corr * dt
        self.step_acc["libmassbffl"] += -melt * dt * (~grounded_cells)
        # grounding-line flux into the first floating cell: upwind facet flux
        # across facets whose two cells differ in grounding, booked to the
        # floating side (its test function)
        g = Function(self.Q_dg); g.dat.data[:] = grounded_cells.astype(float)
        un = fd.dot(u, self._n); un_plus = (un + abs(un)) / 2
        phi = self._phi
        form = ((un_plus("+") * h_dg("+") * g("+") * (1 - g("-")) * phi("-")
                 + un_plus("-") * h_dg("-") * g("-") * (1 - g("+")) * phi("+")) * dS)
        assemble(form, tensor=self._gl_cof)
        self.step_acc["ligroundf"] += self._gl_cof.dat.data_ro / self.cell_area * dt   # m/yr equivalent
        self.step_time += dt

    def book_removal(self, cells, thickness_removed):
        r"""Ice removed at the front this advance (m of thickness per cell)."""
        self.step_acc["licalvf"][cells] -= thickness_removed

    def commit_step(self):
        for k in self.year_acc:
            self.year_acc[k] += self.step_acc[k]
        self.year_time += self.step_time

    # ---- carried across a chained resume ----------------------------------
    def state_fields(self):
        r"""The year in progress as DG0 Functions, for the run's checkpoint.

        Firedrake redistributes these on load, so a chained link may run on a
        different rank count than the one that wrote them."""
        fields = {}
        for k in self.ACCUMULATORS:
            f = Function(self.Q_dg, name=self.STATE_PREFIX + k)
            f.dat.data[:] = self.year_acc[k]
            fields[f.name()] = f
        h0 = Function(self.Q_dg, name=self.STATE_THICKNESS)
        if self.h_year_start is not None:
            h0.dat.data[:] = self.h_year_start
        fields[h0.name()] = h0
        return fields

    def state_attrs(self):
        r"""The year being accumulated and how much of it is in the sums."""
        return {self.STATE_YEAR: int(self.year),
                self.STATE_YEAR_TIME: float(self.year_time)}

    @classmethod
    def read_state(cls, chk, mesh):
        r"""The counterpart of ``state_fields``/``state_attrs``: the resume
        dict to hand back to the constructor, or None if the checkpoint was
        written by a run without ISMIP7 output."""
        if not chk.has_attr("/", cls.STATE_YEAR):
            return None
        return {
            "year": int(chk.get_attr("/", cls.STATE_YEAR)),
            "year_time": float(chk.get_attr("/", cls.STATE_YEAR_TIME)),
            "acc": {k: chk.load_function(mesh, name=cls.STATE_PREFIX + k).dat.data_ro.copy()
                    for k in cls.ACCUMULATORS},
            "h_year_start": chk.load_function(mesh, name=cls.STATE_THICKNESS).dat.data_ro.copy(),
        }

    # ---- year end ---------------------------------------------------------
    def start_year(self, h_dg):
        self.h_year_start = h_dg.dat.data_ro.copy()
        for a in self.year_acc.values():
            a[:] = 0.0
        self.year_time = 0.0

    def year_end(self, h_dg, s, b, u, tau, grounded_cells, ice_cells):
        r"""Write this year's state snapshot and flux means, then start the next."""
        yr = self.year
        T = self.year_time if self.year_time > 0 else 1.0
        fields = {}
        Q = self.Q_dg
        def dg(arr):
            f = Function(Q); f.dat.data[:] = arr; return f
        # thickness only where the mask says ice (h > 1 m): the checker
        # requires lithk == 0 wherever sftgif == 0, and sub-metre inflow in
        # buffer cells is the front's bookkeeping, not ice
        fields["lithk"] = dg(h_dg.dat.data_ro * ice_cells)
        fields["orog"] = Function(Q).interpolate(s)
        fields["topg"] = Function(Q).interpolate(b)
        fields["base"] = dg(fields["orog"].dat.data_ro - h_dg.dat.data_ro)
        fields["sftgif"] = dg(ice_cells.astype(float))
        fields["sftgrf"] = dg((ice_cells & grounded_cells).astype(float))
        fields["sftflf"] = dg((ice_cells & ~grounded_cells).astype(float))
        fields["strbasemag"] = Function(Q).interpolate(fd.sqrt(fd.inner(tau, tau)))   # MPa
        ux = Function(Q).interpolate(u[0]); uy = Function(Q).interpolate(u[1])
        for name in ("xvelmean", "xvelsurf", "xvelbase"):
            fields[name] = ux
        for name in ("yvelmean", "yvelsurf", "yvelbase"):
            fields[name] = uy
        for name in ("acabf", "acabf_correction", "libmassbffl", "licalvf", "ligroundf"):
            fields[name] = dg(self.year_acc[name] / T)                            # m/yr ice
        fields["libmassbfgr"] = dg(np.zeros_like(self.cell_area))
        fields["lifmassbf"] = dg(np.zeros_like(self.cell_area))
        h0 = self.h_year_start if self.h_year_start is not None else h_dg.dat.data_ro
        fields["dlithkdt"] = dg((h_dg.dat.data_ro - h0) / T)                        # m/yr
        # idx is a 0-based position, not the year: Firedrake sizes the
        # timestepping dataset by the largest idx, so idx=2016 allocated 2017
        # rows per field (2.5 GB for two years of a 7543-cell mesh); the
        # years attribute maps position to year.
        with fd.CheckpointFile(self.out_path, self._chk_mode) as chk:
            if self._chk_mode == "w":
                chk.save_mesh(self.mesh)
            for name, f in fields.items():
                chk.save_function(f, name=name, idx=len(self._written_years))
            chk.set_attr("/", "years", ",".join(str(y) for y in self._written_years + [yr]))
        self._chk_mode = "a"
        self._written_years.append(yr)
        # scalars, from the same fields (kg, m2, kg/s)
        area = self.cell_area
        def integ(arr):
            return self.comm.allreduce(float((arr * area).sum()))
        # limnsw is the mass of the ice ABOVE FLOTATION: the request defines it
        # as that volume times the ice density, so the integrand is the
        # thickness above flotation, h - h_f with h_f = max(-b, 0) / rho_ratio,
        # not the height above flotation s - s_float (which is rho_ratio times
        # smaller on marine beds and wrong outright where the bed is dry).
        topg = fields["topg"].dat.data_ro
        haf_thickness = np.maximum(
            h_dg.dat.data_ro - np.maximum(-topg, 0.0) / self.rho_ratio, 0.0)
        row = {
            "lim": integ(h_dg.dat.data_ro) * RHO_I,
            "limnsw": integ(haf_thickness * grounded_cells) * RHO_I,
            "iareagr": integ((ice_cells & grounded_cells).astype(float)),
            "iareafl": integ((ice_cells & ~grounded_cells).astype(float)),
            "tendacabf": integ(fields["acabf"].dat.data_ro) * RHO_I / SECONDS_PER_YEAR,
            "tendlibmassbfgr": 0.0,
            "tendlibmassbffl": integ(fields["libmassbffl"].dat.data_ro) * RHO_I / SECONDS_PER_YEAR,
            "tendlicalvf": integ(fields["licalvf"].dat.data_ro) * RHO_I / SECONDS_PER_YEAR,
            "tendlifmassbf": 0.0,
            "tendligroundf": integ(fields["ligroundf"].dat.data_ro) * RHO_I / SECONDS_PER_YEAR,
        }
        if self._csv is not None:
            self._csv.write(f"{yr}," + ",".join(f"{row[k]:.6e}" for k in SCALARS) + "\n"); self._csv.flush()
        self.log(f"  ISMIP7 output: year {yr} written ({len(fields)} fields; GL flux "
                 f"{row['tendligroundf'] * SECONDS_PER_YEAR / 1e12:+.0f} Gt/yr, calving "
                 f"{row['tendlicalvf'] * SECONDS_PER_YEAR / 1e12:+.0f} Gt/yr)")
        self.year = yr + 1
        self.start_year(h_dg)

    def close(self):
        if self._csv is not None:
            self._csv.close()
