r"""The writer loads the model mesh once per series and writes what loading
it every year wrote.

Loading the mesh was most of the writer's cost, 92 to 95 s of about 104 s a
year on the 2 km production mesh, so ``write_ismip7_output.SeriesMesh``
loads it from the first year and reads a later year onto it when the year
file's ``mesh_signature`` matches the first's, and onto the year's own mesh
when it does not. These run the writer end to end on a series the forward's
own ``AnnualOutput`` banks on a 32-cell mesh of two by two 8 km pixels,
once as it runs by default and once with ``--mesh-per-year``, and hold
every array and attribute of the two trees equal.
"""
import os
import sys

import numpy as np
import pytest

fd = pytest.importorskip("firedrake")
netCDF4 = pytest.importorskip("netCDF4")
sp = pytest.importorskip("scipy.sparse")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "antarctica", "scripts"))
wio = pytest.importorskip("write_ismip7_output")

from firedrake.mesh import plex_from_cell_list  # noqa: E402

from icepack2_tools.ismip7_output import AnnualOutput  # noqa: E402

RHO_RATIO = 917.0 / 1024.0
DX = wio.ISMIP7_DX
# two by two pixels of the 8 km grid, edges on pixel edges, so every cell
# lies inside one pixel and the overlap operator is each cell's area
X_LO = wio.ISMIP7_X0 - DX / 2 + 300 * DX
Y_LO = wio.ISMIP7_Y0 - DX / 2 + 200 * DX
STEM = "hist_test_2000_ismip7_annual.h5"


def _mesh(faces=(4, 4), name="firedrake_default"):
    r"""32 triangles over the four pixels. The topology carries a fixed
    name, as a mesh read from a .msh or checkpoint file does, so two meshes
    built alike are stored alike and two built differently (``faces``) are
    told apart by their content, the names being equal."""
    nx, ny = faces
    X, Y = np.meshgrid(X_LO + 2 * DX * np.arange(nx + 1) / nx,
                       Y_LO + 2 * DX * np.arange(ny + 1) / ny, indexing="ij")
    cells = []
    for i in range(nx):
        for j in range(ny):
            v00, v10 = i * (ny + 1) + j, (i + 1) * (ny + 1) + j
            cells += [[v00, v10, v10 + 1], [v00, v10 + 1, v00 + 1]]
    plex = plex_from_cell_list(2, cells, np.column_stack([X.ravel(), Y.ravel()]),
                               fd.COMM_WORLD, name + "_topology")
    return fd.Mesh(plex, name=name)


def _bank(out, mesh, years, seed, uniform=False):
    r"""Years of a forward on ``mesh`` through ``AnnualOutput``: one booked
    advance per year and its year end. ``uniform`` grounds ice on every
    cell, whose area sums then hold on any mesh of the four pixels."""
    Q = fd.FunctionSpace(mesh, "DG", 0)
    V = fd.VectorFunctionSpace(mesh, "CG", 1)
    rng = np.random.default_rng(seed)
    n = Q.dim()
    bed = rng.uniform(-700.0, 150.0, n) if not uniform else np.full(n, 100.0)
    annual = AnnualOutput(mesh, Q, os.path.join(out, STEM),
                          os.path.join(out, STEM.replace("_annual.h5", "_scalars.csv")),
                          first_year=years[0], rho_ratio=RHO_RATIO, bed=bed)

    def dg(values):
        f = fd.Function(Q)
        f.dat.data[:] = values
        return f

    h = dg(rng.uniform(0.0, 1500.0, n))
    annual.start_year(h)
    x, y = fd.SpatialCoordinate(mesh)
    for _ in years:
        thick = rng.uniform(0.0, 1500.0, n)
        if not uniform:
            thick[rng.choice(n, 4, replace=False)] = 0.0          # ice-free cells
        else:
            thick += 10.0
        h = dg(thick)
        grounded = (thick * RHO_RATIO > -bed) | (bed >= 0.0)
        s = dg(np.where(grounded, bed + thick, (1.0 - RHO_RATIO) * thick))
        u = fd.Function(V).interpolate(fd.as_vector(
            (100.0 + rng.uniform(-50.0, 50.0) * (x - X_LO) / DX,
             -30.0 + rng.uniform(-20.0, 20.0) * (y - Y_LO) / DX)))
        tau = fd.Function(V).interpolate(fd.as_vector((0.05 * (x - X_LO) / DX, 0.02)))
        annual.begin_step()
        annual.book_advance(1.0, dg(rng.uniform(-0.5, 1.0, n)), dg(rng.uniform(0.0, 5.0, n)),
                            None, h, u, grounded)
        annual.commit_step()
        annual.year_end(h, s, dg(bed), u, tau, grounded, thick > AnnualOutput.ICE_THICKNESS)
    annual.close()


def _operator(out, first_year):
    r"""The overlap operator of the first year's mesh, cached where the
    writer looks for it: each cell's area in the pixel holding it."""
    with fd.CheckpointFile(AnnualOutput.year_path(os.path.join(out, STEM), first_year), "r") as chk:
        mesh = chk.load_mesh()
    Q = fd.FunctionSpace(mesh, "DG", 0)
    area = fd.assemble(fd.TestFunction(Q) * fd.dx).dat.data_ro
    xc = [fd.Function(Q).interpolate(fd.SpatialCoordinate(mesh)[i]).dat.data_ro for i in (0, 1)]
    i = np.floor((xc[0] - (wio.ISMIP7_X0 - DX / 2)) / DX).astype(int)
    j = np.floor((xc[1] - (wio.ISMIP7_Y0 - DX / 2)) / DX).astype(int)
    W = sp.csr_matrix((area, (j * wio.ISMIP7_NX + i, np.arange(len(area)))),
                      shape=(wio.ISMIP7_NX * wio.ISMIP7_NY, len(area)))
    np.savez(os.path.join(out, STEM) + ".overlap.npz", data=W.data, indices=W.indices,
             indptr=W.indptr, shape=np.array(W.shape))


def _write(out, tree, monkeypatch, capsys, *extra):
    monkeypatch.setattr(sys, "argv", ["write_ismip7_output.py", os.path.join(out, STEM),
                                      "--out-dir", tree, "--esm", "CESM2-WACCM",
                                      "--scenario", "historical", "--source-id", "IU",
                                      "--ism-id", "icepack2", *extra])
    wio.main()
    return capsys.readouterr().out


def _files(tree):
    found = {}
    for root, _, names in os.walk(tree):
        for name in names:
            found[name] = os.path.join(root, name)
    return found


def _attrs(obj):
    return {k: obj.getncattr(k) for k in obj.ncattrs()}


def _same_tree(a, b):
    r"""Every file, variable, array and attribute of tree ``a`` in ``b``."""
    fa, fb = _files(a), _files(b)
    assert sorted(fa) == sorted(fb) and len(fa) == len(wio.VARIABLES_2D) + len(wio.SCALARS)
    for name in fa:
        with netCDF4.Dataset(fa[name]) as da, netCDF4.Dataset(fb[name]) as db:
            assert _attrs(da) == _attrs(db), name
            assert list(da.variables) == list(db.variables), name
            for var in da.variables:
                va, vb = da.variables[var], db.variables[var]
                va.set_auto_mask(False)
                vb.set_auto_mask(False)
                assert _attrs(va) == _attrs(vb), (name, var)
                assert va.dtype == vb.dtype and va.dimensions == vb.dimensions, (name, var)
                assert np.array_equal(va[:], vb[:]), (name, var)


def _signature(path):
    with fd.CheckpointFile(path, "r") as chk:
        return wio.mesh_signature(chk.h5pyfile)


def test_reading_onto_the_first_years_mesh_writes_what_per_year_meshes_write(
        tmp_path, monkeypatch, capsys):
    out = str(tmp_path / "run")
    _bank(out, _mesh(), [2015, 2016, 2017], seed=1)
    _operator(out, 2015)
    log = _write(out, str(tmp_path / "once"), monkeypatch, capsys)
    assert "mesh: loaded once, from 2015; 3 of 3 years read onto it" in log
    assert "onto their own" not in log
    log = _write(out, str(tmp_path / "per_year"), monkeypatch, capsys, "--mesh-per-year")
    assert "each of the 3 years read onto its own" in log
    _same_tree(str(tmp_path / "once"), str(tmp_path / "per_year"))

    # the series is not trivially empty: ice on the four pixels, another
    # thickness every year, and nothing elsewhere
    lithk, = [p for n, p in _files(str(tmp_path / "once")).items() if n.startswith("lithk_")]
    with netCDF4.Dataset(lithk) as ds:
        v = np.asarray(ds.variables["lithk"][:])
    four = v[:, 200:202, 300:302]
    assert v.shape[0] == 3 and (four > 0.0).all()
    assert not np.array_equal(four[0], four[1]) and not np.array_equal(four[1], four[2])
    assert float(np.abs(v).sum()) == pytest.approx(float(np.abs(four).sum()))


def test_a_year_whose_mesh_differs_is_read_onto_its_own(tmp_path, monkeypatch, capsys):
    r"""2016 comes from a mesh of as many cells under the same names but
    another topology, so its fields cannot be read through the first year's
    mesh. It is read onto its own, the run names it, and 2017, on the first
    year's mesh again, is read onto that."""
    out = str(tmp_path / "run")
    _bank(out, _mesh(), [2015], seed=2, uniform=True)
    _bank(out, _mesh(faces=(2, 8)), [2016], seed=3, uniform=True)
    _bank(out, _mesh(), [2017], seed=4, uniform=True)
    _operator(out, 2015)
    log = _write(out, str(tmp_path / "once"), monkeypatch, capsys)
    assert "year 2016: the file's mesh differs from 2015's" in log
    assert "2 of 3 years read onto it, 1 onto their own, whose files differ: 2016" in log
    _write(out, str(tmp_path / "per_year"), monkeypatch, capsys, "--mesh-per-year")
    _same_tree(str(tmp_path / "once"), str(tmp_path / "per_year"))


def test_the_signature_sees_the_mesh_and_not_the_fields(tmp_path):
    r"""Years of one run, and a run on a mesh built the same way, sign
    alike whatever their fields hold; another topology or a coordinate moved
    by a metre signs otherwise."""
    def signatures(name, mesh, seed, years):
        out = str(tmp_path / name)
        _bank(out, mesh, years, seed=seed)
        return [_signature(AnnualOutput.year_path(os.path.join(out, STEM), y)) for y in years]

    moved = _mesh()
    moved.coordinates.dat.data[0] += 1.0
    a = signatures("a", _mesh(), 5, [2015, 2016])
    b = signatures("b", _mesh(), 6, [2015])
    assert a[0] == a[1] == b[0]
    assert signatures("faces", _mesh(faces=(2, 8)), 7, [2015]) != b
    assert signatures("moved", moved, 8, [2015]) != b
