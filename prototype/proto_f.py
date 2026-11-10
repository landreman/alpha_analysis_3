"""End-to-end prototype of the ridge-based algorithm for the alpha loss-fraction metric f.

Pitch-independent stage: FFT tables per flux surface, ridge network (along-line maxima of B),
chains of extrema along field lines from every ridge node.
Per-pitch stage: active states, exits, J and K, transition detection, cells, sheets, action
bins, graph reachability from s=1, K-weighted reachable weight, ledger check.
Outer stage: f from a coarse pitch grid, for several source profiles h(rho).

Deliberately simple (no Newton polishing, no bisection of cut edges, no adaptive refinement);
it tests the structure of the algorithm, not its final accuracy.
"""
from __future__ import annotations
import sys, time, os
import numpy as np
import contourpy
from scipy.ndimage import map_coordinates, spline_filter
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from alpha_analysis_3.boozer_field import BoozerField

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..') + '/'
FILES = {
 'PCA':  'data/boozmn_20260402-01-038_Ax_PCA_20dofs_allNfp_aspect6_eval000290_low_resolution.nc',
 'TURBO':'data/boozmn_20260402-01-178_TURBO_Garabedian_mpol1_xmin0p1_allNfp_aspect6_eval000155.nc',
 'DMerc':'data/boozmn_20260406-01-262-Ax_nfp4_Garabedian_mpol2_ntor2_minx0_allNfp_aspect10_DMercFail_m0p3_eval000323_low_resolution.nc',
 'd23p4':'data/boozmn_d23p4_tm_ns51_mbooz16_nbooz16.nc',
 'n3are':'data/boozmn_n3are_R7.75B5.7_mbooz18_nbooz12.nc',
}
TWO_PI = 2 * np.pi
SOURCES = {'h=1': lambda rho: np.ones_like(np.asarray(rho, float)),
           'h=(1-rho^2)^4': lambda rho: (1.0 - np.asarray(rho, float) ** 2) ** 4}
trapz = np.trapezoid if hasattr(np, 'trapezoid') else np.trapz


# ----------------------------------------------------------------------------- field tables
class Surface:
    """FFT tables of B, D_par B, D_par^2 B on the (theta, phi) torus at one s; phi = nfp*zeta."""

    def __init__(self, fld: BoozerField, s: float, NT=512, NZ=512):
        self.s = s; self.NT = NT; self.NZ = NZ
        self.nfp = fld.nfp; self.L = TWO_PI / fld.nfp
        self.iota = float(fld.iota(s)); self.G = float(fld.G(s)); self.I = float(fld.I(s))
        self.C = self.G + self.iota * self.I
        self.sgn = 1.0 if self.C >= 0 else -1.0
        c = np.asarray(fld.bmnc(s)).ravel()
        m = fld.xm.astype(int); n = fld.xn.astype(int)
        k = m * self.iota - n                        # along-line wavenumber (per radian of zeta)
        nprime = n // self.nfp
        A = np.zeros((NT, NZ), complex); AD = np.zeros_like(A); AD2 = np.zeros_like(A)
        np.add.at(A, (m % NT, (-nprime) % NZ), c)
        np.add.at(AD, (m % NT, (-nprime) % NZ), 1j * k * c)
        np.add.at(AD2, (m % NT, (-nprime) % NZ), -(k ** 2) * c)
        scale = NT * NZ
        self.B = np.real(np.fft.ifft2(A)) * scale
        self.DB = np.real(np.fft.ifft2(AD)) * scale
        self.D2B = np.real(np.fft.ifft2(AD2)) * scale
        self.Bmin = self.B.min(); self.Bmax = self.B.max()
        self._Bc = spline_filter(self.B, order=3, mode='grid-wrap')
        self._DBc = spline_filter(self.DB, order=3, mode='grid-wrap')
        self._D2Bc = spline_filter(self.D2B, order=3, mode='grid-wrap')
        self.inv_B2_mean = np.mean(1.0 / self.B ** 2)

    def _interp(self, coeffs, theta, phi):
        th = np.asarray(theta, float) / TWO_PI * self.NT
        ph = np.asarray(phi, float) / TWO_PI * self.NZ
        return map_coordinates(coeffs, [th.ravel(), ph.ravel()], order=3, mode='grid-wrap',
                               prefilter=False).reshape(th.shape)

    def b(self, theta, phi): return self._interp(self._Bc, theta, phi)
    def db(self, theta, phi): return self._interp(self._DBc, theta, phi)
    def d2b(self, theta, phi): return self._interp(self._D2Bc, theta, phi)


# ----------------------------------------------------------------------------- ridge network
def _stitch(lines, tol=1e-6):
    """Join contourpy polylines across the periodic seams of the [0,2pi]^2 domain into curves."""
    pieces = [np.array(l) for l in lines if len(l) >= 2]
    closed, open_ = [], []
    for p in pieces:
        if np.hypot(*(p[0] - p[-1])) < tol:
            closed.append(p)
        else:
            open_.append(p)

    def same(a, b):
        d = (a - b + np.pi) % TWO_PI - np.pi
        return np.hypot(*d) < 1e-3

    used = [False] * len(open_)
    for i in range(len(open_)):
        if used[i]: continue
        used[i] = True
        cur = open_[i].copy()
        while True:
            end = cur[-1]
            found = False
            for j in range(len(open_)):
                if used[j]: continue
                q = open_[j]
                for flip in (False, True):
                    qq = q[::-1] if flip else q
                    if same(qq[0], end):
                        shift = np.round((end - qq[0]) / TWO_PI) * TWO_PI
                        cur = np.vstack([cur, (qq + shift)[1:]])
                        used[j] = True; found = True; break
                if found: break
            if not found: break
            if same(cur[-1], cur[0]):
                break
        closed.append(cur)
    return closed


def ridge_arcs(surf: Surface):
    """Maximum arcs of D_par B = 0 on the surface: dict(theta, phi unwrapped, closed, wind=(p,q))."""
    NT, NZ = surf.NT, surf.NZ
    DBe = np.pad(surf.DB, ((0, 1), (0, 1)), mode='wrap')
    th = np.linspace(0, TWO_PI, NT + 1); ph = np.linspace(0, TWO_PI, NZ + 1)
    gen = contourpy.contour_generator(x=ph, y=th, z=DBe)   # x=phi, y=theta
    lines = [np.c_[l[:, 1], l[:, 0]] for l in gen.lines(0.0)]  # -> (theta, phi)
    curves = _stitch(lines)
    arcs = []
    for c in curves:
        th_c, ph_c = c[:, 0], c[:, 1]
        d2 = surf.d2b(th_c, ph_c)
        ismax = d2 < 0
        closed_whole = np.hypot(*((c[-1] - c[0] + np.pi) % TWO_PI - np.pi)) < 1e-3
        if ismax.all() and closed_whole:
            p = int(np.round((th_c[-1] - th_c[0]) / TWO_PI)); q = int(np.round((ph_c[-1] - ph_c[0]) / TWO_PI))
            if p != 0 or q != 0:
                arcs.append(dict(theta=th_c[:-1], phi=ph_c[:-1], closed=True, wind=(p, q)))
                continue
        if not ismax.any():
            continue
        n = len(c)
        if closed_whole:
            k0 = int(np.argmin(ismax))
            idx = np.r_[k0:n - 1, 0:k0]
            th_r, ph_r, m_r = np.unwrap(th_c[idx]), np.unwrap(ph_c[idx]), ismax[idx]
        else:
            th_r, ph_r, m_r = th_c, ph_c, ismax
        start = None
        for i, flag in enumerate(np.r_[m_r, False]):
            if flag and start is None: start = i
            if (not flag) and start is not None:
                if i - start >= 3:
                    arcs.append(dict(theta=th_r[start:i], phi=ph_r[start:i], closed=False, wind=(0, 0)))
                start = None
    return arcs


def place_nodes(arcs, n_closed=256, n_open=48):
    """Resample arcs to nodes: closed curves uniformly in the winding angle, open arcs in arc length."""
    out = []
    for a in arcs:
        th, ph = a['theta'], a['phi']
        if a['closed']:
            p, q = a['wind']
            ang = th / p if p != 0 else ph / q
            per = TWO_PI
            mono_up = np.all(np.diff(ang) > 0); mono_dn = np.all(np.diff(ang) < 0)
            if mono_up or mono_dn:
                u = (ang - ang[0]) % per if mono_up else (ang[0] - ang) % per
                order = np.argsort(u)
                u, th, ph = u[order], th[order], ph[order]
                uu = np.linspace(0, per, n_closed, endpoint=False)
                u_ext = np.r_[u, u[0] + per]
                th_ext = np.r_[th, th[0] + TWO_PI * p]
                ph_ext = np.r_[ph, ph[0] + TWO_PI * q]
                out.append(dict(theta=np.interp(uu, u_ext, th_ext), phi=np.interp(uu, u_ext, ph_ext),
                                closed=True, wind=(p, q), u=uu / per))
            else:
                seg = np.hypot(np.diff(th), np.diff(ph)); s_arc = np.r_[0, np.cumsum(seg)]
                tot = s_arc[-1] + np.hypot(th[0] + TWO_PI * p - th[-1], ph[0] + TWO_PI * q - ph[-1])
                uu = np.linspace(0, tot, n_closed, endpoint=False)
                out.append(dict(theta=np.interp(uu, np.r_[s_arc, tot], np.r_[th, th[0] + TWO_PI * p]),
                                phi=np.interp(uu, np.r_[s_arc, tot], np.r_[ph, ph[0] + TWO_PI * q]),
                                closed=True, wind=(p, q), u=uu / tot))
        else:
            seg = np.hypot(np.diff(th), np.diff(ph)); s_arc = np.r_[0, np.cumsum(seg)]
            if s_arc[-1] < 1e-3: continue
            uu = np.linspace(0, s_arc[-1], n_open)
            out.append(dict(theta=np.interp(uu, s_arc, th), phi=np.interp(uu, s_arc, ph), closed=False,
                            wind=(0, 0), u=uu / s_arc[-1]))
    return out


# ----------------------------------------------------------------------------- chains
def chains(surf: Surface, theta0, phi0, periods=48, nper=128):
    """Follow field lines forward from ridge nodes; return extrema arrays (NaN padded).

    max_z/max_H: maxima k=1,2,.. after the node (array index k-1); min_z/min_H: minima, with
    min index j lying between max_j and max_{j+1} (max_0 is the node itself at offset 0)."""
    N = len(theta0)
    dz = surf.sgn * surf.L / nper
    k = np.arange(periods * nper + 1)
    zeta = k[None, :] * dz
    theta = theta0[:, None] + surf.iota * zeta
    phi = phi0[:, None] + surf.nfp * zeta
    Bl = surf.b(theta, phi)
    d = np.diff(Bl, axis=1)
    is_max = np.zeros_like(Bl, bool); is_min = np.zeros_like(Bl, bool)
    is_max[:, 1:-1] = (d[:, :-1] > 0) & (d[:, 1:] <= 0)
    is_min[:, 1:-1] = (d[:, :-1] < 0) & (d[:, 1:] >= 0)
    is_max[:, :3] = False; is_min[:, :2] = False

    def collect(mask):
        counts = mask.sum(axis=1); KM = int(counts.max()) if N else 0
        Z = np.full((N, KM), np.nan); H = np.full((N, KM), np.nan)
        rows, cols = np.nonzero(mask)
        y0 = Bl[rows, cols - 1]; y1 = Bl[rows, cols]; y2 = Bl[rows, cols + 1]
        denom = (y0 - 2 * y1 + y2)
        with np.errstate(divide='ignore', invalid='ignore'):
            off = np.where(np.abs(denom) > 1e-300, 0.5 * (y0 - y2) / denom, 0.0)
        off = np.clip(off, -1, 1)
        hval = y1 - 0.25 * (y0 - y2) * off
        zval = (cols + off) * dz
        order = np.lexsort((cols, rows))
        rows, zval, hval = rows[order], zval[order], hval[order]
        first = np.r_[True, rows[1:] != rows[:-1]]
        idx = np.arange(len(rows)); starts = np.where(first)[0]
        pos = idx - np.repeat(starts, np.diff(np.r_[starts, len(rows)]))
        Z[rows, pos] = zval; H[rows, pos] = hval
        return Z, H

    max_z, max_H = collect(is_max)
    min_z, min_H = collect(is_min)
    # pad both extremum tables to a common width so index arithmetic between them is safe
    KM = max(max_z.shape[1], min_z.shape[1]) + 1

    def pad(a):
        out = np.full((a.shape[0], KM), np.nan); out[:, :a.shape[1]] = a; return out
    return dict(max_z=pad(max_z), max_H=pad(max_H), min_z=pad(min_z), min_H=pad(min_H), dz=dz)


# ----------------------------------------------------------------------------- per-pitch wells
GL_T, GL_W = np.polynomial.legendre.leggauss(384)
GL_T = 0.5 * np.pi * GL_T; GL_W = 0.5 * np.pi * GL_W      # t in [-pi/2, pi/2]
SIN_T = np.sin(GL_T); COS_T = np.cos(GL_T)


def wells_at_pitch(surf: Surface, nodes_theta, nodes_phi, node_H, ch, b, n_bisect=40):
    """For every ridge node, the well at pitch b (if any): exits, bounce points, J, K."""
    N = len(nodes_theta)
    max_H, max_z, min_H, min_z = ch['max_H'], ch['max_z'], ch['min_H'], ch['min_z']
    m0 = np.where(np.isnan(min_H[:, 0]), np.inf, min_H[:, 0])
    active = (node_H > b) & (m0 < b)
    above = np.where(np.isnan(max_H), False, max_H > b)
    has_exit = above.any(axis=1)
    n_exit = np.where(has_exit, np.argmax(above, axis=1), -1)
    # the minimum preceding the exit maximum must exist (clean alternation); otherwise unknown
    ok_min = has_exit & np.isfinite(min_z[np.arange(N), np.clip(n_exit, 0, None)])
    has_exit = ok_min
    unknown = active & ~has_exit
    act = active & has_exit
    idx = np.nonzero(act)[0]
    out = dict(active=active, unknown=unknown, n_exit=n_exit, J=np.full(N, np.nan), K=np.full(N, np.nan),
               z_in=np.full(N, np.nan), z_out=np.full(N, np.nan), z_exit=np.full(N, np.nan),
               exit_gap=np.full(N, np.nan),
               barrier_H=np.full(N, np.nan), barrier_z=np.full(N, np.nan), barrier_k=np.full(N, -1))
    if len(idx) == 0:
        return out
    ne = n_exit[idx]
    z_exit = max_z[idx, ne]
    # spacing from the exit maximum to its neighbouring maxima along the line (identity scale)
    z_prev = np.where(ne > 0, max_z[idx, np.clip(ne - 1, 0, None)], 0.0)
    nxt = np.clip(ne + 1, 0, max_z.shape[1] - 1)
    z_next = max_z[idx, nxt]
    gap_prev = np.abs(z_exit - z_prev)
    gap_next = np.where(np.isfinite(z_next) & (nxt != ne), np.abs(z_next - z_exit), np.inf)
    out['exit_gap'][idx] = np.minimum(gap_prev, gap_next)
    z_min_last = min_z[idx, ne]
    KM = max_H.shape[1]
    kk = np.arange(KM)[None, :]
    interior = (kk < ne[:, None])
    Hint = np.where(interior, max_H[idx], -np.inf)
    kb = np.argmax(Hint, axis=1)
    has_int = ne > 0
    out['barrier_H'][idx] = np.where(has_int, Hint[np.arange(len(idx)), kb], np.nan)
    out['barrier_z'][idx] = np.where(has_int, max_z[idx, kb], np.nan)
    out['barrier_k'][idx] = np.where(has_int, kb, -1)

    th0, ph0 = nodes_theta[idx], nodes_phi[idx]

    def B_at(zoff):
        zoff = np.asarray(zoff)
        if zoff.ndim == 2:
            return surf.b(th0[:, None] + surf.iota * zoff, ph0[:, None] + surf.nfp * zoff)
        return surf.b(th0 + surf.iota * zoff, ph0 + surf.nfp * zoff)

    lo = np.zeros(len(idx)); hi = min_z[idx, 0].copy()
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi); gt = B_at(mid) > b
        lo = np.where(gt, mid, lo); hi = np.where(gt, hi, mid)
    z_in = 0.5 * (lo + hi)
    lo = z_min_last.copy(); hi = z_exit.copy()
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi); lt = B_at(mid) < b
        lo = np.where(lt, mid, lo); hi = np.where(lt, hi, mid)
    z_out = 0.5 * (lo + hi)
    mid = 0.5 * (z_in + z_out); half = 0.5 * (z_out - z_in)
    zq = mid[:, None] + half[:, None] * SIN_T[None, :]
    Bq = B_at(zq)
    g = np.clip(1.0 - Bq / b, 0.0, None)
    sq = np.sqrt(g)
    dl = np.abs(surf.C) / Bq * np.abs(half[:, None]) * COS_T[None, :]
    J = np.sum(GL_W[None, :] * dl * sq, axis=1)
    with np.errstate(divide='ignore', invalid='ignore'):
        K = np.sum(GL_W[None, :] * dl / np.where(sq > 0, sq, np.inf), axis=1)
    out['J'][idx] = J; out['K'][idx] = K
    out['z_in'][idx] = z_in; out['z_out'][idx] = z_out; out['z_exit'][idx] = z_exit
    return out


# ----------------------------------------------------------------------------- geometry assembly
class Geometry:
    """Pitch-independent data for one equilibrium: surfaces, arcs (families), nodes, chains."""

    def __init__(self, name, n_s=32, s_min=0.02, s_max=1.0, n_closed=256, n_open=48, periods=48,
                 NT=512, NZ=512, verbose=True):
        t0 = time.time()
        self.name = name
        self.fld = BoozerField(ROOT + FILES[name])
        self.nfp = self.fld.nfp; self.L = TWO_PI / self.nfp
        self.s = np.linspace(s_min, s_max, n_s)
        self.ds = self.s[1] - self.s[0]
        self.surfs = []; self.arcs = []
        self.Bmin_g = np.inf; self.Bmax_g = -np.inf
        n_nodes = 0
        for s in self.s:
            surf = Surface(self.fld, s, NT, NZ)
            arcs = place_nodes(ridge_arcs(surf), n_closed, n_open)
            for a in arcs:
                a['H'] = surf.b(a['theta'], a['phi'])
                a['alpha'] = a['theta'] - surf.iota * a['phi'] / self.nfp
                dal = np.diff(a['alpha'])
                if a['closed']:
                    p, q = a['wind']
                    total = TWO_PI * p - surf.iota * self.L * q
                    dal = np.r_[dal, a['alpha'][0] + total - a['alpha'][-1]]
                a['dalpha'] = np.abs(dal)
                a['ch'] = chains(surf, a['theta'], a['phi'], periods=periods)
                n_nodes += len(a['theta'])
            self.surfs.append(surf); self.arcs.append(arcs)
            self.Bmin_g = min(self.Bmin_g, surf.Bmin); self.Bmax_g = max(self.Bmax_g, surf.Bmax)
        self.match_families()
        self.t_geom = time.time() - t0
        self.n_nodes = n_nodes
        if verbose:
            kinds = {}
            for fam in self.families:
                a0 = self.arcs[fam[0][0]][fam[0][1]]
                key = ('closed ' + str(a0['wind'])) if a0['closed'] else 'open'
                kinds[key] = kinds.get(key, 0) + 1
            print(f'[{name}] geometry: {n_s} surfaces (s_min={s_min}), {n_nodes} ridge nodes, {len(self.families)} arc families '
                  f'{kinds}, B in [{self.Bmin_g:.4f},{self.Bmax_g:.4f}], {self.t_geom:.1f} s')

    def match_families(self):
        fam_of = {}; families = []
        for i in range(len(self.s)):
            for j, a in enumerate(self.arcs[i]):
                best = None
                if i > 0:
                    for jp, ap in enumerate(self.arcs[i - 1]):
                        if ap['closed'] != a['closed'] or ap['wind'] != a['wind'] or len(ap['theta']) != len(a['theta']):
                            continue
                        dth = (a['theta'] - ap['theta'] + np.pi) % TWO_PI - np.pi
                        dph = (a['phi'] - ap['phi'] + np.pi) % TWO_PI - np.pi
                        dist = np.mean(np.hypot(dth, dph))
                        if dist < 0.5 and (best is None or dist < best[0]):
                            best = (dist, jp)
                if best is not None and (i - 1, best[1]) in fam_of:
                    fid = fam_of[(i - 1, best[1])]
                    if families[fid][-1][0] == i - 1:
                        fam_of[(i, j)] = fid; families[fid].append((i, j)); continue
                fam_of[(i, j)] = len(families); families.append([(i, j)])
        self.families = families; self.fam_of = fam_of

    def denominator(self, h):
        vals = np.array([h(np.sqrt(sf.s)) * abs(sf.C) * sf.inv_B2_mean * TWO_PI * self.L for sf in self.surfs])
        return trapz(vals, self.s)

    def trapped_weight_volume(self, b, h):
        """Independent ledger: int ds dtheta dzeta h |C| /(B sqrt(1-B/b)) [B<b][Bmax(s)>b] (dense-line)."""
        vals = []
        for sf in self.surfs:
            if sf.Bmax <= b:
                vals.append(0.0); continue
            g = 1.0 - sf.B / b
            m = g > 0
            with np.errstate(divide='ignore', invalid='ignore'):
                integrand = np.where(m, 1.0 / (sf.B * np.sqrt(np.where(m, g, 1.0))), 0.0)
            vals.append(h(np.sqrt(sf.s)) * abs(sf.C) * integrand.mean() * TWO_PI * self.L)
        return trapz(np.array(vals), self.s)


# ----------------------------------------------------------------------------- per-pitch slice
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0)]


def slice_at_pitch(geom: Geometry, b, jump_tol=0.1, bins_per_cell=1.0, verbose=False, sources=SOURCES):
    """Q(b) lower/upper (K-weighted reachable weight), Q_total, Q_unknown per source, plus diagnostics."""
    t0 = time.time()
    n_s = len(geom.s)
    W = {}
    for i in range(n_s):
        sf = geom.surfs[i]
        for j, a in enumerate(geom.arcs[i]):
            W[(i, j)] = wells_at_pitch(sf, a['theta'], a['phi'], a['H'], a['ch'], b)

    # per-surface transition points (u-edge exit jumps) vs Gamma_max crossings; multi-jump edges
    per_surface = []; n_multi = 0
    for i in range(n_s):
        nT = nG = 0
        for j, a in enumerate(geom.arcs[i]):
            w = W[(i, j)]
            good = w['active'] & ~w['unknown'] & np.isfinite(w['z_out'])
            z = w['z_out']; ne = w['n_exit']; hb = a['H'] - b
            if a['closed']:
                gp = good & np.roll(good, -1); dz = np.abs(z - np.roll(z, -1)); dn = np.abs(ne - np.roll(ne, -1))
                nG += int(np.sum(hb * np.roll(hb, -1) < 0))
            else:
                gp = good[:-1] & good[1:]; dz = np.abs(z[:-1] - z[1:]); dn = np.abs(ne[:-1] - ne[1:])
                nG += int(np.sum(hb[:-1] * hb[1:] < 0))
            nT += int(np.sum(gp & (dz > jump_tol)))
            n_multi += int(np.sum(gp & (dn >= 2)))
        per_surface.append((nT, nG))
    n_T_points = sum(t for t, g in per_surface); n_gamma = sum(g for t, g in per_surface)
    mismatch_surfaces = [(round(float(geom.s[i]), 3), t, g) for i, (t, g) in enumerate(per_surface) if t != g]

    # cells
    cells = []; transitions = []; unknown_cells = []
    n_T_u = n_T_s = 0
    for fid, fam in enumerate(geom.families):
        for (i, j), (i2, j2) in zip(fam[:-1], fam[1:]):
            a, a2 = geom.arcs[i][j], geom.arcs[i2][j2]
            nk = len(a['theta']); closed = a['closed']
            nK = nk if closed else nk - 1
            s_mid = 0.5 * (geom.s[i] + geom.s[i2])
            for k in range(nK):
                k1 = (k + 1) % nk
                ids = [(i, j, k), (i, j, k1), (i2, j2, k1), (i2, j2, k)]
                act = np.array([W[(ii, jj)]['active'][kk] for ii, jj, kk in ids])
                unk = np.array([W[(ii, jj)]['unknown'][kk] for ii, jj, kk in ids])
                if not act.any():
                    continue
                Jv = np.array([W[(ii, jj)]['J'][kk] for ii, jj, kk in ids])
                Kv = np.array([W[(ii, jj)]['K'][kk] for ii, jj, kk in ids])
                zex = np.array([W[(ii, jj)]['z_out'][kk] for ii, jj, kk in ids])
                dal = 0.5 * (a['dalpha'][k] + a2['dalpha'][k])
                area = geom.ds * dal
                good = act & ~unk & np.isfinite(Jv)
                base = dict(fam=fid, i=i, k=k, nodes=ids, s_mid=s_mid)
                if unk.any() or good.sum() == 0:
                    wgt = area * (np.mean(Kv[good]) if good.any() else 0.0)
                    cells.append(dict(base, kind='unknown', weight=wgt, Jmin=np.nan, Jmax=np.nan))
                    unknown_cells.append(len(cells) - 1); continue
                cuts = [e for e, (p, q) in enumerate(EDGES) if good[p] and good[q] and abs(zex[p] - zex[q]) > jump_tol]
                frac_active = good.mean()
                if len(cuts) == 0:
                    cells.append(dict(base, kind='regular', weight=area * frac_active * np.mean(Kv[good]),
                                      Jmin=np.min(Jv[good]), Jmax=np.max(Jv[good])))
                elif len(cuts) == 2:
                    e1, e2 = cuts
                    if (e2 - e1) % 4 == 2:
                        sideA = [(e1 + 1) % 4, e2]; sideB = [(e2 + 1) % 4, e1]
                    else:
                        shared = e2 if (e2 - e1) % 4 == 1 else e1
                        sideA = [shared]; sideB = [x for x in range(4) if x != shared]
                    made = []
                    for side_nodes, tag in ((sideA, 'A'), (sideB, 'B')):
                        gs = [p for p in side_nodes if good[p]]
                        if not gs:
                            made.append(None); continue
                        cells.append(dict(base, kind='sub', nodes=[ids[p] for p in gs], side=tag,
                                          weight=area * (len(side_nodes) / 4.0) * np.mean(Kv[gs]),
                                          Jmin=np.min(Jv[gs]), Jmax=np.max(Jv[gs]), side_nodes=side_nodes))
                        made.append(len(cells) - 1)
                    cA, cB = made
                    if cA is not None and cB is not None:
                        for e in cuts:
                            p, q = EDGES[e]
                            longp, shortp = (p, q) if abs(zex[p]) > abs(zex[q]) else (q, p)
                            ii, jj, kk = ids[longp]
                            bar_z = W[(ii, jj)]['barrier_z'][kk]; bar_H = W[(ii, jj)]['barrier_H'][kk]
                            sf = geom.surfs[ii]; arc = geom.arcs[ii][jj]
                            if np.isfinite(bar_z):
                                pth = (arc['theta'][kk] + sf.iota * bar_z) % TWO_PI
                                pph = (arc['phi'][kk] + sf.nfp * bar_z) % TWO_PI
                            else:
                                pth = pph = np.nan
                            cl = cA if longp in cells[cA]['side_nodes'] else cB
                            cs = cB if cl == cA else cA
                            transitions.append(dict(cell_long=cl, cell_short=cs, J_long=Jv[longp], J_short=Jv[shortp],
                                                    partner=(ii, pth, pph), barrier_H=bar_H, J3=Jv[longp] - Jv[shortp]))
                            if e % 2 == 1: n_T_s += 1
                            else: n_T_u += 1
                else:
                    cells.append(dict(base, kind='unknown', weight=area * frac_active * np.mean(Kv[good]),
                                      Jmin=np.nan, Jmax=np.nan, ncuts=len(cuts)))
                    unknown_cells.append(len(cells) - 1)
    n_cells = len(cells)
    empty = dict(n_cells=n_cells, time=time.time() - t0, n_T_points=n_T_points, n_gamma=n_gamma, n_multi=n_multi,
                 mismatch=mismatch_surfaces, n_sheets=0, n_unknown_cells=len(unknown_cells), n_trans_linked=0,
                 n_T_u=n_T_u, n_T_s=n_T_s, cells=cells, transitions=transitions)
    if n_cells == 0:
        for key in sources:
            empty[key] = dict(Q_lower=0.0, Q_upper=0.0, Q_total=0.0, Q_unknown=0.0)
        return empty

    # cell adjacency across uncut shared edges
    pos_index = {}
    for cid, c in enumerate(cells):
        pos_index.setdefault((c['fam'], c['i'], c['k']), []).append(cid)
    node_sets = [set(c['nodes']) for c in cells]

    def arc_of(fam, i):
        return geom.arcs[i][[jj for (ii, jj) in geom.families[fam] if ii == i][0]]

    rows, cols = [], []
    for cid, c in enumerate(cells):
        if c['kind'] == 'unknown': continue
        fam, i, k = c['fam'], c['i'], c['k']
        arc_i = arc_of(fam, i); nkk = len(arc_i['theta'])
        for di, dk in ((0, 1), (1, 0)):
            k2 = (k + dk) % nkk if arc_i['closed'] else k + dk
            for cid2 in pos_index.get((fam, i + di, k2), []):
                if cells[cid2]['kind'] == 'unknown': continue
                if len(node_sets[cid] & node_sets[cid2]) >= 2:
                    rows.append(cid); cols.append(cid2)
    adj = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_cells, n_cells))
    adj = (adj + adj.T).tocsr()
    n_sheets, sheet = connected_components(adj, directed=False)

    dJ_sheet = np.ones(n_sheets)
    for sh in range(n_sheets):
        ranges = [cells[c]['Jmax'] - cells[c]['Jmin'] for c in np.nonzero(sheet == sh)[0] if cells[c]['kind'] != 'unknown']
        ranges = [r for r in ranges if r > 0]
        dJ_sheet[sh] = (np.median(ranges) / bins_per_cell) if ranges else 1.0

    def bin_of(cid, J): return int(np.floor(J / dJ_sheet[sheet[cid]]))
    bin_ranges = {cid: (bin_of(cid, c['Jmin']), bin_of(cid, c['Jmax'])) for cid, c in enumerate(cells) if c['kind'] != 'unknown'}
    gid = {}

    def node_id(cid, bn):
        key = (cid, bn)
        if key not in gid: gid[key] = len(gid)
        return gid[key]
    for cid, (b0, b1) in bin_ranges.items():
        for bn in range(b0, b1 + 1): node_id(cid, bn)
    g_rows, g_cols = [], []
    r_arr, c_arr = adj.nonzero()
    for cid, cid2 in zip(r_arr, c_arr):
        if cid >= cid2 or cid not in bin_ranges or cid2 not in bin_ranges: continue
        lo = max(bin_ranges[cid][0], bin_ranges[cid2][0]); hi = min(bin_ranges[cid][1], bin_ranges[cid2][1])
        for bn in range(lo, hi + 1):
            g_rows.append(node_id(cid, bn)); g_cols.append(node_id(cid2, bn))
    n_trans_linked = 0
    for tr in transitions:
        cl, cs = tr['cell_long'], tr['cell_short']
        if cl not in bin_ranges or cs not in bin_ranges: continue
        bl = int(np.clip(bin_of(cl, tr['J_long']), *bin_ranges[cl])); bs = int(np.clip(bin_of(cs, tr['J_short']), *bin_ranges[cs]))
        g_rows.append(node_id(cl, bl)); g_cols.append(node_id(cs, bs))
        ii, pth, pph = tr['partner']
        if np.isfinite(pth):
            best = None
            for jj, arc in enumerate(geom.arcs[ii]):
                dth = (arc['theta'] % TWO_PI - pth + np.pi) % TWO_PI - np.pi
                dph = (arc['phi'] % TWO_PI - pph + np.pi) % TWO_PI - np.pi
                dd = np.hypot(dth, dph); kk = int(np.argmin(dd))
                if best is None or dd[kk] < best[0]: best = (dd[kk], jj, kk)
            if best is not None and best[0] < 0.2:
                _, jj, kk = best
                fam = geom.fam_of[(ii, jj)]
                cands = []
                for (fi, ci, ck) in [(fam, ii, kk), (fam, ii, kk - 1), (fam, ii - 1, kk), (fam, ii - 1, kk - 1)]:
                    cands += pos_index.get((fi, ci, ck), [])
                cands = [c for c in cands if c in bin_ranges and (ii, jj, kk) in node_sets[c]]
                for c3 in cands:
                    b3 = int(np.clip(bin_of(c3, tr['J3']), *bin_ranges[c3]))
                    g_rows.append(node_id(cl, bl)); g_cols.append(node_id(c3, b3))
                    n_trans_linked += 1
    # axis rule: innermost ring, same sheet and bin -> linked
    by_key = {}
    for cid, c in enumerate(cells):
        if c['i'] == 0 and cid in bin_ranges:
            for bn in range(bin_ranges[cid][0], bin_ranges[cid][1] + 1):
                by_key.setdefault((sheet[cid], bn), []).append(cid)
    for key, lst in by_key.items():
        for c1, c2 in zip(lst[:-1], lst[1:]):
            g_rows.append(node_id(c1, key[1])); g_cols.append(node_id(c2, key[1]))
    n_g = len(gid)
    seeds = [node_id(cid, bn) for cid, c in enumerate(cells) if c['i'] == n_s - 2 and cid in bin_ranges
             for bn in range(bin_ranges[cid][0], bin_ranges[cid][1] + 1)]

    def reach(extra_rows, extra_cols, extra_n):
        G = coo_matrix((np.ones(len(g_rows) + len(extra_rows)), (g_rows + extra_rows, g_cols + extra_cols)),
                       shape=(n_g + extra_n, n_g + extra_n))
        ncomp, lab = connected_components(G, directed=False)
        seed_labels = set(lab[s] for s in seeds)
        return np.array([lab[x] in seed_labels for x in range(n_g)])

    reachable_lower = reach([], [], 0) if seeds else np.zeros(n_g, bool)
    ex_r, ex_c = [], []
    for wi, ucid in enumerate(unknown_cells):
        c = cells[ucid]; fam, i, k = c['fam'], c['i'], c['k']
        arc_i = arc_of(fam, i); nkk = len(arc_i['theta'])
        for di, dk in ((0, 1), (1, 0), (0, -1), (-1, 0)):
            k2 = (k + dk) % nkk if arc_i['closed'] else k + dk
            for cid2 in pos_index.get((fam, i + di, k2), []):
                if cid2 in bin_ranges:
                    for bn in range(bin_ranges[cid2][0], bin_ranges[cid2][1] + 1):
                        ex_r.append(n_g + wi); ex_c.append(node_id(cid2, bn))
                elif cells[cid2]['kind'] == 'unknown':
                    ex_r.append(n_g + wi); ex_c.append(n_g + unknown_cells.index(cid2))
        if i == n_s - 2 and seeds:
            ex_r.append(n_g + wi); ex_c.append(seeds[0])
    reachable_upper = reach(ex_r, ex_c, len(unknown_cells)) if seeds else np.zeros(n_g, bool)

    # reachable fraction per cell (lower/upper)
    frac_lower = np.zeros(n_cells); frac_upper = np.zeros(n_cells)
    for cid, c in enumerate(cells):
        if cid not in bin_ranges: continue
        b0, b1 = bin_ranges[cid]; dJ = dJ_sheet[sheet[cid]]; Jr = c['Jmax'] - c['Jmin']
        for bn in range(b0, b1 + 1):
            if Jr <= 0:
                part = 1.0
            else:
                lo = max(bn * dJ, c['Jmin']); hi = min((bn + 1) * dJ, c['Jmax']); part = max(hi - lo, 0.0) / Jr
            g = gid[(cid, bn)]
            if reachable_lower[g]: frac_lower[cid] += part
            if reachable_upper[g]: frac_upper[cid] += part
    is_unknown = np.array([c['kind'] == 'unknown' for c in cells])
    weight0 = np.array([c['weight'] for c in cells]); s_mid = np.array([c['s_mid'] for c in cells])
    res = dict(n_cells=n_cells, n_sheets=n_sheets, n_T_u=n_T_u, n_T_s=n_T_s, n_gamma=n_gamma, n_multi=n_multi,
               n_T_points=n_T_points, mismatch=mismatch_surfaces, n_trans=len(transitions), n_trans_linked=n_trans_linked,
               n_unknown_cells=len(unknown_cells), time=time.time() - t0, cells=cells, transitions=transitions,
               frac_lower=frac_lower, frac_upper=frac_upper, is_unknown=is_unknown)
    for key, h in sources.items():
        w = weight0 * h(np.sqrt(s_mid))
        res[key] = dict(Q_lower=float(np.sum(w * frac_lower)), Q_upper=float(np.sum(w * frac_upper) + np.sum(w[is_unknown])),
                        Q_total=float(np.sum(w)), Q_unknown=float(np.sum(w[is_unknown])))
    if verbose:
        parts = []
        for key in sources:
            r = res[key]; T = r['Q_total']
            parts.append(f"{key}: {r['Q_lower']/T if T else 0:.3f}..{r['Q_upper']/T if T else 0:.3f}")
        print(f"  b={b:.4f}: cells={n_cells} sheets={n_sheets} unknown={len(unknown_cells)} T-points={n_T_points} "
              f"Gamma_max={n_gamma} (mismatch on {len(mismatch_surfaces)} surf, {n_multi} multi-jump edges) "
              f"trans_linked={n_trans_linked}  reachable: {' | '.join(parts)}  {res['time']:.1f}s")
    return res


def theta_map(geom, res, fname):
    """Plot reachable fraction (lower) per cell for each arc family as (u, s) images."""
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cells = res['cells']; fl = res['frac_lower']
    fams = sorted(set(c['fam'] for c in cells))
    fig, axes = plt.subplots(1, len(fams), figsize=(4 * len(fams) + 1.5, 4), squeeze=False)
    im = None
    for ax, fid in zip(axes[0], fams):
        fam = geom.families[fid]
        nk = max(len(geom.arcs[i][j]['theta']) for i, j in fam)
        img = np.full((len(geom.s) - 1, nk), np.nan)
        for cid, c in enumerate(cells):
            if c['fam'] != fid: continue
            val = -0.2 if c['kind'] == 'unknown' else fl[cid]
            prev = img[c['i'], c['k']]
            img[c['i'], c['k']] = val if np.isnan(prev) else 0.5 * (prev + val)
        im = ax.imshow(img, origin='lower', aspect='auto', vmin=-0.2, vmax=1, cmap='viridis',
                       extent=[0, 1, geom.s[0], geom.s[-1]])
        a0 = geom.arcs[fam[0][0]][fam[0][1]]
        ax.set_title(f"family {fid}: {'closed ' + str(a0['wind']) if a0['closed'] else 'open arc'}, "
                     f"s in [{geom.s[fam[0][0]]:.2f},{geom.s[fam[-1][0]]:.2f}]", fontsize=8)
        ax.set_xlabel('ridge parameter u'); ax.set_ylabel('s')
    if im is not None:
        fig.colorbar(im, ax=axes[0].tolist(), label='reachable fraction (lower); -0.2 = unknown')
    fig.suptitle(os.path.basename(fname), fontsize=9)
    fig.savefig(fname, dpi=110); plt.close(fig)


def run(name, n_s, n_closed, n_lam, s_min=0.02, outdir=None, levels=(0.05, 0.1, 0.5, 0.8, 0.9, 0.95), plots=True):
    geom = Geometry(name, n_s=n_s, s_min=s_min, n_closed=n_closed, n_open=max(16, n_closed // 5))
    Vh = {key: geom.denominator(h) for key, h in SOURCES.items()}
    print('  V_h =', {k: round(v, 5) for k, v in Vh.items()})
    for lam in levels:
        b = geom.Bmin_g + lam * (geom.Bmax_g - geom.Bmin_g)
        res = slice_at_pitch(geom, b, verbose=True)
        led = {key: geom.trapped_weight_volume(b, h) for key, h in SOURCES.items()}
        ratios = {key: (res[key]['Q_total'] / led[key] if led[key] else float('nan')) for key in SOURCES}
        print(f'      ledger ratio (cells / volume integral): ' + ', '.join(f'{k}: {v:.4f}' for k, v in ratios.items())
              + (f'; mismatch (s, T, Gamma): {res["mismatch"][:5]}' if res['mismatch'] else ''))
        if plots and outdir and res['n_cells']:
            theta_map(geom, res, f'{outdir}/{name}_ns{n_s}_nu{n_closed}_smin{s_min}_lam{lam}.png')
    out = {}
    if n_lam:
        t0 = time.time()
        edges = np.linspace(geom.Bmin_g, geom.Bmax_g, n_lam + 1); bm = 0.5 * (edges[1:] + edges[:-1]); db = edges[1] - edges[0]
        acc = {key: dict(QL=0.0, QU=0.0, QT=0.0, QUnk=0.0) for key in SOURCES}
        n_mismatch = 0
        for b in bm:
            r = slice_at_pitch(geom, b)
            n_mismatch += len(r['mismatch'])
            for key in SOURCES:
                acc[key]['QL'] += r[key]['Q_lower'] * db / b ** 2; acc[key]['QU'] += r[key]['Q_upper'] * db / b ** 2
                acc[key]['QT'] += r[key]['Q_total'] * db / b ** 2; acc[key]['QUnk'] += r[key]['Q_unknown'] * db / b ** 2
        dt = time.time() - t0
        for key in SOURCES:
            a = acc[key]; D = 2 * Vh[key]
            out[key] = dict(f_lower=a['QL'] / D, f_upper=a['QU'] / D, trapped=a['QT'] / D, unknown=a['QUnk'] / D)
            print(f'  [{key}] f_lower={a["QL"]/D:.4f} f_upper={a["QU"]/D:.4f} trapped fraction={a["QT"]/D:.4f} unknown weight={a["QUnk"]/D:.4f}')
        print(f'  ({n_lam} slices in {dt:.1f}s, {dt/n_lam:.2f}s/slice; geometry {geom.t_geom:.1f}s; {n_mismatch} surface mismatches over all slices)')
    return geom, out


if __name__ == '__main__':
    name = sys.argv[1] if len(sys.argv) > 1 else 'DMerc'
    n_s = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    n_closed = int(sys.argv[3]) if len(sys.argv) > 3 else 128
    n_lam = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    s_min = float(sys.argv[5]) if len(sys.argv) > 5 else 0.02
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(outdir, exist_ok=True)
    run(name, n_s, n_closed, n_lam, s_min=s_min, outdir=outdir)
