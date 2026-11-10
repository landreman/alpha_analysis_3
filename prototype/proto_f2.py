"""Prototype v2: ridge-based algorithm with node-level radial matching.

Differences from proto_f.py:
* nodes on open arcs are placed with the same linear density as on closed curves;
* cells between adjacent surfaces are built from each ridge node and its nearest ridge node on
  the next surface (any arc), so merging/splitting arcs are handled; nodes without a partner,
  and partner nodes not covered by any cell, become explicit unknown weight;
* cells may be quads or triangles; adjacency is 'two shared nodes';
* seeds are cells touching s = 1; the axis rule links same-bin cells touching s_min;
* Theta maps are scatter plots in the physical (alpha, s) plane.
"""
from __future__ import annotations
import sys, time, os
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

MAX_BINS_PER_CELL = int(os.environ.get('MAX_BINS_PER_CELL', '128'))   # dyadic coarsening above this
MEM_CAP_GB = float(os.environ.get('MEM_CAP_GB', '3'))                 # abort a run instead of swapping


def _rss_gb():
    try:
        return int(os.popen(f'ps -o rss= -p {os.getpid()}').read()) / 1048576
    except Exception:
        return 0.0

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from proto_f import (Surface, ridge_arcs, chains, wells_at_pitch, SOURCES, FILES, ROOT, TWO_PI, trapz,
                     BoozerField)


def place_nodes(arcs, density):
    """Resample arcs to nodes with `density` nodes per radian (angle for closed curves, arc length for open)."""
    out = []
    for a in arcs:
        th, ph = a['theta'], a['phi']
        if a['closed']:
            p, q = a['wind']
            seg = np.hypot(np.diff(th), np.diff(ph)); s_arc = np.r_[0, np.cumsum(seg)]
            tot = s_arc[-1] + np.hypot(th[0] + TWO_PI * p - th[-1], ph[0] + TWO_PI * q - ph[-1])
            n = max(16, int(round(density * TWO_PI)))
            uu = np.linspace(0, tot, n, endpoint=False)
            out.append(dict(theta=np.interp(uu, np.r_[s_arc, tot], np.r_[th, th[0] + TWO_PI * p]),
                            phi=np.interp(uu, np.r_[s_arc, tot], np.r_[ph, ph[0] + TWO_PI * q]),
                            closed=True, wind=(p, q)))
        else:
            seg = np.hypot(np.diff(th), np.diff(ph)); s_arc = np.r_[0, np.cumsum(seg)]
            if s_arc[-1] < 1e-3: continue
            n = max(6, int(round(density * s_arc[-1])) + 1)
            uu = np.linspace(0, s_arc[-1], n)
            out.append(dict(theta=np.interp(uu, s_arc, th), phi=np.interp(uu, s_arc, ph), closed=False, wind=(0, 0)))
    return out


class Geometry:
    """Pitch-independent data: surfaces, arcs with nodes and chains, nearest-node partners between surfaces."""

    def __init__(self, name, s_grid, density=128 / TWO_PI, periods=48, NT=512, NZ=512, partner_tol=0.4, verbose=True):
        t0 = time.time()
        self.name = name
        self.fld = BoozerField(ROOT + FILES[name])
        self.nfp = self.fld.nfp; self.L = TWO_PI / self.nfp
        self.s = np.asarray(s_grid, float)
        self.surfs = []; self.arcs = []
        self.Bmin_g = np.inf; self.Bmax_g = -np.inf
        n_nodes = 0
        for s in self.s:
            surf = Surface(self.fld, s, NT, NZ)
            arcs = place_nodes(ridge_arcs(surf), density)
            surf.DB = surf.D2B = surf._DBc = surf._D2Bc = None   # derivative tables only serve ridge extraction
            for a in arcs:
                a['H'] = surf.b(a['theta'], a['phi'])
                a['alpha'] = a['theta'] - surf.iota * a['phi'] / self.nfp
                dal = np.diff(a['alpha'])
                if a['closed']:
                    p, q = a['wind']
                    total = TWO_PI * p - surf.iota * self.L * q
                    dal = np.r_[dal, a['alpha'][0] + total - a['alpha'][-1]]
                    a['alpha_total'] = total
                a['dalpha'] = np.abs(dal)
                a['ch'] = chains(surf, a['theta'], a['phi'], periods=periods)
                n_nodes += len(a['theta'])
            self.surfs.append(surf); self.arcs.append(arcs)
            self.Bmin_g = min(self.Bmin_g, surf.Bmin); self.Bmax_g = max(self.Bmax_g, surf.Bmax)
        self.n_nodes = n_nodes
        # flat node tables per surface and nearest partners on the next surface
        self.flat = []
        for i in range(len(self.s)):
            labels = []; th = []; ph = []
            for j, a in enumerate(self.arcs[i]):
                for k in range(len(a['theta'])):
                    labels.append((j, k))
                th.append(a['theta'] % TWO_PI); ph.append(a['phi'] % TWO_PI)
            self.flat.append(dict(labels=labels, theta=np.concatenate(th) if th else np.zeros(0),
                                  phi=np.concatenate(ph) if ph else np.zeros(0)))
        self.partner = []      # partner[i][(j,k)] -> (j',k') on surface i+1 or None
        self.n_unpartnered = 0
        for i in range(len(self.s) - 1):
            A, B = self.flat[i], self.flat[i + 1]
            part = {}
            if len(A['labels']) and len(B['labels']):
                dth = (A['theta'][:, None] - B['theta'][None, :] + np.pi) % TWO_PI - np.pi
                dph = (A['phi'][:, None] - B['phi'][None, :] + np.pi) % TWO_PI - np.pi
                D = np.hypot(dth, dph)
                nn = np.argmin(D, axis=1); dmin = D[np.arange(len(nn)), nn]
                for idx, lab in enumerate(A['labels']):
                    part[lab] = B['labels'][nn[idx]] if dmin[idx] < partner_tol else None
                    if part[lab] is None: self.n_unpartnered += 1
            else:
                for lab in A['labels']: part[lab] = None
            self.partner.append(part)
        self.t_geom = time.time() - t0
        if verbose:
            n_closed = sum(1 for arcs in self.arcs for a in arcs if a['closed'])
            n_open = sum(1 for arcs in self.arcs for a in arcs if not a['closed'])
            print(f'[{name}] geometry: {len(self.s)} surfaces (s in [{self.s[0]},{self.s[-1]}]), {n_nodes} ridge nodes, '
                  f'{n_closed} closed + {n_open} open arcs, {self.n_unpartnered} unpartnered nodes, '
                  f'B in [{self.Bmin_g:.4f},{self.Bmax_g:.4f}], {self.t_geom:.1f} s')

    def denominator(self, h):
        vals = np.array([h(np.sqrt(sf.s)) * abs(sf.C) * sf.inv_B2_mean * TWO_PI * self.L for sf in self.surfs])
        return trapz(vals, self.s)

    def trapped_weight_volume(self, b, h):
        vals = []
        for sf in self.surfs:
            if sf.Bmax <= b:
                vals.append(0.0); continue
            g = 1.0 - sf.B / b; m = g > 0
            with np.errstate(divide='ignore', invalid='ignore'):
                integrand = np.where(m, 1.0 / (sf.B * np.sqrt(np.where(m, g, 1.0))), 0.0)
            vals.append(h(np.sqrt(sf.s)) * abs(sf.C) * integrand.mean() * TWO_PI * self.L)
        return trapz(np.array(vals), self.s)


EDGE_RANGE_LINKS = os.environ.get('EDGE_RANGE_LINKS', '1') == '1'


def _dalpha_between(arc, k0, k1):
    """|alpha(k1) - alpha(k0)| along an arc, accounting for the seam of closed curves."""
    d = abs(arc['alpha'][k1] - arc['alpha'][k0])
    if arc['closed']:
        t = abs(arc['alpha_total'])
        d = min(d, abs(d - t))
    return d


def slice_at_pitch(geom: Geometry, b, jump_frac=0.5, bins_per_cell=1.0, verbose=False, sources=SOURCES):
    """jump_frac: an edge is cut when the exit maximum moves by more than this fraction of the
    spacing between neighbouring maxima along the line, i.e. the exit is a different maximum."""
    t0 = time.time()
    n_s = len(geom.s)
    W = {(i, j): wells_at_pitch(geom.surfs[i], a['theta'], a['phi'], a['H'], a['ch'], b)
         for i in range(n_s) for j, a in enumerate(geom.arcs[i])}

    def good_of(n):
        i, j, k = n; w = W[(i, j)]
        return bool(w['active'][k] and not w['unknown'][k] and np.isfinite(w['J'][k]))

    def val(n, key):
        i, j, k = n; return W[(i, j)][key][k]

    def is_transition(na, nb, frac=0.3):
        """True when the exit maximum of the shorter well is an interior maximum of the longer well,
        i.e. an interior barrier crossed b between the two states."""
        za, zb = val(na, 'z_out'), val(nb, 'z_out')
        if not (np.isfinite(za) and np.isfinite(zb)): return False
        long_, short = (na, nb) if abs(za) > abs(zb) else (nb, na)
        i, j, k = long_; ne = W[(i, j)]['n_exit'][k]
        if ne <= 0: return False
        interior = geom.arcs[i][j]['ch']['max_z'][k, :ne]
        zs = val(short, 'z_exit'); gs = val(short, 'exit_gap')
        tol = frac * (gs if np.isfinite(gs) else abs(zs))
        return bool(np.any(np.abs(interior - zs) < tol))

    # per-surface transition points vs Gamma_max crossings
    per_surface = []; n_multi = 0
    for i in range(n_s):
        nT = nG = 0
        for j, a in enumerate(geom.arcs[i]):
            w = W[(i, j)]
            good = w['active'] & ~w['unknown'] & np.isfinite(w['z_exit'])
            ne = w['n_exit']; hb = a['H'] - b; nk = len(a['theta'])
            pairs = [(k, (k + 1) % nk) for k in range(nk if a['closed'] else nk - 1)]
            if a['closed']:
                nG += int(np.sum(hb * np.roll(hb, -1) < 0))
            else:
                nG += int(np.sum(hb[:-1] * hb[1:] < 0))
            for k, k1 in pairs:
                if good[k] and good[k1]:
                    if is_transition((i, j, k), (i, j, k1)): nT += 1
                    if abs(ne[k] - ne[k1]) >= 2: n_multi += 1
        per_surface.append((nT, nG))
    n_T_points = sum(t for t, g in per_surface); n_gamma = sum(g for t, g in per_surface)
    mismatch = [(round(float(geom.s[i]), 3), t, g) for i, (t, g) in enumerate(per_surface) if t != g]

    # cells
    cells = []; transitions = []; unknown_cells = []
    n_frontier = 0
    for i in range(n_s - 1):
        ds = geom.s[i + 1] - geom.s[i]; s_mid = 0.5 * (geom.s[i] + geom.s[i + 1])
        part = geom.partner[i]
        used_top = set()
        for j, a in enumerate(geom.arcs[i]):
            nk = len(a['theta']); nK = nk if a['closed'] else nk - 1
            w = W[(i, j)]
            for k in range(nK):
                k1 = (k + 1) % nk
                n0, n1 = (i, j, k), (i, j, k1)
                act0, act1 = w['active'][k], w['active'][k1]
                p0, p1 = part[(j, k)], part[(j, k1)]
                top_active = False
                if p0 is not None: top_active |= bool(W[(i + 1, p0[0])]['active'][p0[1]])
                if p1 is not None: top_active |= bool(W[(i + 1, p1[0])]['active'][p1[1]])
                if not (act0 or act1 or top_active):
                    continue
                dal_b = a['dalpha'][k]
                if p0 is None or p1 is None:
                    gb = [n for n in (n0, n1) if good_of(n)]
                    Km = np.mean([val(n, 'K') for n in gb]) if gb else 0.0
                    cells.append(dict(kind='unknown', nodes=[n0, n1], i=i, s_mid=s_mid, weight=ds * dal_b * Km,
                                      Jmin=np.nan, Jmax=np.nan, why='frontier'))
                    unknown_cells.append(len(cells) - 1); n_frontier += 1
                    continue
                t0n, t1n = (i + 1, p0[0], p0[1]), (i + 1, p1[0], p1[1])
                used_top.add(p0); used_top.add(p1)
                nodes = [n0, n1]
                if t1n != n1 and t1n not in nodes: nodes.append(t1n)
                if t0n not in nodes: nodes.append(t0n)
                if p0[0] == p1[0]:
                    arc_t = geom.arcs[i + 1][p0[0]]
                    dal_t = _dalpha_between(arc_t, p0[1], p1[1]) if p0 != p1 else 0.0
                    nkt = len(arc_t['theta'])
                    dk = abs(p0[1] - p1[1]); dk = min(dk, nkt - dk) if arc_t['closed'] else dk
                    if dk > 4:
                        dal_t = dal_b   # suspicious long top edge; treat as unknown below
                else:
                    dal_t = dal_b; dk = 99
                area = ds * 0.5 * (dal_b + dal_t)
                act = np.array([val(n, 'active') for n in nodes]); unk = np.array([val(n, 'unknown') for n in nodes])
                good = np.array([good_of(n) for n in nodes])
                Jv = np.array([val(n, 'J') for n in nodes]); Kv = np.array([val(n, 'K') for n in nodes])
                zex = np.array([val(n, 'z_exit') for n in nodes]); gaps = np.array([val(n, 'exit_gap') for n in nodes])
                zout = np.array([val(n, 'z_out') for n in nodes])
                base = dict(nodes=nodes, i=i, s_mid=s_mid)
                if unk.any() or good.sum() == 0 or dk > 4:
                    wgt = area * (np.mean(Kv[good]) if good.any() else 0.0)
                    cells.append(dict(base, kind='unknown', weight=wgt, Jmin=np.nan, Jmax=np.nan,
                                      why='unknown-node' if unk.any() else ('top-span' if dk > 4 else 'no-good')))
                    unknown_cells.append(len(cells) - 1); continue
                m = len(nodes)
                edges = [(e, (e + 1) % m) for e in range(m)]
                cuts = [e for e, (p, q) in enumerate(edges) if good[p] and good[q] and is_transition(nodes[p], nodes[q])]
                frac_active = good.mean()
                if len(cuts) == 0:
                    cells.append(dict(base, kind='regular', weight=area * frac_active * np.mean(Kv[good]),
                                      Jmin=np.min(Jv[good]), Jmax=np.max(Jv[good])))
                elif len(cuts) == 2:
                    e1, e2 = cuts
                    # side A: nodes strictly after edge e1 up to and including the start of e2 (cyclic)
                    sideA = []
                    x = (e1 + 1) % m
                    while True:
                        sideA.append(x)
                        if x == e2: break
                        x = (x + 1) % m
                    sideB = [x for x in range(m) if x not in sideA]
                    made = []
                    for side_nodes in (sideA, sideB):
                        gs = [p for p in side_nodes if good[p]]
                        if not gs:
                            made.append(None); continue
                        cells.append(dict(base, kind='sub', nodes=[nodes[p] for p in gs], side_idx=side_nodes,
                                          weight=area * (len(side_nodes) / m) * np.mean(Kv[gs]),
                                          Jmin=np.min(Jv[gs]), Jmax=np.max(Jv[gs])))
                        made.append(len(cells) - 1)
                    cA, cB = made
                    if cA is not None and cB is not None:
                        ends = []
                        for e in cuts:
                            p, q = edges[e]
                            longp, shortp = (p, q) if abs(zout[p]) > abs(zout[q]) else (q, p)
                            ii, jj, kk = nodes[longp]
                            bar_z = W[(ii, jj)]['barrier_z'][kk]
                            sf = geom.surfs[ii]; arc = geom.arcs[ii][jj]
                            if np.isfinite(bar_z):
                                pth = (arc['theta'][kk] + sf.iota * bar_z) % TWO_PI
                                pph = (arc['phi'][kk] + sf.nfp * bar_z) % TWO_PI
                            else:
                                pth = pph = np.nan
                            cl = cA if longp in cells[cA]['side_idx'] else cB
                            ends.append(dict(cl=cl, J_long=Jv[longp], J_short=Jv[shortp], partner=(ii, pth, pph)))
                        if ends[0]['cl'] == ends[1]['cl']:
                            cl = ends[0]['cl']; cs = cB if cl == cA else cA
                            transitions.append(dict(cell_long=cl, cell_short=cs, ends=ends))   # one T segment
                        else:
                            for en in ends:   # two different barriers cut this cell: two point transitions
                                cl = en['cl']; cs = cB if cl == cA else cA
                                transitions.append(dict(cell_long=cl, cell_short=cs, ends=[en]))
                else:
                    cells.append(dict(base, kind='unknown', weight=area * frac_active * np.mean(Kv[good]),
                                      Jmin=np.nan, Jmax=np.nan, why=f'{len(cuts)}-cuts'))
                    unknown_cells.append(len(cells) - 1)
        # top nodes with active wells that no cell from below touched: unknown weight
        for jt, at in enumerate(geom.arcs[i + 1]):
            wt = W[(i + 1, jt)]
            for kt in range(len(at['theta'])):
                if (jt, kt) in used_top or not wt['active'][kt]: continue
                n = (i + 1, jt, kt)
                nd = len(at['dalpha'])
                e1 = min(kt, nd - 1); e0 = (kt - 1) % nd if at['closed'] else max(kt - 1, 0)
                dal = 0.5 * (at['dalpha'][e1] + at['dalpha'][e0])
                Km = wt['K'][kt] if np.isfinite(wt['K'][kt]) else 0.0
                cells.append(dict(kind='unknown', nodes=[n], i=i, s_mid=s_mid, weight=0.5 * ds * dal * Km,
                                  Jmin=np.nan, Jmax=np.nan, why='uncovered-top'))
                unknown_cells.append(len(cells) - 1)

    n_cells = len(cells)
    res = dict(n_cells=n_cells, n_T_points=n_T_points, n_gamma=n_gamma, n_multi=n_multi, mismatch=mismatch,
               n_unknown_cells=len(unknown_cells), n_frontier=n_frontier, n_trans=len(transitions),
               cells=cells, transitions=transitions, n_sheets=0, n_trans_linked=0)
    if n_cells == 0:
        for key in sources: res[key] = dict(Q_lower=0.0, Q_upper=0.0, Q_total=0.0, Q_unknown=0.0)
        res['time'] = time.time() - t0
        return res

    node_to_cells = {}
    for cid, c in enumerate(cells):
        for n in c['nodes']:
            node_to_cells.setdefault(n, []).append(cid)
    # adjacency: two non-unknown cells sharing two nodes
    rows, cols = [], []
    edge_J = {}   # (cid, cid2) -> J range of the shared edge(s): contours cross only at those levels
    for cid, c in enumerate(cells):
        if c['kind'] == 'unknown': continue
        ns = c['nodes']; m = len(ns)
        for e in range(m):
            na, nb = ns[e], ns[(e + 1) % m]
            if na == nb: continue
            for cid2 in node_to_cells.get(na, []):
                if cid2 <= cid or cells[cid2]['kind'] == 'unknown': continue
                if nb in cells[cid2]['nodes']:
                    fin = [x for x in (val(na, 'J'), val(nb, 'J')) if np.isfinite(x)]
                    if not fin: continue          # edge at the Gamma boundary with no action data: no link
                    rows.append(cid); cols.append(cid2)
                    lo, hi = edge_J.get((cid, cid2), (np.inf, -np.inf))
                    edge_J[(cid, cid2)] = (min([lo] + fin), max([hi] + fin))
    adj = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n_cells, n_cells))
    adj = (adj + adj.T).tocsr()
    n_sheets, sheet = connected_components(adj, directed=False)
    dJ_sheet = np.ones(n_sheets)
    for sh in range(n_sheets):
        ranges = [cells[c]['Jmax'] - cells[c]['Jmin'] for c in np.nonzero(sheet == sh)[0] if cells[c]['kind'] != 'unknown']
        ranges = [r for r in ranges if r > 0]
        dJ_sheet[sh] = (np.median(ranges) / bins_per_cell) if ranges else 1.0

    # Dyadic multi-resolution bins: a cell whose J range would span more than MAX_BINS_PER_CELL sheet bins is
    # binned at the coarser width 2**k * dJ_sheet.  Coarse bins are exact unions of fine bins, so a link between
    # cells of different levels mixes contours only inside the coarse bin; the J interval never diffuses along
    # a path of cells.  Without the cap one outlier cell can create 1e5+ graph nodes and exhaust memory.
    lev = np.zeros(n_cells, int)
    n_bins_uncapped = 0
    for cid, c in enumerate(cells):
        if c['kind'] == 'unknown' or not (np.isfinite(c['Jmin']) and np.isfinite(c['Jmax'])): continue
        dJ0 = dJ_sheet[sheet[cid]]
        def nbins(k): return int(np.floor(c['Jmax'] / (dJ0 * 2.0 ** k))) - int(np.floor(c['Jmin'] / (dJ0 * 2.0 ** k))) + 1
        n_bins_uncapped += nbins(0)
        k = 0
        while nbins(k) > MAX_BINS_PER_CELL and k < 60: k += 1
        lev[cid] = k
    n_capped = int(np.sum(lev > 0))

    def dJ_of(cid): return dJ_sheet[sheet[cid]] * 2.0 ** lev[cid]

    def bin_of(cid, J): return int(np.floor(J / dJ_of(cid)))
    bin_ranges = {cid: (bin_of(cid, c['Jmin']), bin_of(cid, c['Jmax'])) for cid, c in enumerate(cells) if c['kind'] != 'unknown'}
    gid = {}

    def node_id(cid, bn):
        key = (cid, bn)
        if key not in gid: gid[key] = len(gid)
        return gid[key]
    for cid, (b0, b1) in bin_ranges.items():
        for bn in range(b0, b1 + 1): node_id(cid, bn)
    g_rows, g_cols = [], []

    def link_cells(ca, cb, J0, J1):
        """Link the bins of cells ca and cb over the J interval [J0, J1]; fine bins map onto the coarser cell's bins."""
        fine, coarse = (ca, cb) if lev[ca] <= lev[cb] else (cb, ca)
        d = int(lev[coarse] - lev[fine])
        lo = max(bin_ranges[fine][0], bin_ranges[coarse][0] << d)
        hi = min(bin_ranges[fine][1], ((bin_ranges[coarse][1] + 1) << d) - 1)
        if np.isfinite(J0): lo = max(lo, bin_of(fine, J0))
        if np.isfinite(J1): hi = min(hi, bin_of(fine, J1))
        for bn in range(lo, hi + 1):
            g_rows.append(node_id(fine, bn)); g_cols.append(node_id(coarse, bn >> d))
    r_arr, c_arr = adj.nonzero()
    for cid, cid2 in zip(r_arr, c_arr):
        if cid >= cid2 or cid not in bin_ranges or cid2 not in bin_ranges: continue
        eJ = edge_J.get((cid, cid2))
        if eJ is None: continue
        if EDGE_RANGE_LINKS:
            link_cells(cid, cid2, eJ[0], eJ[1])
        else:
            link_cells(cid, cid2, -np.inf, np.inf)
    # transition links
    n_trans_linked = 0

    def partner_cells(ii, pth, pph):
        if not np.isfinite(pth): return []
        F = geom.flat[ii]
        if not len(F['labels']): return []
        dth = (F['theta'] - pth + np.pi) % TWO_PI - np.pi; dph = (F['phi'] - pph + np.pi) % TWO_PI - np.pi
        dd = np.hypot(dth, dph); idx = int(np.argmin(dd))
        if dd[idx] >= 0.2: return []
        jj, kk = F['labels'][idx]
        return [c3 for c3 in node_to_cells.get((ii, jj, kk), []) if c3 in bin_ranges]

    for tr in transitions:
        cl, cs = tr['cell_long'], tr['cell_short']
        if cl not in bin_ranges or cs not in bin_ranges: continue
        ends = tr['ends']
        Jl = np.array([e['J_long'] for e in ends]); Js = np.array([e['J_short'] for e in ends])
        c3s = set()
        for e in ends: c3s.update(partner_cells(*e['partner']))
        nb = abs(bin_of(cl, Jl.max()) - bin_of(cl, Jl.min())) + abs(bin_of(cs, Js.max()) - bin_of(cs, Js.min()))
        ts = np.linspace(0.0, 1.0, 2 * nb + 2) if len(ends) == 2 else np.array([0.0])
        for t in ts:
            jl = Jl[0] + t * (Jl[-1] - Jl[0]); js = Js[0] + t * (Js[-1] - Js[0]); j3 = jl - js
            bl = int(np.clip(bin_of(cl, jl), *bin_ranges[cl])); bs = int(np.clip(bin_of(cs, js), *bin_ranges[cs]))
            g_rows.append(node_id(cl, bl)); g_cols.append(node_id(cs, bs))
            for c3 in c3s:
                b3 = int(np.clip(bin_of(c3, j3), *bin_ranges[c3]))
                g_rows.append(node_id(cl, bl)); g_cols.append(node_id(c3, b3))
        n_trans_linked += len(c3s)
    # axis rule: a contour entering the omitted disk s < s_min leaves it where J takes the same value.
    # Link inner-ring cells of the same sheet where the J ranges of their ring edges intersect.
    ring = {}
    for cid, c in enumerate(cells):
        if cid not in bin_ranges: continue
        Jr = [val(n, 'J') for n in c['nodes'] if n[0] == 0 and np.isfinite(val(n, 'J'))]
        if Jr:
            ring.setdefault(sheet[cid], []).append((cid, min(Jr), max(Jr)))
    for sh, lst in ring.items():
        lst.sort(key=lambda t: t[1])
        for a_i in range(len(lst)):
            c1, lo1, hi1 = lst[a_i]
            for b_i in range(a_i + 1, len(lst)):
                c2, lo2, hi2 = lst[b_i]
                if lo2 > hi1: break
                link_cells(c1, c2, max(lo1, lo2), min(hi1, hi2))
    n_g = len(gid)
    seeds = [node_id(cid, bn) for cid, c in enumerate(cells) if cid in bin_ranges and any(n[0] == n_s - 1 for n in c['nodes'])
             for bn in range(bin_ranges[cid][0], bin_ranges[cid][1] + 1)]

    def reach(extra_rows, extra_cols, extra_n):
        G = coo_matrix((np.ones(len(g_rows) + len(extra_rows)), (g_rows + extra_rows, g_cols + extra_cols)),
                       shape=(n_g + extra_n, n_g + extra_n))
        ncomp, lab = connected_components(G, directed=False)
        seed_labels = set(lab[s] for s in seeds)
        return np.array([lab[x] in seed_labels for x in range(n_g)])

    reachable_lower = reach([], [], 0) if seeds else np.zeros(n_g, bool)
    ex_r, ex_c = [], []
    uidx = {ucid: wi for wi, ucid in enumerate(unknown_cells)}
    for ucid, wi in uidx.items():
        c = cells[ucid]
        for n in c['nodes']:
            for cid2 in node_to_cells.get(n, []):
                if cid2 == ucid: continue
                if cid2 in bin_ranges:
                    for bn in range(bin_ranges[cid2][0], bin_ranges[cid2][1] + 1):
                        ex_r.append(n_g + wi); ex_c.append(node_id(cid2, bn))
                elif cid2 in uidx:
                    ex_r.append(n_g + wi); ex_c.append(n_g + uidx[cid2])
            if n[0] == n_s - 1 and seeds:
                ex_r.append(n_g + wi); ex_c.append(seeds[0])
    reachable_upper = reach(ex_r, ex_c, len(unknown_cells)) if seeds else np.zeros(n_g, bool)

    frac_lower = np.zeros(n_cells); frac_upper = np.zeros(n_cells)
    for cid, c in enumerate(cells):
        if cid not in bin_ranges: continue
        b0, b1 = bin_ranges[cid]; dJ = dJ_of(cid); Jr = c['Jmax'] - c['Jmin']
        for bn in range(b0, b1 + 1):
            part = 1.0 if Jr <= 0 else max(min((bn + 1) * dJ, c['Jmax']) - max(bn * dJ, c['Jmin']), 0.0) / Jr
            g = gid[(cid, bn)]
            if reachable_lower[g]: frac_lower[cid] += part
            if reachable_upper[g]: frac_upper[cid] += part
    is_unknown = np.array([c['kind'] == 'unknown' for c in cells])
    weight0 = np.array([c['weight'] for c in cells]); s_mid = np.array([c['s_mid'] for c in cells])
    res.update(n_g=n_g, n_capped=n_capped, max_level=int(lev.max()) if n_cells else 0, n_bins_uncapped=n_bins_uncapped)
    if verbose and n_capped:
        top = sorted(bin_ranges, key=lambda cid: -lev[cid])[:3]
        print(f'      graph: {n_g} nodes ({n_bins_uncapped} without the cap), {n_capped} capped cells, max level {lev.max()}; worst: '
              + '; '.join(f"s={cells[cid]['s_mid']:.3f} {cells[cid]['kind']} J=[{cells[cid]['Jmin']:.3g},{cells[cid]['Jmax']:.3g}] "
                          f"sheet dJ={dJ_sheet[sheet[cid]]:.2g} level={lev[cid]} n_exit={[val(n, 'n_exit') for n in cells[cid]['nodes']]}" for cid in top))
    res.update(n_sheets=n_sheets, n_trans_linked=n_trans_linked, frac_lower=frac_lower, frac_upper=frac_upper,
               is_unknown=is_unknown, time=time.time() - t0)
    why = {}
    for ucid in unknown_cells:
        why[cells[ucid].get('why', '?')] = why.get(cells[ucid].get('why', '?'), 0) + 1
    res['unknown_reasons'] = why
    for key, h in sources.items():
        w = weight0 * h(np.sqrt(s_mid))
        res[key] = dict(Q_lower=float(np.sum(w * frac_lower)), Q_upper=float(np.sum(w * frac_upper) + np.sum(w[is_unknown])),
                        Q_total=float(np.sum(w)), Q_unknown=float(np.sum(w[is_unknown])))
    if verbose:
        parts = []
        for key in sources:
            r = res[key]; T = r['Q_total']
            parts.append(f"{key}: {r['Q_lower']/T if T else 0:.3f}..{r['Q_upper']/T if T else 0:.3f}")
        print(f"  b={b:.4f}: cells={n_cells} sheets={n_sheets} unknown={len(unknown_cells)} {why} T-points={n_T_points} "
              f"Gamma_max={n_gamma} (mismatch {len(mismatch)} surf, {n_multi} multi-jump edges) "
              f"trans={len(transitions)}/{n_trans_linked}  reachable: {' | '.join(parts)}  {res['time']:.1f}s")
    return res


def theta_map(geom, res, fname, title=''):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    cells = res['cells']
    acc = {}
    for cid, c in enumerate(cells):
        v = -0.2 if c['kind'] == 'unknown' else res['frac_lower'][cid]
        vu = -0.2 if c['kind'] == 'unknown' else res['frac_upper'][cid]
        for n in c['nodes']:
            acc.setdefault(n, []).append((v, vu))
    xs, ys, vl, vu = [], [], [], []
    for n, lst in acc.items():
        i, j, k = n
        xs.append(geom.arcs[i][j]['alpha'][k] % TWO_PI); ys.append(geom.s[i])
        vl.append(np.mean([a for a, _ in lst])); vu.append(np.mean([b for _, b in lst]))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, v, lab in zip(axes, (vl, vu), ('lower', 'upper')):
        sc = ax.scatter(xs, ys, c=v, s=6, vmin=-0.2, vmax=1, cmap='viridis', linewidths=0)
        ax.set_xlabel('field-line label alpha of the state (mod 2pi)'); ax.set_title(f'reachable fraction ({lab}); -0.2 = unknown', fontsize=9)
        ax.set_xlim(0, TWO_PI); ax.set_ylim(0, 1)
    axes[0].set_ylabel('s')
    fig.colorbar(sc, ax=axes.tolist())
    fig.suptitle(title or os.path.basename(fname), fontsize=9)
    fig.savefig(fname, dpi=110); plt.close(fig)


def run(name, s_grid, density, n_lam, outdir=None, levels=(0.05, 0.1, 0.5, 0.8, 0.9, 0.95), plots=True, tag='',
        bins_per_cell=1.0):
    geom = Geometry(name, s_grid, density=density)
    Vh = {key: geom.denominator(h) for key, h in SOURCES.items()}
    for lam in levels:
        b = geom.Bmin_g + lam * (geom.Bmax_g - geom.Bmin_g)
        res = slice_at_pitch(geom, b, verbose=True, bins_per_cell=bins_per_cell)
        if MEM_CAP_GB and _rss_gb() > MEM_CAP_GB:
            print(f'  ABORT: resident memory {_rss_gb():.2f} GB exceeds MEM_CAP_GB={MEM_CAP_GB}', flush=True); sys.exit(2)
        led = {key: geom.trapped_weight_volume(b, h) for key, h in SOURCES.items()}
        ratios = {key: (res[key]['Q_total'] / led[key] if led[key] else float('nan')) for key in SOURCES}
        unk = {key: (res[key]['Q_unknown'] / res[key]['Q_total'] if res[key]['Q_total'] else 0.0) for key in SOURCES}
        print('      ledger ratio: ' + ', '.join(f'{k}: {v:.4f}' for k, v in ratios.items())
              + ' | unknown weight fraction: ' + ', '.join(f'{k}: {v:.4f}' for k, v in unk.items())
              + (f' | mismatch (s,T,Gamma): {res["mismatch"][:4]}' if res['mismatch'] else ''))
        if plots and outdir and res['n_cells']:
            theta_map(geom, res, f'{outdir}/{name}{tag}_lam{lam}.png', title=f'{name} lambda_n={lam} b={b:.4f} {tag}')
    out = {}
    if n_lam:
        t0 = time.time()
        edges = np.linspace(geom.Bmin_g, geom.Bmax_g, n_lam + 1); bm = 0.5 * (edges[1:] + edges[:-1]); db = edges[1] - edges[0]
        acc = {key: dict(QL=0.0, QU=0.0, QT=0.0, QUnk=0.0) for key in SOURCES}
        n_mismatch = 0; profile = []
        for b in bm:
            r = slice_at_pitch(geom, b, bins_per_cell=bins_per_cell)
            if MEM_CAP_GB and _rss_gb() > MEM_CAP_GB:
                print(f'  ABORT: resident memory {_rss_gb():.2f} GB exceeds MEM_CAP_GB={MEM_CAP_GB}', flush=True); sys.exit(2)
            n_mismatch += len(r['mismatch'])
            profile.append((b, {key: (r[key]['Q_lower'] / b ** 2, r[key]['Q_upper'] / b ** 2, r[key]['Q_total'] / b ** 2) for key in SOURCES}))
            for key in SOURCES:
                acc[key]['QL'] += r[key]['Q_lower'] * db / b ** 2; acc[key]['QU'] += r[key]['Q_upper'] * db / b ** 2
                acc[key]['QT'] += r[key]['Q_total'] * db / b ** 2; acc[key]['QUnk'] += r[key]['Q_unknown'] * db / b ** 2
        dt = time.time() - t0
        for key in SOURCES:
            a = acc[key]; D = 2 * Vh[key]
            out[key] = dict(f_lower=a['QL'] / D, f_upper=a['QU'] / D, trapped=a['QT'] / D, unknown=a['QUnk'] / D)
            print(f'  [{key}] f_lower={a["QL"]/D:.4f} f_upper={a["QU"]/D:.4f} trapped fraction={a["QT"]/D:.4f} unknown weight={a["QUnk"]/D:.4f}')
        print(f'  ({n_lam} slices in {dt:.1f}s, {dt/n_lam:.2f}s/slice; geometry {geom.t_geom:.1f}s; total {dt+geom.t_geom:.1f}s; '
              f'{n_mismatch} surface mismatches over all slices)')
        out['profile'] = profile
    return geom, out


if __name__ == '__main__':
    name = sys.argv[1] if len(sys.argv) > 1 else 'DMerc'
    n_s = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    n_closed = int(sys.argv[3]) if len(sys.argv) > 3 else 128
    n_lam = int(sys.argv[4]) if len(sys.argv) > 4 else 0
    s_min = float(sys.argv[5]) if len(sys.argv) > 5 else 0.02
    extra_inner = int(sys.argv[6]) if len(sys.argv) > 6 else 0   # extra surfaces inserted between s_min and the first grid point
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out')
    os.makedirs(outdir, exist_ok=True)
    base = np.linspace(0.02, 1.0, n_s)
    if extra_inner:
        inner = np.geomspace(s_min, base[0], extra_inner + 1)[:-1]
        s_grid = np.r_[inner, base]
    else:
        s_grid = np.linspace(s_min, 1.0, n_s)
    bpc = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    tag = f'_ns{len(s_grid)}_nu{n_closed}_smin{s_min}' + (f'_bpc{bpc}' if bpc != 1.0 else '')
    run(name, s_grid, n_closed / TWO_PI, n_lam, outdir=outdir, tag=tag, bins_per_cell=bpc)
