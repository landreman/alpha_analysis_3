"""Per-surface Fourier tables and periodic cubic interpolation (DESIGN §4.1).

B has the input file's field units; theta and phi=N_fp*zeta are radians.
DB and D2B are derivatives with respect to zeta along a field line.
The FFT construction and interpolation are adapted from prototype/proto_f.py;
sine modes, resolution validation and broadcast semantics are added here.
"""

import numpy as np
from scipy.ndimage import map_coordinates, spline_filter
from scipy.optimize import minimize

from .conventions import field_line_direction

TWO_PI = 2 * np.pi


def global_extrema(field, radial_samples=65, angular_samples=64):
    """Numerical global extrema on s in [0,1], including the excluded axis.

    A coarse Fourier grid seeds exact-series optimization in all three
    coordinates. These are numerical estimates, not rigorous enclosures.
    """
    m, n = np.asarray(field.xm, int), np.asarray(field.xn, int) // field.nfp
    if angular_samples <= 2 * max(np.max(np.abs(m)), np.max(np.abs(n))):
        raise ValueError("extrema angular resolution is below the Fourier bandwidth")
    candidates = [[], []]
    for s in np.linspace(0, 1, radial_samples):
        spectrum = np.zeros((angular_samples, angular_samples), complex)
        np.add.at(
            spectrum,
            (m % angular_samples, (-n) % angular_samples),
            field.bmnc(s) - 1j * field.bmns(s),
        )
        B = np.fft.ifft2(spectrum).real * angular_samples**2
        for j, idx in enumerate((np.argmin(B), np.argmax(B))):
            a, c = np.unravel_index(idx, B.shape)
            candidates[j].append(
                (
                    B[a, c],
                    np.array(
                        [s, TWO_PI * a / angular_samples, TWO_PI * c / angular_samples]
                    ),
                )
            )

    def exact(x):
        s, theta, phi = x
        phase = m * theta - n * phi
        return float(np.cos(phase) @ field.bmnc(s) + np.sin(phase) @ field.bmns(s))

    result = {}
    for j, sign in enumerate((1, -1)):
        seeds = sorted(candidates[j], key=lambda r: sign * r[0])[:8]
        best = None
        for value, point in seeds:
            opt = minimize(
                lambda x: sign * exact(x),
                point,
                method="L-BFGS-B",
                bounds=[(0, 1), (None, None), (None, None)],
                options={"ftol": 1e-14, "gtol": 1e-9},
            )
            candidate = (sign * float(opt.fun), opt.x)
            if best is None or sign * candidate[0] < sign * best[0]:
                best = candidate
        result["min" if j == 0 else "max"] = {
            "B": best[0],
            "s_theta_phi": best[1].tolist(),
        }
    return result


class Surface:
    """A flux surface at s=rho²; no pitch-dependent state is stored."""

    def __init__(self, field, s, nt=512, nphi=512):
        self.s, self.nfp = float(s), int(field.nfp)
        self.nt, self.nphi = int(nt), int(nphi)
        self.L = TWO_PI / self.nfp
        self.iota = float(field.iota(s))
        self.G, self.I = float(field.G(s)), float(field.I(s))
        self.C = self.G + self.iota * self.I
        self.sign = float(field_line_direction(self.G, self.I, self.iota))
        self.m, self.n = np.asarray(field.xm, int), np.asarray(field.xn, int)
        if np.any(self.n % self.nfp):
            raise ValueError("Fourier n must be a multiple of N_fp")
        if nt <= 2 * np.max(np.abs(self.m)) or nphi <= 2 * np.max(
            np.abs(self.n // self.nfp)
        ):
            raise ValueError("table resolution must exceed twice every Fourier mode")
        self.c = np.asarray(field.bmnc(s), float).ravel()
        self.d = np.asarray(field.bmns(s), float).ravel()
        k = self.iota * self.m - self.n
        amplitude = self.c - 1j * self.d
        arrays = []
        corrections = []
        for factor in (np.ones_like(k), 1j * k, -k * k):
            spectrum = np.zeros((nt, nphi), complex)
            np.add.at(
                spectrum,
                (self.m % nt, (-self.n // self.nfp) % nphi),
                amplitude * factor,
            )
            arrays.append(np.fft.ifft2(spectrum).real * nt * nphi)
            # Remove the leading cubic-spline defect using exact fourth
            # derivatives. The correction vanishes at table grid points.
            fm = np.fft.fftfreq(nt) * nt
            fn = np.fft.fftfreq(nphi) * nphi
            corrections.append(
                [
                    spline_filter(
                        np.fft.ifft2(spectrum * axis**4).real * nt * nphi,
                        order=3,
                        mode="grid-wrap",
                    )
                    for axis in (fm[:, None], fn[None, :])
                ]
            )
        self.B, self.DB, self.D2B = arrays
        if np.min(self.B) <= 0:
            raise ValueError("B must be strictly positive")
        self.Bmin, self.Bmax = float(self.B.min()), float(self.B.max())
        self._coeff = [spline_filter(a, order=3, mode="grid-wrap") for a in arrays]
        self._corrections = corrections

    def _interp(self, which, theta, phi):
        theta, phi = np.broadcast_arrays(theta, phi)
        coords = [theta.ravel() * self.nt / TWO_PI, phi.ravel() * self.nphi / TWO_PI]
        value = map_coordinates(
            self._coeff[which], coords, order=3, mode="grid-wrap", prefilter=False
        )
        for axis, size in enumerate((self.nt, self.nphi)):
            fraction = coords[axis] % 1
            defect = (TWO_PI / size) ** 4 / 24 * fraction**2 * (1 - fraction) ** 2
            value += defect * map_coordinates(
                self._corrections[which][axis],
                coords,
                order=3,
                mode="grid-wrap",
                prefilter=False,
            )
        return value.reshape(theta.shape)

    def b(self, theta, phi):
        return self._interp(0, theta, phi)

    def db(self, theta, phi):
        return self._interp(1, theta, phi)

    def d2b(self, theta, phi):
        return self._interp(2, theta, phi)

    def exact(self, theta, phi, derivative=0):
        """Exact series, bounded batches; derivative is the along-zeta order."""
        theta, phi = np.broadcast_arrays(theta, phi)
        result = np.empty(theta.size)
        k = self.iota * self.m - self.n
        amp = (self.c - 1j * self.d) * (1j * k) ** derivative
        batch = max(1, 100_000 // len(self.m))
        for start in range(0, theta.size, batch):
            sl = slice(start, start + batch)
            phase = (
                theta.ravel()[sl, None] * self.m
                - phi.ravel()[sl, None] * self.n / self.nfp
            )
            result[sl] = np.real(np.exp(1j * phase) @ amp)
        return result.reshape(theta.shape)

    def release_derivatives(self):
        """Release extraction-only tables after nodes and chains are built."""
        self.DB = self.D2B = None
        self._coeff[1:] = [None, None]
        self._corrections[1:] = [None, None]
