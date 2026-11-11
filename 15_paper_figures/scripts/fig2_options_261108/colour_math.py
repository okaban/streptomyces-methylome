"""Colour measurements: sRGB -> linear -> XYZ(D65) -> CIELAB; CVD simulation
(Machado, Oliveira & Fernandes 2009, severity 1.0, applied in linear RGB);
colour difference CIEDE2000 (Sharma, Wu & Dalal 2005 formulation)."""
import numpy as np

def hex2rgb(h):
    h = h.lstrip('#')
    return np.array([int(h[i:i+2], 16) / 255 for i in (0, 2, 4)])

def rgb2hex(c):
    c = np.clip(c, 0, 1)
    return '#%02X%02X%02X' % tuple(int(round(v * 255)) for v in c)

def srgb_to_lin(c):
    c = np.asarray(c, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)

def lin_to_srgb(c):
    c = np.clip(np.asarray(c, float), 0, 1)
    return np.where(c <= 0.0031308, c * 12.92, 1.055 * c ** (1 / 2.4) - 0.055)

M_XYZ = np.array([[0.4124564, 0.3575761, 0.1804375],
                  [0.2126729, 0.7151522, 0.0721750],
                  [0.0193339, 0.1191920, 0.9503041]])
WP = np.array([0.95047, 1.00000, 1.08883])

def lab(rgb):
    xyz = M_XYZ @ srgb_to_lin(rgb) / WP
    f = np.where(xyz > (6 / 29) ** 3, np.cbrt(xyz), xyz / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])

def luminance(rgb):          # Rec.709 relative luminance Y, 0-1
    return float(M_XYZ[1] @ srgb_to_lin(rgb))

CVD = {  # Machado et al. 2009, severity 1.0
 'deutan': np.array([[0.367322, 0.860646, -0.227968],
                     [0.280085, 0.672501,  0.047413],
                     [-0.011820, 0.042940, 0.968881]]),
 'protan': np.array([[0.152286, 1.052583, -0.204868],
                     [0.114503, 0.786281,  0.099216],
                     [-0.003882, -0.048116, 1.051998]]),
 'tritan': np.array([[1.255528, -0.076749, -0.178779],
                     [-0.078411, 0.930809,  0.147602],
                     [0.004733,  0.691367,  0.303900]])}

def simulate(rgb, kind):
    return lin_to_srgb(CVD[kind] @ srgb_to_lin(rgb))

def greyscale(rgb):          # print-style conversion: keep relative luminance
    y = luminance(rgb)
    return lin_to_srgb(np.array([y, y, y]))

def de2000(lab1, lab2):
    L1, a1, b1 = lab1; L2, a2, b2 = lab2
    C1, C2 = np.hypot(a1, b1), np.hypot(a2, b2)
    Cb = (C1 + C2) / 2
    G = 0.5 * (1 - np.sqrt(Cb**7 / (Cb**7 + 25.0**7)))
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = np.hypot(a1p, b1), np.hypot(a2p, b2)
    h1p = np.degrees(np.arctan2(b1, a1p)) % 360 if (abs(a1p) + abs(b1)) else 0.0
    h2p = np.degrees(np.arctan2(b2, a2p)) % 360 if (abs(a2p) + abs(b2)) else 0.0
    dLp = L2 - L1
    dCp = C2p - C1p
    if C1p * C2p == 0:
        dhp = 0.0
    elif abs(h2p - h1p) <= 180:
        dhp = h2p - h1p
    elif h2p - h1p > 180:
        dhp = h2p - h1p - 360
    else:
        dhp = h2p - h1p + 360
    dHp = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dhp) / 2)
    Lbp = (L1 + L2) / 2
    Cbp = (C1p + C2p) / 2
    if C1p * C2p == 0:
        hbp = h1p + h2p
    elif abs(h1p - h2p) <= 180:
        hbp = (h1p + h2p) / 2
    elif h1p + h2p < 360:
        hbp = (h1p + h2p + 360) / 2
    else:
        hbp = (h1p + h2p - 360) / 2
    T = (1 - 0.17 * np.cos(np.radians(hbp - 30)) + 0.24 * np.cos(np.radians(2 * hbp))
         + 0.32 * np.cos(np.radians(3 * hbp + 6)) - 0.20 * np.cos(np.radians(4 * hbp - 63)))
    dtheta = 30 * np.exp(-(((hbp - 275) / 25) ** 2))
    Rc = 2 * np.sqrt(Cbp**7 / (Cbp**7 + 25.0**7))
    Sl = 1 + 0.015 * (Lbp - 50) ** 2 / np.sqrt(20 + (Lbp - 50) ** 2)
    Sc = 1 + 0.045 * Cbp
    Sh = 1 + 0.015 * Cbp * T
    Rt = -np.sin(np.radians(2 * dtheta)) * Rc
    return float(np.sqrt((dLp / Sl) ** 2 + (dCp / Sc) ** 2 + (dHp / Sh) ** 2
                         + Rt * (dCp / Sc) * (dHp / Sh)))
