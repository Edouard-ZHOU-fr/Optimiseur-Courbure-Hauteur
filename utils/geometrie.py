# -*- coding: utf-8 -*-
"""Geometrie plane : courbure discrete, cap, polygones, polylignes."""

import numpy as np


def menger(p1, p2, p3):
    """Courbure (1/m) du cercle passant par les triplets (p1, p2, p3)."""
    v1 = p2 - p1
    v2 = p3 - p2
    v3 = p3 - p1
    den = (np.hypot(v1[:, 0], v1[:, 1]) * np.hypot(v2[:, 0], v2[:, 1])
           * np.hypot(v3[:, 0], v3[:, 1]))
    aire2 = np.abs(v1[:, 0] * v3[:, 1] - v1[:, 1] * v3[:, 0])
    k = np.zeros(len(p1))
    ok = den > 1e-9
    k[ok] = 2.0 * aire2[ok] / den[ok]
    return k


def courbure_menger(pts):
    """Courbure discrete (1/m) aux points interieurs. Renvoie un tableau de
    taille len(pts), nul aux deux extremites."""
    k = np.zeros(len(pts))
    if len(pts) >= 3:
        k[1:-1] = menger(pts[:-2], pts[1:-1], pts[2:])
    return k


def cap_depuis_points(pts):
    """Cap (rad) par differences centrees sur les points interieurs."""
    n = len(pts)
    yaw = np.zeros(n)
    if n < 2:
        return yaw
    d = np.empty_like(pts)
    d[1:-1] = pts[2:] - pts[:-2]
    d[0] = pts[1] - pts[0]
    d[-1] = pts[-1] - pts[-2]
    return np.arctan2(d[:, 1], d[:, 0])


def dans_polygone(pts, poly):
    """Test pair-impair, vectorise : pts (n, 2), poly (m, 2) -> (n,) bool."""
    x, y = pts[:, :1], pts[:, 1:2]
    x1, y1 = poly[:, 0], poly[:, 1]
    x2, y2 = np.roll(x1, -1), np.roll(y1, -1)
    dy = np.where(y2 == y1, 1.0, y2 - y1)
    croise = ((y1 > y) != (y2 > y)) & (x < (x2 - x1) * (y - y1) / dy + x1)
    return croise.sum(axis=1) % 2 == 1


def polygones_se_touchent(p1, p2):
    """Vrai si les polygones se recouvrent ou si leurs bords se coupent."""
    if dans_polygone(p1, p2).any() or dans_polygone(p2, p1).any():
        return True
    return _bords_se_coupent(p1, p2)


def polygone_contient(poly, pts):
    """Vrai si l'anneau `pts` est entierement dans `poly` : tous ses sommets
    dedans et aucun de ses cotes ne traverse le bord (poly non convexe)."""
    return bool(dans_polygone(pts, poly).all()) and not _bords_se_coupent(pts, poly)


def _bords_se_coupent(p1, p2):
    a, b = p1[:, None], p2[None, :]
    r, s = np.roll(p1, -1, axis=0)[:, None] - a, np.roll(p2, -1, axis=0)[None, :] - b
    qp = b - a
    den = r[..., 0] * s[..., 1] - r[..., 1] * s[..., 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (qp[..., 0] * s[..., 1] - qp[..., 1] * s[..., 0]) / den
        u = (qp[..., 0] * r[..., 1] - qp[..., 1] * r[..., 0]) / den
    return bool(np.any((den != 0) & (t >= 0) & (t <= 1) & (u >= 0) & (u <= 1)))


def longueur(pts):
    return float(np.hypot(*np.diff(pts, axis=0).T).sum())


def abscisses(pts, t):
    """Pour chaque fraction t de la longueur d'une polyligne : indice du
    segment porteur et poids de son second sommet."""
    s = np.concatenate([[0.0], np.cumsum(np.hypot(*np.diff(pts, axis=0).T))])
    s /= s[-1]
    j = np.clip(np.searchsorted(s, t, side="right") - 1, 0, len(pts) - 2)
    ds = s[j + 1] - s[j]
    alpha = np.where(ds > 0, (t - s[j]) / np.where(ds > 0, ds, 1.0), 0.0)
    return j, np.clip(alpha, 0.0, 1.0)


def echantillonner(pts, t):
    j, alpha = abscisses(pts, t)
    return pts[j] * (1.0 - alpha)[:, None] + pts[j + 1] * alpha[:, None]
