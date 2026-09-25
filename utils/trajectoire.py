# -*- coding: utf-8 -*-
"""Mode trajectoire : nettoyage, decoupage en troncons, reechantillonnage
et lissage de l'altitude."""

import math

import numpy as np
from scipy.interpolate import CubicSpline


def decimer(points, pas_min):
    """Supprime les points trop proches du precedent retenu (vehicule a l'arret).
    Ces points produisent des courbures aberrantes (division par ~0)."""
    garde = [0]
    ref = points[0, :2]
    for i in range(1, len(points)):
        if math.dist(points[i, :2], ref) >= pas_min:
            garde.append(i)
            ref = points[i, :2]
    return points[garde]


def detecter_rebroussements(points, seuil_deg):
    """Indices ou le cap s'inverse brutalement : marche arriere, manoeuvre en
    plusieurs points, ou saut de localisation. La courbure n'a aucun sens sur
    un tel point de rebroussement, il faut couper la trajectoire."""
    if len(points) < 3:
        return np.empty(0, dtype=int)
    v = np.diff(points[:, :2], axis=0)
    th = np.arctan2(v[:, 1], v[:, 0])
    dth = np.abs((np.diff(th) + np.pi) % (2 * np.pi) - np.pi)
    return np.flatnonzero(dth > math.radians(seuil_deg)) + 1


def reechantillonner(troncon, pas):
    """Spline cubique parametree par l'abscisse curviligne, echantillonnee a
    pas constant. Indispensable : sans cela la courbure est incomparable d'un
    point a l'autre (le pas brut varie ici de 0.4 mm a 23 m)."""
    seg = np.hypot(*np.diff(troncon[:, :2], axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    garde = np.concatenate([[True], seg > 1e-9])
    troncon, s = troncon[garde], s[garde]
    if len(troncon) < 4 or s[-1] < 2 * pas:
        return None
    su = np.arange(0.0, s[-1], pas)
    return np.column_stack([CubicSpline(s, troncon[:, c])(su) for c in range(3)])


def preparer_trajectoire(points, pas_min, seuil_rebrous, pas_reech, jrn):
    jrn.titre("Etape 2-3 : nettoyage et reechantillonnage de la trajectoire")
    n_brut = len(points)
    points = decimer(points, pas_min)
    jrn("  decimation (pas min %.2f m)      : %d -> %d points (%d supprimes)"
        % (pas_min, n_brut, len(points), n_brut - len(points)))

    coupures = detecter_rebroussements(points, seuil_rebrous)
    jrn("  rebroussements detectes (>%d deg): %d" % (seuil_rebrous, len(coupures)))
    if len(coupures):
        for i in coupures:
            jrn("      x=%.2f y=%.2f  (marche arriere ou saut de localisation)"
                % (points[i, 0], points[i, 1]))

    bornes = np.concatenate([[0], coupures, [len(points)]])
    troncons, rejetes = [], 0
    for a, b in zip(bornes[:-1], bornes[1:]):
        r = reechantillonner(points[a:b], pas_reech)
        if r is None:
            rejetes += 1
        else:
            troncons.append(r)
    total = sum(len(t) for t in troncons)
    jrn("  troncons retenus                 : %d (%d trop courts, ignores)"
        % (len(troncons), rejetes))
    jrn("  reechantillonnage (pas %.2f m)   : %d points" % (pas_reech, total))
    if not troncons:
        raise SystemExit("Aucun troncon exploitable dans la trajectoire.")
    return troncons


def lisser_altitude(z, noyau):
    """Moyenne glissante centree, a phase nulle.

    L'ancienne implementation utilisait np.convolve(..., 'full') puis un
    decoupage errone qui decalait la sortie d'un echantillon. Ici le
    remplissage par les valeurs de bord suivi d'un mode 'valid' donne
    exactement len(z) points, correctement centres."""
    noyau = max(3, int(noyau) | 1)          # impair, >= 3
    if len(z) <= noyau:
        return np.asarray(z, dtype=float)
    demi = noyau // 2
    rempli = np.concatenate([np.full(demi, z[0]), z, np.full(demi, z[-1])])
    return np.convolve(rempli, np.ones(noyau) / noyau, mode="valid")
