# -*- coding: utf-8 -*-
"""Correction de courbure d'une polyligne sous contrainte d'ecart lateral
(mode trajectoire), et briques communes au mode AOI."""

import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

from .geometrie import courbure_menger


LAMBDA_MIN = 1e-2


LAMBDA_MAX = 1e12


N_BISSECTIONS = 40


N_GRILLE = 60


def _lisser_bande(ref, lam):
    """Minimise  sum lam_i |p_{i-1}-2p_i+p_{i+1}|^2 + sum |p_i - ref_i|^2.
    Le terme d'ancrage sur `ref` est essentiel : un lissage laplacien pur est
    un flot de raccourcissement de courbe, qui *augmente* la courbure d'un arc
    soutenu au lieu de la reduire."""
    n = len(ref)
    D = sp.diags([1.0, -2.0, 1.0], [0, 1, 2], shape=(n - 2, n))
    A = (D.T @ sp.diags(lam[1:n - 1]) @ D + sp.eye(n, format="dia")).tocsc()
    lu = spl.splu(A)
    return np.column_stack([lu.solve(ref[:, d]) for d in range(2)])


def _zones_violation(k, k_max, marge_pts, n):
    """Regroupe les points en exces de courbure en zones dilatees et fusionnees."""
    idx = np.flatnonzero(k > k_max)
    if not len(idx):
        return []
    zones, deb, fin = [], idx[0], idx[0]
    for i in idx[1:]:
        if i - fin <= 2 * marge_pts:
            fin = i
        else:
            zones.append((deb, fin))
            deb = fin = i
    zones.append((deb, fin))
    return [(max(0, a - marge_pts), min(n, b + marge_pts + 1)) for a, b in zones]


def bissection(lisser, predicat, sens, lo=LAMBDA_MIN, hi=LAMBDA_MAX):
    """Lambda extremal de [lo, hi] verifiant `predicat`, monotone en lambda.
    `lisser(lam)` renvoie la solution lissee pour ce lambda.
    sens=+1 : predicat vrai pour les grands lambda -> renvoie le plus petit
    sens=-1 : predicat vrai pour les petits lambda -> renvoie le plus grand"""
    meilleur = None
    for _ in range(N_BISSECTIONS):
        mid = math.sqrt(lo * hi)
        cand = lisser(mid)
        if predicat(cand):
            meilleur = cand
            lo, hi = (lo, mid) if sens > 0 else (mid, hi)
        else:
            lo, hi = (mid, hi) if sens > 0 else (lo, mid)
        if hi / lo < 1.01:
            break
    return meilleur


def lissage_admissible(lisser, depart, k_zone, ecart, ecart_max):
    """Cas infaisable dans le budget d'ecart. On retient, parmi les lissages
    admissibles, celui qui minimise la courbure maximale : un lissage plus
    fort n'est pas toujours meilleur, il peut pincer l'arc et aggraver la
    courbure. Le point de depart est la zone non modifiee, donc la correction
    ne peut jamais degrader le trace."""
    meilleur, k_best = depart, k_zone(depart)
    for lam in np.logspace(math.log10(LAMBDA_MIN),
                           math.log10(LAMBDA_MAX), N_GRILLE):
        cand = lisser(lam)
        if ecart(cand) > ecart_max:
            continue
        k_cand = k_zone(cand)
        if k_cand < k_best:
            meilleur, k_best = cand, k_cand
    return meilleur


def residu(x, y, k, empattement):
    return {"x": float(x), "y": float(y), "k": float(k),
            "R": 1.0 / max(k, 1e-9),
            "delta": math.degrees(math.atan(k * empattement))}


def journal_residus(jrn, residus):
    if not residus:
        jrn("  toutes les zones ont ete ramenees sous la limite cinematique.")
        return
    jrn("")
    jrn("  /!\\ %d zone(s) NON corrigeables dans le budget d'ecart :" % len(residus))
    jrn("      %-14s %-14s %-10s %-10s %-16s %s"
        % ("x", "y", "R reel(m)", "k(1/m)", "braquage requis", "lanelets"))
    for r in sorted(residus, key=lambda z: -z["k"]):
        jrn("      %-14.2f %-14.2f %-10.2f %-10.4f %-16s %s"
            % (r["x"], r["y"], r["R"], r["k"], "%.1f deg" % r["delta"],
               " ".join(str(i) for i in r.get("lanelets", ()))))
    jrn("      -> soit augmenter courbure.ecart_max, soit revoir le trace de la carte,")
    jrn("         soit verifier les parametres vehicule.")


def corriger_courbure(troncon, k_max, ecart_max, pas, empattement, jrn):
    """Reduit la courbure sous k_max, sans deplacer un point de plus de
    `ecart_max`. Le lambda est cherche par bissection donc la contrainte
    d'ecart est respectee par construction (aucun ecretage a posteriori, qui
    recreerait des angles vifs)."""
    ref = troncon[:, :2].copy()
    cur = ref.copy()
    k = courbure_menger(cur)
    marge_pts = max(5, int(round(3.0 / (k_max * pas))))
    zones = _zones_violation(k, k_max, marge_pts, len(cur))
    if not zones:
        return cur, []

    residus = []
    for a, b in zones:
        sub = ref[a:b]
        if len(sub) < 5:
            continue
        # champ de poids en cosinus sureleve : bords adoucis, coeur a 1
        m = b - a
        rampe = min(marge_pts, m // 3)
        poids = np.ones(m)
        if rampe > 0:
            t = np.linspace(0.0, math.pi, 2 * rampe)
            taper = 0.5 * (1.0 - np.cos(t))
            poids[:rampe] = taper[:rampe]
            poids[-rampe:] = taper[rampe:]

        def k_zone(p):
            return courbure_menger(p)[1:-1].max(initial=0.0)

        def ecart(p):
            return np.hypot(*(p - sub).T).max(initial=0.0)

        def lisser(lam):
            return _lisser_bande(sub, poids * lam)

        # plus petit lissage qui ramene la courbure sous la limite
        sol = bissection(lisser, lambda p: k_zone(p) <= k_max, +1)
        if sol is not None and ecart(sol) <= ecart_max:
            cur[a:b] = sol
            continue
        cur[a:b] = lissage_admissible(lisser, sub, k_zone, ecart, ecart_max)
        k_reste = k_zone(cur[a:b])
        if k_reste > k_max:
            c = (a + b) // 2
            residus.append(residu(troncon[c, 0], troncon[c, 1], k_reste,
                                   empattement))
    return cur, residus


def appliquer_correction_courbure(troncons, k_max, ecart_max, pas, empattement, jrn):
    jrn.titre("Etape 4 : correction de courbure (limite cinematique)")
    jrn("  courbure max admissible : %.4f 1/m  (rayon min %.2f m)"
        % (k_max, 1.0 / k_max))
    jrn("  ecart lateral autorise  : %.2f m" % ecart_max)

    sortie, residus, n_av, n_ap, ecarts = [], [], 0, 0, []
    for tr in troncons:
        n_av += int((courbure_menger(tr[:, :2]) > k_max).sum())
        xy, res = corriger_courbure(tr, k_max, ecart_max, pas, empattement, jrn)
        n_ap += int((courbure_menger(xy) > k_max).sum())
        ecarts.append(np.hypot(*(xy - tr[:, :2]).T))
        residus += res
        sortie.append(np.column_stack([xy, tr[:, 2]]))

    ecarts = np.concatenate(ecarts) if ecarts else np.zeros(1)
    jrn("  points en exces de courbure : %d -> %d" % (n_av, n_ap))
    jrn("  deplacement lateral         : moyen %.3f m, max %.3f m"
        % (ecarts.mean(), ecarts.max()))
    journal_residus(jrn, residus)
    return sortie, residus
