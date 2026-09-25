# -*- coding: utf-8 -*-
"""Entrees / sorties : trajectoire CSV, carte Lanelet2 OSM, appariement
trajectoire <-> noeuds."""

import csv
import xml.etree.ElementTree as ET

import numpy as np
from scipy.spatial import cKDTree


def charger_trajectoire(chemin, jrn):
    pts = []
    with open(chemin, encoding="utf-8-sig") as f:
        lecteur = csv.reader(f, skipinitialspace=True)
        entete = next(lecteur, None)
        for n_ligne, ligne in enumerate(lecteur, start=2):
            if len(ligne) < 3 or not ligne[0].strip():
                continue
            try:
                pts.append([float(ligne[0]), float(ligne[1]), float(ligne[2])])
            except ValueError:
                jrn("  ligne %d ignoree (valeur non numerique)" % n_ligne)
    if len(pts) < 4:
        raise SystemExit("Trajectoire trop courte : %d points." % len(pts))
    jrn("  trajectoire : %s" % chemin)
    jrn("    entete=%s, %d points" % (entete, len(pts)))
    return np.asarray(pts, dtype=float)


def charger_carte(chemin, jrn):
    arbre = ET.parse(chemin)
    noeuds = arbre.getroot().findall("node")
    points, ignores = [], 0
    for nd in noeuds:
        # lecture par cle : l'ordre des tags varie selon la version de
        # VectorMap (mgrs_code en tete avant oct/2024, tag color optionnel)
        tags = {t.get("k"): t.get("v") for t in nd.findall("tag")}
        try:
            points.append([int(nd.get("id")), float(tags["local_x"]),
                           float(tags["local_y"]), float(tags["ele"])])
        except (KeyError, TypeError, ValueError):
            ignores += 1
    if not points:
        raise SystemExit("Aucun noeud exploitable : format Lanelet2 inattendu.")
    jrn("  carte Lanelet2 : %s" % chemin)
    jrn("    %d noeuds lus, %d ignores (tags local_x/local_y/ele absents)"
        % (len(points), ignores))
    return arbre, noeuds, np.asarray(points, dtype=float)


def ecrire_trajectoire(chemin, pts, yaw):
    with open(chemin, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["x", "y", "z", "yaw"])
        for (x, y, z), a in zip(pts, yaw):
            w.writerow(["%.4f" % x, "%.4f" % y, "%.4f" % z, "%.6f" % a])


def apparier(pts_traj, pts_carte, distance_max, offset, ecart_min,
             ele_nul_mediane, jrn):
    jrn.titre("Etape 7 : appariement trajectoire <-> noeuds Lanelet2")
    arbre = cKDTree(pts_traj[:, :2])
    dist, idx = arbre.query(pts_carte[:, 1:3], workers=-1)
    z_opt = pts_traj[idx, 2] - offset
    z_med = float(np.median(pts_traj[:, 2])) - offset

    voisin = dist < distance_max
    ecart = np.abs(z_opt - pts_carte[:, 3])
    a_modifier = voisin & (ecart >= ecart_min)

    dico = {int(i): float(z) for i, z in
            zip(pts_carte[a_modifier, 0], z_opt[a_modifier])}
    n_median = 0
    if ele_nul_mediane:
        # noeuds sans voisin de trajectoire et dont l'altitude est restee a 0
        secours = (~voisin) & (pts_carte[:, 3] == 0.0)
        for i in pts_carte[secours, 0]:
            dico[int(i)] = z_med
        n_median = int(secours.sum())

    jrn("  rayon de recherche       : %.1f m" % distance_max)
    jrn("  offset capteur soustrait : %.3f m" % offset)
    jrn("  noeuds avec voisin       : %d / %d (%.1f %%)"
        % (voisin.sum(), len(pts_carte), 100.0 * voisin.mean()))
    jrn("  distance d'appariement   : mediane %.2f m, p95 %.2f m"
        % (np.median(dist[voisin]) if voisin.any() else 0.0,
           np.percentile(dist[voisin], 95) if voisin.any() else 0.0))
    jrn("  ecart d'altitude corrige : mediane %.3f m, max %.3f m"
        % (np.median(ecart[a_modifier]) if a_modifier.any() else 0.0,
           ecart[a_modifier].max(initial=0.0)))
    jrn("  noeuds modifies          : %d (dont %d par mediane de secours)"
        % (len(dico), n_median))
    return dico


def ecrire_carte(arbre, noeuds, dico, chemin, precision, jrn, dico_xy=None,
                 titre="Etape 8 : ecriture des resultats"):
    """dico : {id: ele}. dico_xy : {id: (local_x, local_y, d(lat, lon))},
    d(lat, lon) a None si lat/lon ne peuvent etre mis a jour."""
    jrn.titre(titre)
    dico_xy = dico_xy or {}
    n_mod, n_xy = 0, 0
    for nd in noeuds:
        ident = int(nd.get("id"))
        if ident not in dico and ident not in dico_xy:
            continue
        valeurs = {}
        if ident in dico:
            valeurs["ele"] = dico[ident]
        if ident in dico_xy:
            x, y, dgeo = dico_xy[ident]
            valeurs.update(local_x=x, local_y=y)
            if dgeo is not None:
                for att, dv in zip(("lat", "lon"), dgeo):
                    nd.set(att, "%.11f" % (float(nd.get(att)) + dv))
            n_xy += 1
        for t in nd.findall("tag"):
            if t.get("k") in valeurs:
                t.set("v", "%.*f" % (precision, valeurs[t.get("k")]))
        n_mod += 1
    arbre.write(chemin, encoding="utf-8", xml_declaration=True)
    jrn("  carte ecrite : %s (%d noeuds sur %d modifies%s)"
        % (chemin, n_mod, len(noeuds),
           ", dont %d deplaces en x/y" % n_xy if dico_xy else ""))
    return n_mod
