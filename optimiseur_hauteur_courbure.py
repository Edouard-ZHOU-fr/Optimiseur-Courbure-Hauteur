#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Optimiseur d'altitude (axe Z) d'une carte Lanelet2 a partir d'une trajectoire
enregistree, avec correction de courbure de la trajectoire bornee par la
cinematique du vehicule (angle de braquage maximum des roues avant).

Remplace optimiseur_hauteur_KNN.py. Corrections apportees :
  * alignement du lissage par convolution (l'ancien code decalait le signal
    d'un echantillon), fenetre exprimee en metres et non en nombre de points ;
  * appariement par KD-tree au lieu d'une recherche exhaustive O(N*M) ;
  * lecture des tags OSM par cle et non par position ;
  * chemins et parametres dans un fichier de configuration ;
  * precision d'ecriture controlee, plus de flottants a 17 chiffres ;
  * pas d'affichage bloquant, les figures sont enregistrees en PNG.

Chaine de traitement :
  1. chargement carte Lanelet2 + trajectoire CSV
  2. nettoyage : decimation au pas minimum, detection des rebroussements
     (marche arriere / saut de localisation) et decoupage en troncons
  3. reechantillonnage a pas constant par spline cubique
  4. correction de courbure sous contrainte d'ecart lateral maximum
  5. recalcul du cap (yaw)
  6. lissage de l'altitude
  7. appariement KD-tree trajectoire <-> noeuds Lanelet2
  8. ecriture de la carte, de la trajectoire corrigee et du rapport

Tous les parametres sont lus dans config/config.yaml ; le script ne prend
aucun argument en ligne de commande.

Mode AOI (aoi: true) : aucune trajectoire n'est lue. Les lanelets entierement
couverts par une zone AOI de la carte (way subtype=AOI) sont corriges
directement, raccordes en douceur a leurs voisins hors AOI :
  A1. lecture des zones AOI, selection des lanelets, insertion eventuelle de
      noeuds dans leurs bornes (mode_aoi.noeuds_ajoutes_max)
  A2. par deplacement des noeuds des bornes, sous contrainte d'ecart
      lateral maximum : lissage lateral global des centrelignes (courbure a
      variation continue), puis correction locale des depassements de k_max
  A3. lissage de l'altitude des bornes (pas de recalage sur une trajectoire)
  A4. ecriture de la carte et du rapport
"""

import math
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

from utils.config import CONFIG_DEFAUT, RACINE, charger_config
from utils.journal import Journal
from utils.vehicule import (charger_vehicule, courbure_max_vehicule,
                            lister_vehicules)
from utils.geometrie import (cap_depuis_points, courbure_menger, longueur,
                             polygone_contient, polygones_se_touchent)
from utils.trajectoire import lisser_altitude, preparer_trajectoire
from utils.courbure import appliquer_correction_courbure
from utils.aoi import (Centrelignes, SOUS_TYPES_NON_ROULABLES,
                       corriger_courbure_aoi, densifier_bornes, jacobien_geo,
                       lire_topologie, lisser_altitude_graphe, orienter_bornes)
from utils.entrees_sorties import (apparier, charger_carte, charger_trajectoire,
                                   ecrire_carte, ecrire_trajectoire)
from utils.figures import tracer, tracer_aoi


# --------------------------------------------------------------------------- #
#  Mode AOI
# --------------------------------------------------------------------------- #

def traiter_aoi(a, jrn, arbre, noeuds, pts_carte, k_max, sortie_carte):
    """Mode AOI : pas de trajectoire. Renvoie la liste des residus."""
    jrn.titre("Etape A1 : zones AOI et selection des lanelets")
    ways, tags_ways, lanelets = lire_topologie(arbre)
    xy_de = {int(i): p for i, p in zip(pts_carte[:, 0], pts_carte[:, 1:3])}
    z_de = {int(i): z for i, z in zip(pts_carte[:, 0], pts_carte[:, 3])}

    def xy(ids):
        return np.array([xy_de[i] for i in ids])

    zones = []
    for wid, refs in ways.items():
        if tags_ways[wid].get("subtype", "").upper() != "AOI":
            continue
        refs = refs[:-1] if len(refs) > 1 and refs[0] == refs[-1] else refs
        if len(refs) >= 3 and all(i in xy_de for i in refs):
            zones.append((wid, xy(refs)))
    if not zones:
        raise SystemExit("Mode AOI : aucune zone AOI (way subtype=AOI) dans %s."
                         % a.carte)

    roulables, voies = {}, {}
    for lid, gw, dw, st in lanelets:
        g, d = np.array(ways[gw]), np.array(ways[dw])
        if (st in SOUS_TYPES_NON_ROULABLES or len(g) < 2 or len(d) < 2
                or not all(i in xy_de for i in np.concatenate([g, d]))
                or longueur(xy(g)) < 1e-3 or longueur(xy(d)) < 1e-3):
            continue
        roulables[lid], voies[lid] = (g, d), (gw, dw)

    # seuls les lanelets ENTIEREMENT couverts par une AOI sont traites ; ceux
    # qui la traversent ou y debordent sont signales et laisses tels quels
    selection, n_partiels = set(), 0
    for wid, poly in zones:
        anneaux = {lid: np.vstack([xy(g), xy(d)[::-1]])
                   for lid, (g, d) in roulables.items()}
        touches = [lid for lid, r in anneaux.items() if polygones_se_touchent(r, poly)]
        dedans = [lid for lid in touches if polygone_contient(poly, anneaux[lid])]
        selection.update(dedans)
        n_partiels += len(touches) - len(dedans)
        aire = 0.5 * abs(np.dot(poly[:, 0], np.roll(poly[:, 1], -1))
                         - np.dot(poly[:, 1], np.roll(poly[:, 0], -1)))
        jrn("  AOI way %-8d : %d sommets, %.0f m2, %d lanelet(s) entierement "
            "couvert(s), %d partiellement (ignores)"
            % (wid, len(poly), aire, len(dedans), len(touches) - len(dedans)))
    if not selection:
        tous = np.vstack([xy(np.concatenate(b)) for b in roulables.values()])
        dmin = min(cKDTree(tous).query(poly)[0].min() for _, poly in zones)
        raise SystemExit(
            "Mode AOI : aucun lanelet roulable n'est entierement couvert par une "
            "zone AOI (%d partiellement couvert(s) ; noeud de lanelet le plus "
            "proche a %.1f m d'un sommet AOI). Agrandir l'AOI pour englober "
            "les lanelets a traiter." % (n_partiels, dmin))

    bornes = {lid: orienter_bornes(*roulables[lid], xy) for lid in roulables}
    debuts, fins = {}, {}
    for lid, (g, d) in bornes.items():
        debuts.setdefault((g[0], d[0]), []).append(lid)
        fins.setdefault((g[-1], d[-1]), []).append(lid)
    contexte = set()
    for lid in selection:
        g, d = bornes[lid]
        contexte.update(fins.get((g[0], d[0]), []))
        contexte.update(debuts.get((g[-1], d[-1]), []))
    contexte -= selection

    if a.noeuds_ajoutes_max > 0:
        n_ins, n_b, n_part = densifier_bornes(
            arbre, noeuds, ways, voies, xy_de, z_de, bornes, selection,
            lanelets, debuts, fins, a.noeuds_ajoutes_max, a.pas_reech,
            a.precision)
        jrn("  noeuds inseres          : %d sur %d borne(s) (max %d par borne)%s"
            % (n_ins, n_b, a.noeuds_ajoutes_max,
               ", %d borne(s) partagee(s) hors AOI non densifiee(s)" % n_part
               if n_part else ""))

    # un noeud partage avec un lanelet roulable hors selection reste fixe
    fixes = {int(i) for lid, b in roulables.items() if lid not in selection
             for i in np.concatenate(b)}
    cl = Centrelignes({lid: bornes[lid] for lid in sorted(selection | contexte)},
                      xy, a.pas_reech)
    a_bouger = {int(i) for lid in selection for i in np.concatenate(bornes[lid])} - fixes
    mobile = np.array([i in a_bouger for i in cl.ids])
    jrn("  lanelets selectionnes   : %d  (%s)"
        % (len(selection), " ".join(str(i) for i in sorted(selection))))
    jrn("  lanelets voisins figes  : %d (predecesseurs/successeurs hors AOI)"
        % len(contexte))
    jrn("  noeuds de bornes        : %d mobiles, %d figes (partages hors selection)"
        % (mobile.sum(), (~mobile).sum()))

    P0 = xy(cl.ids)
    z0 = np.array([z_de[i] for i in cl.ids])
    fenetre_lat = 0.0 if a.sans_lissage_lateral else a.fenetre_laterale
    if a.sans_courbure and fenetre_lat <= 0:
        jrn.titre("Etape A2 : lissage lateral et correction de courbure DESACTIVES")
        P, residus, k_av, k_ap = P0, [], None, None
    else:
        P, residus, k_av, k_ap = corriger_courbure_aoi(
            cl, P0, mobile, selection, contexte, k_max, a.ecart_max, a.pas_reech,
            a.base_courbure, a.empattement, fenetre_lat, not a.sans_courbure,
            jrn)

    jrn.titre("Etape A3 : lissage de l'altitude des bornes")
    if a.sans_lissage_z:
        jrn("  lissage de l'altitude desactive")
        z = z0
    else:
        aretes = {tuple(sorted((cl.col[int(u)], cl.col[int(v)])))
                  for lid in selection | contexte for b in bornes[lid]
                  for u, v in zip(b[:-1], b[1:]) if u != v}
        z = lisser_altitude_graphe(np.array(sorted(aretes)), P, z0, mobile,
                                   a.fenetre_z)
        dz = np.abs(z - z0)[mobile]
        jrn("  lissage le long des bornes : fenetre %.1f m, aucune donnee de "
            "trajectoire" % a.fenetre_z)
        jrn("  ecart d'altitude           : moyen %.3f m, max %.3f m"
            % (dz.mean() if len(dz) else 0.0, dz.max(initial=0.0)))

    seuil = 0.5 * 10.0 ** -a.precision
    dico = {int(i): float(zz) for i, zz, zz0 in zip(cl.ids, z, z0)
            if abs(zz - zz0) >= seuil}
    dxy = P - P0
    bouge = np.hypot(*dxy.T) >= seuil
    dico_xy = {}
    if bouge.any():
        J, res, n_ecart = jacobien_geo(noeuds, cl.ids)
        if J is None:
            jrn("  lat/lon non mis a jour (attributs absents)")
        else:
            jrn("  lat/lon : jacobien local ajuste sur les bornes, residu max "
                "%.3f m%s" % (res, ", %d noeud(s) incoherent(s) ecarte(s)"
                              % n_ecart if n_ecart else ""))
        for i, (x, y), d in zip(cl.ids[bouge], P[bouge], dxy[bouge]):
            dico_xy[int(i)] = (x, y, None if J is None else J @ d)

    ecrire_carte(arbre, noeuds, dico, str(sortie_carte), a.precision, jrn,
                 dico_xy, titre="Etape A4 : ecriture des resultats")
    if a.graphique:
        tracer_aoi(str(sortie_carte.with_suffix(".png")), zones, bornes,
                   selection, xy_de, cl, P, k_av, k_ap, k_max, z0[mobile],
                   z[mobile], residus, jrn)
    return residus


# --------------------------------------------------------------------------- #
#  Programme principal
# --------------------------------------------------------------------------- #

def traiter_trajectoire(a, jrn, pts_traj, arbre, noeuds, pts_carte, k_max,
                        sortie_carte, sortie_traj):
    """Mode trajectoire, etapes 2 a 8. Renvoie la liste des residus."""
    # 2-3 ------------------------------------------------------------------
    troncons = preparer_trajectoire(pts_traj, a.pas_min, a.seuil_rebroussement,
                                    a.pas_reech, jrn)

    # 4 --------------------------------------------------------------------
    if a.sans_courbure:
        jrn.titre("Etape 4 : correction de courbure DESACTIVEE")
        residus = []
    else:
        troncons, residus = appliquer_correction_courbure(
            troncons, k_max, a.ecart_max, a.pas_reech, a.empattement, jrn)

    # 5-6 ------------------------------------------------------------------
    jrn.titre("Etape 5-6 : cap et lissage de l'altitude")
    noyau = max(3, int(round(a.fenetre_z / a.pas_reech)) | 1)
    z_brut = np.concatenate([t[:, 2] for t in troncons])
    finaux, yaws = [], []
    for t in troncons:
        z = t[:, 2] if a.sans_lissage_z else lisser_altitude(t[:, 2], noyau)
        finaux.append(np.column_stack([t[:, 0], t[:, 1], z]))
        yaws.append(cap_depuis_points(t[:, :2]))
    pts_final = np.vstack(finaux)
    yaw = np.concatenate(yaws)
    if a.sans_lissage_z:
        jrn("  lissage de l'altitude desactive")
    else:
        jrn("  lissage altitude : fenetre %.1f m = %d echantillons (phase nulle)"
            % (a.fenetre_z, noyau))
        jrn("  amplitude Z : %.2f m avant, %.2f m apres"
            % (np.ptp(z_brut), np.ptp(pts_final[:, 2])))
    jrn("  cap recalcule sur %d points" % len(yaw))

    # 7 --------------------------------------------------------------------
    dico = apparier(pts_final, pts_carte, a.distance_max, a.offset_capteur,
                    a.ecart_min_z, a.ele_nul_mediane, jrn)

    # 8 --------------------------------------------------------------------
    ecrire_carte(arbre, noeuds, dico, str(sortie_carte), a.precision, jrn)
    ecrire_trajectoire(sortie_traj, pts_final, yaw)
    jrn("  trajectoire ecrite : %s (%d points)" % (sortie_traj, len(pts_final)))

    if a.graphique:
        k = np.concatenate([courbure_menger(t[:, :2]) for t in troncons])
        tracer(str(sortie_carte.with_suffix(".png")), z_brut,
               pts_final[:, 2], k, k_max, jrn)

    return residus


def main(chemin_config=CONFIG_DEFAUT):
    a = charger_config(chemin_config)
    if a.lister_vehicules:
        return lister_vehicules(a.vehicule)
    jrn = Journal(verbeux=not a.silencieux)
    debut = time.time()

    sortie_carte = a.sortie_carte or a.carte.with_name(a.carte.stem + "_OUT.osm")
    sortie_traj = a.sortie_traj or a.traj.with_name(a.traj.stem + "_OUT.csv")
    # rapport horodate, dans log/ par defaut : un lancement n'ecrase jamais
    # le rapport d'un lancement precedent
    horodatage = time.strftime("%Y%m%d_%H%M%S", time.localtime(debut))
    rapport = a.rapport or RACINE / "log" / (a.carte.stem + "_rapport.txt")
    rapport = rapport.with_name("%s_%s%s" % (rapport.stem, horodatage,
                                             rapport.suffix or ".txt"))
    rapport.parent.mkdir(parents=True, exist_ok=True)

    for f in ((a.carte,) if a.aoi else (a.traj, a.carte)):
        if not f.is_file():
            raise SystemExit("Fichier introuvable : %s" % f)

    # 1 --------------------------------------------------------------------
    jrn.titre("Etape 1 : chargement des donnees")
    jrn("  date          : %s" % time.strftime("%Y-%m-%d %H:%M:%S",
                                               time.localtime(debut)))
    jrn("  configuration : %s" % a.config)
    if a.aoi:
        jrn("  mode AOI : aucune trajectoire lue")
        pts_traj = None
    else:
        pts_traj = charger_trajectoire(a.traj, jrn)
    arbre_xml, noeuds, pts_carte = charger_carte(a.carte, jrn)

    # modele vehicule : surcharges de config.yaml > fichier vehicule > valeurs integrees
    veh, etiquette = charger_vehicule(a.vehicule, a.modele, jrn)
    surcharges = []
    for cli, cle in (("angle_braquage_max", "angle_braquage_max"),
                     ("empattement", "empattement"),
                     ("marge_courbure", "marge_courbure"),
                     ("offset_capteur", "hauteur_capteur")):
        val = getattr(a, cli)
        if val is not None:
            veh[cle] = val
            surcharges.append(cli)
    a.angle_braquage_max = veh["angle_braquage_max"]
    a.empattement = veh["empattement"]
    a.marge_courbure = veh["marge_courbure"]
    a.offset_capteur = veh["hauteur_capteur"]

    k_max = courbure_max_vehicule(a.angle_braquage_max, a.empattement,
                                  a.marge_courbure)
    jrn("")
    jrn("  modele vehicule : %s" % etiquette)
    if surcharges:
        jrn("    surcharge config.yaml : %s" % ", ".join(surcharges))
    jrn("    empattement %.3f m, braquage max %.2f deg, marge %.2f"
        % (a.empattement, a.angle_braquage_max, a.marge_courbure))
    jrn("    hauteur capteur %.3f m%s" % (a.offset_capteur,
                                    " (non utilisee en mode AOI)" if a.aoi else ""))
    r_phys = a.empattement / math.tan(math.radians(a.angle_braquage_max))
    jrn("    -> rayon de braquage physique %.2f m" % r_phys)
    jrn("    -> rayon impose %.2f m, courbure max %.4f 1/m"
        % (1.0 / k_max, k_max))

    if a.aoi:
        if a.base_courbure is None:
            a.base_courbure = a.empattement
        residus = traiter_aoi(a, jrn, arbre_xml, noeuds, pts_carte, k_max,
                              sortie_carte)
    else:
        residus = traiter_trajectoire(a, jrn, pts_traj, arbre_xml, noeuds,
                                      pts_carte, k_max, sortie_carte,
                                      sortie_traj)

    jrn("")
    jrn("  duree totale : %.2f s" % (time.time() - debut))
    jrn.ecrire(rapport)
    if not a.silencieux:
        print("\n  rapport ecrit : %s" % rapport)
    return 0 if not residus else 2


if __name__ == "__main__":
    if len(sys.argv) > 1:
        sys.exit("Ce script ne prend aucun argument : les parametres se "
                 "reglent dans %s" % CONFIG_DEFAUT)
    sys.exit(main())
