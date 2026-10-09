# -*- coding: utf-8 -*-
"""Mode AOI : topologie Lanelet2, centrelignes, lissage lateral et
correction de courbure des lanelets, lissage de l'altitude des bornes,
insertion de noeuds."""

import copy
import math
import xml.etree.ElementTree as ET

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .courbure import bissection, journal_residus, lissage_admissible, residu
from .geometrie import abscisses, echantillonner, longueur, menger


POIDS_LARGEUR = 0.3     # mode AOI : poids de la conservation de largeur de voie (correction)


POIDS_LARGEUR_LISSAGE = 1.0     # idem pour le lissage lateral global


N_ANCRAGE = 30          # mode AOI : iterations du budget d'ecart par noeud


POIDS_RACCORD = 10.0    # mode AOI : poids du lissage a travers les raccords
N_PASSES_ELAGAGE = 5    # mode AOI : passes max A2 -> A3 -> elagage


# lanelets exclus de la correction : la contrainte de braquage n'y a pas de sens
SOUS_TYPES_NON_ROULABLES = {"crosswalk", "walkway", "stairs"}


def noeuds_proteges(arbre, types):
    """Ids des noeuds des ways dont le tag `type` est dans `types`, et nombre
    de ways concernes."""
    ids, n_ways = set(), 0
    for w in arbre.getroot().findall("way"):
        tags = {t.get("k"): t.get("v") for t in w.findall("tag")}
        if tags.get("type") in types:
            ids.update(int(nd.get("ref")) for nd in w.findall("nd"))
            n_ways += 1
    return ids, n_ways


def lire_topologie(arbre):
    """Ways {id: [refs]}, leurs tags, et les lanelets (id, gauche, droite,
    subtype)."""
    racine = arbre.getroot()
    ways, tags_ways = {}, {}
    for w in racine.findall("way"):
        i = int(w.get("id"))
        ways[i] = [int(n.get("ref")) for n in w.findall("nd")]
        tags_ways[i] = {t.get("k"): t.get("v") for t in w.findall("tag")}
    lanelets = []
    for r in racine.findall("relation"):
        tags = {t.get("k"): t.get("v") for t in r.findall("tag")}
        if tags.get("type") != "lanelet":
            continue
        roles = {m.get("role"): int(m.get("ref"))
                 for m in r.findall("member") if m.get("type") == "way"}
        if roles.get("left") in ways and roles.get("right") in ways:
            lanelets.append((int(r.get("id")), roles["left"], roles["right"],
                             tags.get("subtype", "")))
    return ways, tags_ways, lanelets


def orienter_bornes(g, d, xy):
    """Ids des bornes gauche/droite, orientees dans le sens de circulation.
    OSM n'impose pas que les deux ways aillent dans le meme sens (bord
    partage entre deux voies opposees) : on retourne la droite si besoin,
    puis le sens de circulation est celui ou la gauche est bien a gauche."""
    pg, pd = xy(g), xy(d)
    if (math.dist(pg[0], pd[-1]) + math.dist(pg[-1], pd[0])
            < math.dist(pg[0], pd[0]) + math.dist(pg[-1], pd[-1])):
        d, pd = d[::-1], pd[::-1]
    t = np.linspace(0.0, 1.0, 11)
    cg, cd = echantillonner(pg, t), echantillonner(pd, t)
    tang = np.diff(cg + cd, axis=0)
    lat = (cg - cd)[:-1]
    if np.sum(tang[:, 0] * lat[:, 1] - tang[:, 1] * lat[:, 0]) < 0:
        g, d = g[::-1], d[::-1]
    return g, d


class Centrelignes:
    """Centrelignes des lanelets, combinaisons lineaires figees des noeuds de
    leurs bornes : C = M @ P. Chaque centreligne est echantillonnee par
    fractions egales de la longueur des deux bornes, comme le fait
    lanelet2_extension::generateFineCenterline dans Autoware. Les poids sont
    figes sur la geometrie d'origine, le probleme reste donc lineaire en P.

    Les echantillons forment un graphe oriente : un lanelet est relie a ses
    successeurs par le point milieu des noeuds d'extremite qu'ils partagent.
    L = M_l @ P donne en chaque echantillon le vecteur gauche - droite, qui
    sert a conserver la largeur de voie."""

    def __init__(self, bornes, xy, pas):
        self.col = {}                   # id noeud -> colonne de M
        self.suivants, self.proprio = [], []
        lignes, cols, vals, signes = [], [], [], []
        extremites = {}
        for lid, (g, d) in bornes.items():
            pg, pd = xy(g), xy(d)
            n = max(3, int(math.ceil(0.5 * (longueur(pg) + longueur(pd)) / pas)) + 1)
            t = np.linspace(0.0, 1.0, n)
            jg, ag = abscisses(pg, t)
            jd, ad = abscisses(pd, t)
            chaine = []
            for k in range(n):
                if k in (0, n - 1):
                    cle = (g[0], d[0]) if k == 0 else (g[-1], d[-1])
                    v = extremites.get(cle)
                    if v is None:
                        v = extremites[cle] = self._sommet()
                        termes = [(cle[0], 1.0, 1.0), (cle[1], 1.0, -1.0)]
                    else:
                        termes = []
                else:
                    v = self._sommet()
                    termes = [(g[jg[k]], 1 - ag[k], 1.0), (g[jg[k] + 1], ag[k], 1.0),
                              (d[jd[k]], 1 - ad[k], -1.0), (d[jd[k] + 1], ad[k], -1.0)]
                for nid, w, sg in termes:
                    lignes.append(v)
                    cols.append(self.col.setdefault(int(nid), len(self.col)))
                    vals.append(w)
                    signes.append(sg)
                self.proprio[v].add(lid)
                chaine.append(v)
            for u, w in zip(chaine[:-1], chaine[1:]):
                self.suivants[u].append(w)
        forme = (len(self.suivants), len(self.col))
        vals, signes = np.asarray(vals), np.asarray(signes)
        self.M = sp.csr_matrix((0.5 * vals, (lignes, cols)), shape=forme)
        self.M_l = sp.csr_matrix((signes * vals, (lignes, cols)), shape=forme)
        self.ids = np.array(sorted(self.col, key=self.col.get))

    def _sommet(self):
        self.suivants.append([])
        self.proprio.append(set())
        return len(self.suivants) - 1

    def triplets(self, pas):
        """Triplets (a, b, c) : a est `pas` echantillons avant b, c `pas`
        apres, en suivant toutes les branches du graphe."""
        precedents = [[] for _ in self.suivants]
        for u, ws in enumerate(self.suivants):
            for w in ws:
                precedents[w].append(u)

        def atteints(v, adj):
            front = {v}
            for _ in range(pas):
                front = {w for u in front for w in adj[u]}
            return front

        out = [(a, v, c) for v in range(len(self.suivants))
               for a in atteints(v, precedents) for c in atteints(v, self.suivants)]
        return np.array(out, dtype=int).reshape(-1, 3)

    def chemins(self, n):
        """Tous les chemins orientes de n echantillons consecutifs."""
        out = [[v] for v in range(len(self.suivants))]
        for _ in range(n - 1):
            out = [c + [w] for c in out for w in self.suivants[c[-1]]]
        return np.array(out, dtype=int).reshape(-1, n)


def _grappes(pts, seuil):
    """Regroupement par lien simple : liste des indices de chaque grappe."""
    paires = cKDTree(pts).query_pairs(seuil, output_type="ndarray")
    g = sp.coo_matrix((np.ones(len(paires)), (paires[:, 0], paires[:, 1])),
                      shape=(len(pts), len(pts)))
    n, lab = connected_components(g, directed=False)
    return [np.flatnonzero(lab == i) for i in range(n)]


def _operateur(chemins, coefs, n):
    """Matrice creuse des differences finies `coefs` le long des chemins."""
    m = len(chemins)
    return sp.csr_matrix((np.tile(coefs, m),
                          (np.repeat(np.arange(m), len(coefs)), chemins.ravel())),
                         shape=(m, n))


def _operateur_normal(C, chemins, ordre, pas, n):
    """Derivee d'ordre `ordre` de la centreligne le long des chemins, projetee
    sur la normale, en differences divisees sur les abscisses curvilignes de
    C (geometrie d'origine). Mise a l'echelle par pas**ordre : sur un pas
    uniforme egal a `pas`, on retrouve les coefficients usuels (-1, 3, -3, 1
    a l'ordre 3).

    Une difference a coefficients fixes suppose un pas uniforme ; or les
    echantillons ne le sont pas (pas propre a chaque lanelet, bornes de
    formes differentes) : son residu est non nul meme sur une droite, et
    dirige le long de la route. De plus d3r/ds3 = k' N - k^2 T : la
    composante tangentielle ne mesure pas la variation de courbure. La
    penaliser fait glisser les noeuds le long de la route, ce qui consomme le
    budget d'ecart sans rien lisser. Seule la composante normale (k' a
    l'ordre 3) est donc retenue ; la normale (celle de la corde du chemin)
    est figee sur la geometrie d'origine pour que le probleme reste lineaire.
    Renvoie (Nx, Ny), creuses (chemins x n) : residu = Nx @ Cx + Ny @ Cy."""
    m, k = chemins.shape
    Cc = C[chemins]
    seg = np.linalg.norm(np.diff(Cc, axis=1), axis=2)
    s = np.concatenate([np.zeros((m, 1)), np.cumsum(seg, axis=1)], axis=1)
    ok = (seg > 1e-6).all(axis=1)          # echantillons confondus : ligne nulle
    coef = np.zeros((m, k))
    for i in range(k):
        prod = np.ones(m)
        for j in range(k):
            if j != i:
                prod = prod * (s[:, i] - s[:, j])
        coef[:, i] = np.where(ok, math.factorial(ordre) * pas ** ordre
                              / np.where(ok, prod, 1.0), 0.0)
    corde = Cc[:, -1] - Cc[:, 0]
    t = corde / np.maximum(np.linalg.norm(corde, axis=1), 1e-9)[:, None]
    lignes = np.repeat(np.arange(m), k)

    def projete(composante):
        return sp.csr_matrix(((coef * composante[:, None]).ravel(),
                              (lignes, chemins.ravel())), shape=(m, n))
    return projete(-t[:, 1]), projete(t[:, 0])


def _taux_courbure(C, ts, suivants, sommets=None):
    """|dk/ds| (1/m2) le long des aretes du graphe, k etant pris en chaque
    echantillon central des triplets controles. `sommets` : ne garder que
    les aretes qui en touchent un."""
    k = np.full(len(C), np.nan)
    np.fmax.at(k, ts[:, 1], menger(C[ts[:, 0]], C[ts[:, 1]], C[ts[:, 2]]))
    u = np.array([a for a, ws in enumerate(suivants) for _ in ws], dtype=int)
    w = np.array([b for ws in suivants for b in ws], dtype=int)
    ok = ~np.isnan(k[u]) & ~np.isnan(k[w])
    if sommets is not None:
        ok &= np.isin(u, sommets) | np.isin(w, sommets)
    ds = np.maximum(np.hypot(*(C[w] - C[u]).T), 1e-6)
    return np.abs(k[w] - k[u])[ok] / ds[ok]


def _cassures(C, t1):
    """Angle de cassure (deg) de la centreligne aux triplets t1."""
    v1, v2 = C[t1[:, 1]] - C[t1[:, 0]], C[t1[:, 2]] - C[t1[:, 1]]
    return np.degrees(np.abs(np.arctan2(v1[:, 0] * v2[:, 1] - v1[:, 1] * v2[:, 0],
                                        (v1 * v2).sum(axis=1))))


def controle_courbure(bornes, xy, pas, base, selection, contexte, k_max):
    """Courbure des centrelignes de la geometrie finale, re-echantillonnees
    sur celle-ci (comme Autoware au chargement de la carte), et non avec les
    fractions figees sur la geometrie d'origine utilisees par l'optimisation.
    Renvoie (k max, nb de points en exces, |dk/ds| max)."""
    cl = Centrelignes({lid: bornes[lid] for lid in sorted(selection | contexte)},
                      xy, pas)
    C = cl.M @ xy(cl.ids)
    ts = cl.triplets(max(1, int(round(base / pas))))
    ts = ts[[bool(cl.proprio[b] & selection) for b in ts[:, 1]]]
    k = menger(C[ts[:, 0]], C[ts[:, 1]], C[ts[:, 2]])
    return (k.max(initial=0.0), int((k > k_max).sum()),
            _taux_courbure(C, ts, cl.suivants).max(initial=0.0))


def corriger_courbure_aoi(cl, P0, mobile, selection, contexte, k_max,
                          ecart_max, pas, base, empattement, fenetre_lat,
                          corriger, jrn):
    """Correction laterale des centrelignes en deplacant les noeuds mobiles
    des bornes, en deux temps, sous un budget d'ecart `ecart_max` commun,
    mesure depuis la geometrie d'origine : c'est une borne stricte de l'ecart
    entre la carte ecrite et la carte d'origine, independante du nombre
    d'etapes internes. (Un budget par etape a ete essaye : son seul gain
    venait d'un deplacement total plus grand, que l'on obtient aussi bien,
    garanti, en augmentant ecart_max.)

    1. lissage lateral global (fenetre_lat > 0), pour une courbure a
       variation continue :
           min  lam sum |N.D3 C|^2 + lam mu sum |L - L0|^2 + sum a_i |X_i - X0_i|^2
       N.D3 C : derivee troisieme de la centreligne projetee sur la normale
       (differences divisees, cf. _operateur_normal), soit la derivee de la
       courbure k' : on penalise les sauts de courbure et non la courbure, et
       aucune force ne pousse les noeuds le long de la route.
       lam = (fenetre / 2 pi pas)^6 place la coupure du filtre a une longueur
       d'onde `fenetre_lat`. Budget d'ecart et non-aggravation de la courbure
       sont tenus noeud par noeud, par renforcement de l'ancrage a_i.
    2. correction locale des depassements de k_max (corriger=True), meme
       formulation que pour la trajectoire, en difference seconde :
           min  sum lam_i |D2 C|^2 + mu sum lam_i |L - L0|^2 + sum |X - X0|^2
       grappe par grappe de violations, lam par bissection.

    Les variables sont les noeuds X, pas les echantillons C = M P.
    L = M_l P est le vecteur gauche - droite : sans ce terme une borne
    partagee par deux voies opposees se deforme pour servir les deux
    centrelignes et la largeur de voie derive.
    Renvoie (P corrige, residus, k avant, k apres)."""
    jrn.titre("Etape A2 : lissage lateral et correction de courbure")
    s = max(1, int(round(base / pas)))
    ts = cl.triplets(s)
    # seuls les points de centreligne portes par un lanelet selectionne sont
    # controles ; les voisins hors AOI ne servent que de conditions aux limites
    ts = ts[[bool(cl.proprio[b] & selection) for b in ts[:, 1]]]
    jrn("  courbure max admissible : %.4f 1/m  (rayon min %.2f m)"
        % (k_max, 1.0 / k_max))
    jrn("  ecart lateral autorise  : %.2f m depuis la geometrie d'origine "
        "(commun aux deux passes)" % ecart_max)
    jrn("  centrelignes            : %d echantillons au pas %.2f m"
        % (cl.M.shape[0], pas))
    jrn("  base de courbure        : +/- %.2f m (%d echantillons)" % (s * pas, s))

    im, jf = np.flatnonzero(mobile), np.flatnonzero(~mobile)
    n_ech = cl.M.shape[0]
    Al = cl.M_l[:, im].tocsr()
    Bl = cl.M_l[:, jf] @ P0[jf] - cl.M_l @ P0
    X_orig = P0[im].copy()

    def centres(X):
        P = P0.copy()
        P[im] = X
        return cl.M @ P

    def k_de(C, t):
        return menger(C[t[:, 0]], C[t[:, 1]], C[t[:, 2]])

    def ecart(Xc):
        return np.hypot(*(Xc - X_orig).T).max(initial=0.0)

    nm = len(im)
    Al2 = sp.block_diag([Al, Al]).tocsr()
    Bl2 = np.concatenate([Bl[:, 0], Bl[:, 1]])

    def solveur(Ox, Oy, W, Wl, X_dep, ancre=None):
        """X(lam) minimisant
            lam (|Ox Px + Oy Py|_W^2 + |L - L0|_Wl^2) + sum ancre_i |X_i - X_dep_i|^2,
        P = noeuds (mobiles X, figes P0). Ox, Oy couplent x et y (projection
        sur la normale) : inconnue unique [X_x ; X_y]."""
        A = sp.hstack([Ox[:, im], Oy[:, im]]).tocsr()
        B = Ox[:, jf] @ P0[jf, 0] + Oy[:, jf] @ P0[jf, 1]
        Wl2 = sp.block_diag([Wl, Wl])
        H = (A.T @ W @ A + Al2.T @ Wl2 @ Al2).tocsc()
        g = A.T @ (W @ B) + Al2.T @ (Wl2 @ Bl2)
        a = np.tile(np.ones(nm) if ancre is None else ancre, 2)
        x0 = np.concatenate([X_dep[:, 0], X_dep[:, 1]])

        def lisser(lam):
            z = spl.splu((sp.diags(a) + lam * H).tocsc()).solve(a * x0 - lam * g)
            return np.column_stack([z[:nm], z[nm:]])
        return lisser

    k_avant = k_de(centres(X_orig), ts)
    taux_avant = _taux_courbure(centres(X_orig), ts, cl.suivants)
    X, residus = X_orig.copy(), []

    # raccords : echantillons communs a un lanelet traite et a un voisin fige
    raccords = np.array([v for v, p in enumerate(cl.proprio)
                         if p & selection and p & contexte], dtype=int)
    t1 = cl.triplets(1)
    t1_racc = t1[np.isin(t1[:, 1], raccords)]

    def bilan_raccords(C):
        cas = _cassures(C, t1_racc)
        taux = _taux_courbure(C, ts, cl.suivants, raccords)
        return cas.max(initial=0.0), taux.max(initial=0.0)
    racc_avant = bilan_raccords(centres(X_orig))

    # 1. lissage lateral global -------------------------------------------
    if fenetre_lat > 0 and len(im):
        lam = (fenetre_lat / (2.0 * math.pi * pas)) ** 6
        chem4 = cl.chemins(4)
        Nx3, Ny3 = _operateur_normal(centres(X_orig), chem4, 3, pas, n_ech)
        O3 = (Nx3 @ cl.M, Ny3 @ cl.M)
        # raccord lisse avec les voisins figes : les differences troisiemes a
        # cheval sur un lanelet traite et un voisin (continuite de position,
        # de tangente et de courbure) sont renforcees
        sel_v = np.array([bool(p & selection) for p in cl.proprio])
        ctx_v = np.array([bool(p & contexte) for p in cl.proprio])
        a_cheval = sel_v[chem4].any(axis=1) & ctx_v[chem4].any(axis=1)
        W3 = sp.diags(np.where(a_cheval, POIDS_RACCORD, 1.0))
        Wl3 = POIDS_LARGEUR_LISSAGE * sp.identity(n_ech)
        # Contraintes tenues noeud par noeud : l'ancrage des seuls noeuds
        # fautifs est renforce, iterativement, pour que le lissage garde sa
        # force partout ailleurs (une reduction globale de lambda le
        # sacrifierait sur toute l'AOI pour un seul endroit contraint).
        #  - budget d'ecart : noeuds deplaces de plus de ecart_max ;
        #  - courbure : le lissage de la derivee de courbure etale les pics
        #    et peut en remonter d'autres. Un point de centreligne ne doit
        #    pas finir au-dessus de k_max et de sa courbure d'origine : les
        #    noeuds dont il depend sont ancres.
        pos_im = np.full(cl.M.shape[1], -1)
        pos_im[im] = np.arange(len(im))

        def noeuds_de(echantillons):
            cols = cl.M[np.unique(echantillons)].indices
            p = pos_im[cols]
            return np.unique(p[p >= 0])

        def degrades(Xc):
            k = k_de(centres(Xc), ts)
            return (k > k_max) & (k > k_avant + 1e-9)

        ancre = np.ones(len(im))
        retenus_ecart = np.zeros(len(im), dtype=bool)
        retenus_courbure = np.zeros(len(im), dtype=bool)
        for _ in range(N_ANCRAGE):
            X = solveur(*O3, W3, Wl3, X_orig, ancre)(lam)
            d = np.hypot(*(X - X_orig).T)
            trop = d > ecart_max
            fautifs = noeuds_de(ts[degrades(X)].ravel())
            if not trop.any() and not len(fautifs):
                break
            ancre[trop] *= 2.0 * (d[trop] / ecart_max) ** 2
            ancre[fautifs] *= 4.0
            retenus_ecart |= trop
            retenus_courbure[fautifs] = True

        def admissible(Xc):
            return ecart(Xc) <= ecart_max and not degrades(Xc).any()
        reduit = not admissible(X)
        if reduit:
            # repli si l'ancrage local n'a pas converge : reduction globale de
            # lambda, ancrage conserve
            X = bissection(solveur(*O3, W3, Wl3, X_orig, ancre),
                           admissible, -1, hi=lam)
            X = X_orig.copy() if X is None else X
        retenus = []
        if retenus_ecart.any():
            retenus.append("%d par le budget d'ecart" % retenus_ecart.sum())
        if retenus_courbure.any():
            retenus.append("%d pour ne pas aggraver la courbure"
                           % retenus_courbure.sum())
        jrn("  lissage lateral         : fenetre %.1f m%s%s"
            % (fenetre_lat,
               ", noeuds ancres : %s" % ", ".join(retenus) if retenus else "",
               "\n    /!\\ ancrage local non converge : lissage attenue sur "
               "toute la zone" if reduit else ""))
        C1 = centres(X)
        jrn("    ecart max %.3f m ; |dk/ds| max %.4f -> %.4f 1/m2, "
            "p95 %.4f -> %.4f 1/m2"
            % (ecart(X), taux_avant.max(initial=0.0),
               _taux_courbure(C1, ts, cl.suivants).max(initial=0.0),
               np.percentile(taux_avant, 95) if len(taux_avant) else 0.0,
               np.percentile(_taux_courbure(C1, ts, cl.suivants), 95)
               if len(taux_avant) else 0.0))
    elif fenetre_lat <= 0:
        jrn("  lissage lateral desactive")

    # 2. correction locale des depassements --------------------------------
    C0 = centres(X)
    k_mi = k_de(C0, ts)
    viol = np.flatnonzero(k_mi > k_max)
    marge = 3.0 / k_max
    # difference seconde par coordonnee : lignes x puis lignes y
    D2 = _operateur(t1, [1.0, -2.0, 1.0], n_ech) @ cl.M
    Z2 = sp.csr_matrix(D2.shape)
    O2 = (sp.vstack([D2, Z2]).tocsr(), sp.vstack([Z2, D2]).tocsr())
    grappes = _grappes(C0[ts[viol, 1]], 2.0 * marge) if corriger and len(viol) else []
    if not corriger:
        jrn("  correction des depassements desactivee")

    def cloche(d):
        # champ de poids en cosinus sureleve autour des violations
        return 0.5 * (1.0 + np.cos(math.pi * np.minimum(d / marge, 1.0)))

    X1 = X.copy()                   # resultat du lissage lateral
    ecart_corrections = 0.0
    for grappe in grappes:
        arbre = cKDTree(C0[ts[viol[grappe], 1]])
        w2 = cloche(arbre.query(C0[t1[:, 1]])[0])
        W = sp.diags(np.concatenate([w2, w2]))
        Wl = sp.diags(POIDS_LARGEUR * cloche(arbre.query(C0)[0]))
        zone = ts[arbre.query(C0[ts[:, 1]])[0] <= marge]
        X_dep = X.copy()
        lisser = solveur(*O2, W, Wl, X_dep)

        def k_zone(Xc):
            return k_de(centres(Xc), zone).max(initial=0.0)

        if len(im):
            # budget commun, mesure depuis la geometrie d'origine
            sol = bissection(lisser, lambda Xc: k_zone(Xc) <= k_max, +1)
            if sol is not None and ecart(sol) <= ecart_max:
                X = sol
            else:
                X = lissage_admissible(lisser, X_dep, k_zone, ecart, ecart_max)
            ecart_corrections = max(ecart_corrections,
                                    np.hypot(*(X - X_dep).T).max(initial=0.0))
        kz = k_de(centres(X), zone)
        if kz.max(initial=0.0) > k_max:
            b = zone[np.argmax(kz), 1]
            r = residu(*centres(X)[b], kz.max(), empattement)
            r["lanelets"] = sorted(cl.proprio[b])
            residus.append(r)
    if not corriger:
        for i in viol:
            residus.append(dict(residu(*C0[ts[i, 1]], k_mi[i], empattement),
                                lanelets=sorted(cl.proprio[ts[i, 1]])))

    P = P0.copy()
    P[im] = X
    C = cl.M @ P
    k_apres = k_de(C, ts)
    taux_apres = _taux_courbure(C, ts, cl.suivants)
    ecarts = np.hypot(*(X - X_orig).T)
    bouges = ecarts > 1e-6
    jrn("")
    jrn("  bilan (origine -> final)")
    jrn("  points en exces de courbure : %d -> %d  (k max %.4f -> %.4f 1/m)"
        % (int((k_avant > k_max).sum()), int((k_apres > k_max).sum()),
           k_avant.max(initial=0.0), k_apres.max(initial=0.0)))
    jrn("  variation de courbure       : |dk/ds| max %.4f -> %.4f 1/m2"
        % (taux_avant.max(initial=0.0), taux_apres.max(initial=0.0)))
    if len(raccords):
        racc_apres = bilan_raccords(C)
        jrn("  raccords aux voisins figes  : %d, cassure max %.2f -> %.2f deg, "
            "|dk/ds| max %.4f -> %.4f 1/m2"
            % (len(raccords), racc_avant[0], racc_apres[0], racc_avant[1],
               racc_apres[1]))
    jrn("  noeuds deplaces             : %d / %d mobiles, moyen %.3f m, max %.3f m"
        % (bouges.sum(), len(im), ecarts[bouges].mean() if bouges.any() else 0.0,
           ecarts.max(initial=0.0)))
    jrn("  deplacement max par etape   : lissage %.3f m, corrections %.3f m ; "
        "cumule %.3f m (budget commun %.2f m)"
        % (np.hypot(*(X1 - X_orig).T).max(initial=0.0), ecart_corrections,
           ecarts.max(initial=0.0), ecart_max))
    journal_residus(jrn, residus)
    return P, residus, k_avant, k_apres


def lisser_altitude_graphe(aretes, P, z0, mobile, fenetre):
    """Lissage de l'altitude le long des bornes (jamais en travers : le devers
    est conserve). Minimise
        sum m_i (z_i - z0_i)^2 + l^2 sum_(ij) (z_i - z_j)^2 / d_ij,
    discretisation de  int (z - z0)^2 + l^2 z'^2 ds  dont le noyau
    exponentiel a la variance d'une moyenne glissante de largeur `fenetre`
    pour l = fenetre / sqrt(24). Les noeuds non mobiles sont des conditions
    aux limites : l'altitude reste continue avec le reste de la carte."""
    n = len(P)
    i, j = aretes.T
    d = np.maximum(np.hypot(*(P[i] - P[j]).T), 0.05)
    masse = np.maximum(np.bincount(i, d / 2, n) + np.bincount(j, d / 2, n), 1e-6)
    w = 1.0 / d
    lap = sp.coo_matrix((np.concatenate([w, w, -w, -w]),
                         (np.concatenate([i, j, i, j]), np.concatenate([i, j, j, i]))),
                        shape=(n, n)).tocsr()
    l2 = fenetre ** 2 / 24.0
    K = (sp.diags(masse) + l2 * lap).tocsr()
    im, jf = np.flatnonzero(mobile), np.flatnonzero(~mobile)
    z = z0.copy()
    if len(im):
        rhs = masse[im] * z0[im] - l2 * (lap[im][:, jf] @ z0[jf])
        z[im] = spl.spsolve(K[im][:, im].tocsc(), rhs)
    return z


def jacobien_geo(noeuds, ids):
    """d(lat, lon)/d(local_x, local_y), ajuste par moindres carres sur les
    noeuds `ids`. Independant de la projection (MGRS, UTM, locale) ; pour des
    deplacements de quelques decimetres l'erreur de linearisation est
    negligeable. Une seconde passe ecarte les noeuds dont lat/lon ne
    concordent pas avec local_x/local_y (edition manuelle, par exemple).
    Renvoie (J 2x2, residu max en m, nb ecartes) ou (None, None, 0)."""
    ids = set(int(i) for i in ids)
    rows = []
    for nd in noeuds:
        if int(nd.get("id")) not in ids:
            continue
        tags = {t.get("k"): t.get("v") for t in nd.findall("tag")}
        try:
            rows.append((float(tags["local_x"]), float(tags["local_y"]),
                         float(nd.get("lat")), float(nd.get("lon"))))
        except (KeyError, TypeError, ValueError):
            continue
    r = np.asarray(rows).reshape(-1, 4)
    garde = np.ones(len(r), dtype=bool)
    for _ in range(2):
        if garde.sum() < 3:
            return None, None, 0
        centre = r[garde, :2].mean(axis=0)
        X = np.column_stack([np.ones(len(r)), r[:, :2] - centre])
        if np.linalg.matrix_rank(X[garde]) < 3:
            return None, None, 0
        coef = np.linalg.lstsq(X[garde], r[garde, 2:], rcond=None)[0]
        J = coef[1:].T                  # [[dlat/dx, dlat/dy], [dlon/dx, dlon/dy]]
        res_m = np.hypot(*np.linalg.solve(J, (r[:, 2:] - X @ coef).T))
        garde = res_m <= max(0.05, 5.0 * np.median(res_m))
    return J, float(res_m[garde].max()), int((~garde).sum())


def _angles_sommets(pts, avant=None, apres=None):
    """Angle de deviation (rad) en chaque sommet d'une polyligne ; aux
    extremites, calcule avec le point qui precede / suit sur la voie voisine
    s'il est connu, nul sinon."""
    p = np.vstack([q for q in (avant, *pts, apres) if q is not None])
    v = np.diff(p, axis=0)
    th = np.arctan2(v[:, 1], v[:, 0])
    dev = np.abs((np.diff(th) + np.pi) % (2 * np.pi) - np.pi)
    ang = np.zeros(len(pts))
    deb = 0 if avant is not None else 1
    ang[deb:deb + len(dev)] = dev
    return ang


def densifier_bornes(arbre, noeuds, ways, voies, xy_de, z_de, bornes,
                     selection, lanelets, debuts, fins, n_max, long_min,
                     precision):
    """Insere au plus `n_max` noeuds par borne des lanelets selectionnes,
    pour donner au lissage la liberte de former une courbe. Chaque noeud est
    place au milieu du segment de plus fort score  longueur x (deviation aux
    deux extremites), les deviations aux bouts de borne etant mesurees avec
    la voie precedente / suivante : on densifie les virages et les raccords,
    jamais les lignes droites. Un segment n'est coupe que s'il mesure au
    moins 2 * long_min. Seules les bornes utilisees exclusivement par des
    lanelets selectionnes sont modifiees (une borne partagee avec un lanelet
    hors AOI changerait aussi celui-ci).
    Met a jour l'arbre XML, `noeuds`, `ways`, `xy_de`, `z_de`, `bornes`.
    Renvoie (ids inseres, bornes densifiees, bornes partagees ignorees)."""
    usage = {}
    for lid, g, d, st in lanelets:
        usage.setdefault(g, set()).add(lid)
        usage.setdefault(d, set()).add(lid)
    racine = arbre.getroot()
    elem_noeud = {int(nd.get("id")): nd for nd in noeuds}
    elem_way = {int(w.get("id")): w for w in racine.findall("way")}
    prochain = 1 + max(int(e.get("id")) for e in racine if e.get("id") is not None)
    pos_noeud = max(i for i, e in enumerate(racine) if e.tag == "node") + 1

    # extremites de chaque borne dans son sens de stockage
    voisins = {}
    for lid in sorted(selection):
        for cote in (0, 1):
            wid = voies[lid][cote]
            o = bornes[lid][cote]
            g, d = bornes[lid]
            av = next((bornes[p][cote][-2] for p in fins.get((g[0], d[0]), [])
                       if len(bornes[p][cote]) > 1), None)
            ap = next((bornes[q][cote][1] for q in debuts.get((g[-1], d[-1]), [])
                       if len(bornes[q][cote]) > 1), None)
            if o[0] != ways[wid][0]:
                av, ap = ap, av
            voisins.setdefault(wid, (av, ap))

    inseres, n_bornes, partagees = set(), 0, 0
    for wid in sorted(voisins):
        if not usage.get(wid, set()) <= selection:
            partagees += 1
            continue
        refs = list(ways[wid])
        av, ap = voisins[wid]
        pts = [xy_de[i] for i in refs]
        ang = list(_angles_sommets(np.array(pts),
                                   None if av is None else xy_de[int(av)],
                                   None if ap is None else xy_de[int(ap)]))
        nouveaux = 0
        for _ in range(n_max):
            L = [math.dist(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
            score = [L[i] * (ang[i] + ang[i + 1]) if L[i] >= 2 * long_min else 0.0
                     for i in range(len(L))]
            i = int(np.argmax(score))
            if score[i] <= 0.0:
                break
            a_id, b_id = refs[i], refs[i + 1]
            m = 0.5 * (np.asarray(pts[i]) + np.asarray(pts[i + 1]))
            zm = 0.5 * (z_de[a_id] + z_de[b_id])
            src_a, src_b = elem_noeud[a_id], elem_noeud[b_id]
            e = copy.deepcopy(src_a)
            e.set("id", str(prochain))
            for att in ("lat", "lon"):
                if src_a.get(att) is not None and src_b.get(att) is not None:
                    e.set(att, "%.11f" % (0.5 * (float(src_a.get(att))
                                                 + float(src_b.get(att)))))
            for t in e.findall("tag"):
                v = {"local_x": m[0], "local_y": m[1], "ele": zm}.get(t.get("k"))
                if v is not None:
                    t.set("v", "%.*f" % (precision, v))
            racine.insert(pos_noeud, e)
            pos_noeud += 1
            noeuds.append(e)
            elem_noeud[prochain] = e
            xy_de[prochain], z_de[prochain] = m, zm
            refs.insert(i + 1, prochain)
            pts.insert(i + 1, m)
            ang.insert(i + 1, 0.0)
            inseres.add(prochain)
            prochain += 1
            nouveaux += 1
        if not nouveaux:
            continue
        _reecrire_nds(elem_way[wid], refs)
        ways[wid] = refs
        n_bornes += 1

    for lid in selection:
        g, d = voies[lid]
        bornes[lid] = orienter_bornes(np.array(ways[g]), np.array(ways[d]),
                                      lambda ids: np.array([xy_de[int(i)] for i in ids]))
    return inseres, n_bornes, partagees


def _reecrire_nds(w, refs):
    """Remplace la liste des <nd> du way `w`, mise en forme conservee."""
    nds = w.findall("nd")
    queue, derniere = nds[0].tail, nds[-1].tail
    for nd in nds:
        w.remove(nd)
    for k, ref in enumerate(refs):
        nd = ET.Element("nd", ref=str(ref))
        nd.tail = derniere if k == len(refs) - 1 else queue
        w.insert(k, nd)


def _ecart_corde(pts, p, q):
    """Ecarts (lateral, altitude) des points pts[p+1:q] a la corde pts[p] ->
    pts[q] ; l'altitude de reference est interpolee le long de la corde."""
    a, b, m = pts[p], pts[q], pts[p + 1:q]
    ab = b[:2] - a[:2]
    L2 = float(ab @ ab)
    t = np.zeros(len(m)) if L2 == 0 else np.clip((m[:, :2] - a[:2]) @ ab / L2, 0, 1)
    proj = a + t[:, None] * (b - a)
    return np.hypot(*(m[:, :2] - proj[:, :2]).T), np.abs(m[:, 2] - proj[:, 2])


def elaguer_noeuds_inseres(arbre, noeuds, ways, inseres, position, tol):
    """Retire les noeuds inseres restes inutiles : ceux qui, apres lissage,
    sont encore alignes (a `tol` pres, en plan et en altitude) avec leurs
    voisins conserves. Retrait glouton, du plus aligne au moins aligne ; un
    candidat n'est retire que si TOUS les points deja retires entre ses
    voisins restent a moins de `tol` de la nouvelle corde : l'ecart cumule
    entre la borne lissee et la borne elaguee est donc borne par `tol`.
    `position` : {id: (x, y, z)} final. Met a jour l'arbre XML, `noeuds` et
    `ways`. Renvoie (ids retires, ecart max introduit en m)."""
    racine = arbre.getroot()
    elem_way = {int(w.get("id")): w for w in racine.findall("way")}
    retires, ecart_max = set(), 0.0
    for wid, refs in ways.items():
        if not inseres.intersection(refs):
            continue
        pts = np.array([position[i] for i in refs], dtype=float)
        garde = np.ones(len(refs), dtype=bool)
        while True:
            idx = np.flatnonzero(garde)
            meilleur = None
            for j in range(1, len(idx) - 1):
                k = idx[j]
                if refs[k] not in inseres:
                    continue
                dl, dz = _ecart_corde(pts, idx[j - 1], idx[j + 1])
                e = max(dl.max(initial=0.0), dz.max(initial=0.0))
                if e <= tol and (meilleur is None or e < meilleur[0]):
                    meilleur = (e, k)
            if meilleur is None:
                break
            garde[meilleur[1]] = False
            ecart_max = max(ecart_max, meilleur[0])
        if garde.all():
            continue
        retires.update(refs[k] for k in np.flatnonzero(~garde))
        refs = [r for r, g in zip(refs, garde) if g]
        _reecrire_nds(elem_way[wid], refs)
        ways[wid] = refs
    if retires:
        for nd in [nd for nd in noeuds if int(nd.get("id")) in retires]:
            racine.remove(nd)
        noeuds[:] = [nd for nd in noeuds if int(nd.get("id")) not in retires]
    return retires, ecart_max
