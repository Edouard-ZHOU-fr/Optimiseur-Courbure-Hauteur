# -*- coding: utf-8 -*-
"""Lecture et validation de config/config.yaml.

Le fichier est la seule source des parametres : toutes les cles sont
obligatoires, une cle inconnue est refusee (faute de frappe) et les types
sont verifies. Le resultat expose les memes attributs plats qu'un
argparse.Namespace (a.carte, a.ecart_max, ...)."""

from pathlib import Path
from types import SimpleNamespace

try:
    import yaml
except ImportError:
    yaml = None

RACINE = Path(__file__).resolve().parent.parent
CONFIG_DEFAUT = RACINE / "config" / "config.yaml"

# section -> {cle: (type, nullable)} ; type "chemin" = Path relatif a RACINE
SCHEMA = {
    None: {"aoi": (bool, False)},
    # ---- communs aux deux modes
    "entrees_sorties": {
        "carte": ("chemin", False),
        "sortie_carte": ("chemin", True),
        "rapport": ("chemin", True),
        "graphique": (bool, False),
        "silencieux": (bool, False),
    },
    "vehicule": {
        "fichier": ("chemin", False),
        "modele": (str, True),
        "lister_vehicules": (bool, False),
        "angle_braquage_max": (float, True),
        "empattement": (float, True),
        "marge_courbure": (float, True),
    },
    "courbure": {
        "ecart_max": (float, False),
        "pas_reech": (float, False),
        "sans_courbure": (bool, False),
    },
    "altitude": {
        "fenetre_z": (float, False),
        "sans_lissage_z": (bool, False),
        "precision": (int, False),
    },
    # ---- mode AOI seulement
    "mode_aoi": {
        "base_courbure": (float, True),
        "fenetre_laterale": (float, False),
        "sans_lissage_lateral": (bool, False),
        "noeuds_ajoutes_max": (int, False),
    },
    # ---- mode trajectoire seulement (lu mais sans effet en mode AOI)
    "mode_trajectoire": {
        "traj": ("chemin", False),
        "sortie_traj": ("chemin", True),
        "pas_min": (float, False),
        "seuil_rebroussement": (float, False),
        "distance_max": (float, False),
        "ecart_min_z": (float, False),
        "offset_capteur": (float, True),
        "ele_nul_mediane": (bool, False),
    },
}

# nom de la cle YAML -> nom de l'attribut attendu par le reste du code
RENOMMAGES = {("vehicule", "fichier"): "vehicule"}


def _convertir(val, typ, nullable, ou):
    if val is None:
        if nullable:
            return None
        raise ValueError("%s : valeur obligatoire (null interdit)" % ou)
    if typ == "chemin":
        if not isinstance(val, str) or not val.strip():
            raise ValueError("%s : chemin attendu, lu %r" % (ou, val))
        p = Path(val).expanduser()
        return p if p.is_absolute() else RACINE / p
    if typ is bool:
        if not isinstance(val, bool):
            raise ValueError("%s : true ou false attendu, lu %r" % (ou, val))
        return val
    if typ is int:
        if isinstance(val, bool) or not isinstance(val, int):
            raise ValueError("%s : entier attendu, lu %r" % (ou, val))
        return val
    if typ is float:
        if isinstance(val, bool) or not isinstance(val, (int, float)):
            raise ValueError("%s : nombre attendu, lu %r" % (ou, val))
        return float(val)
    if not isinstance(val, str):
        raise ValueError("%s : texte attendu, lu %r" % (ou, val))
    return val


def charger_config(chemin=CONFIG_DEFAUT):
    chemin = Path(chemin)
    if yaml is None:
        raise SystemExit("PyYAML est requis pour lire %s (pip install pyyaml)"
                         % chemin)
    if not chemin.is_file():
        raise SystemExit("Fichier de configuration introuvable : %s" % chemin)
    try:
        doc = yaml.safe_load(chemin.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as err:
        raise SystemExit("%s : YAML invalide\n%s" % (chemin, err))
    if not isinstance(doc, dict):
        raise SystemExit("%s : dictionnaire attendu a la racine" % chemin)

    attendues = set(SCHEMA[None]) | {s for s in SCHEMA if s}
    inconnues = sorted(set(doc) - attendues)
    erreurs = ["cle inconnue : %s" % k for k in inconnues]
    params = {}
    for section, cles in SCHEMA.items():
        bloc = doc if section is None else doc.get(section)
        if section is not None:
            if bloc is None:
                erreurs.append("section manquante : %s" % section)
                continue
            if not isinstance(bloc, dict):
                erreurs.append("%s : section attendue" % section)
                continue
            erreurs += ["cle inconnue : %s.%s" % (section, k)
                        for k in sorted(set(bloc) - set(cles))]
        for cle, (typ, nullable) in cles.items():
            ou = cle if section is None else "%s.%s" % (section, cle)
            if cle not in bloc:
                erreurs.append("cle manquante : %s" % ou)
                continue
            try:
                params[RENOMMAGES.get((section, cle), cle)] = _convertir(
                    bloc[cle], typ, nullable, ou)
            except ValueError as err:
                erreurs.append(str(err))
    if erreurs:
        raise SystemExit("%s :\n  %s" % (chemin, "\n  ".join(erreurs)))
    params["config"] = chemin
    return SimpleNamespace(**params)
