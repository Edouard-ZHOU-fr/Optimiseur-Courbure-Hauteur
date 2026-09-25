# -*- coding: utf-8 -*-
"""Modele de vehicule : lecture du YAML projet ou Autoware, limite
cinematique de courbure (modele bicyclette)."""

import math
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


# fichier de modeles, dans le dossier config de l'outil
VEHICULE_DEFAUT = Path(__file__).resolve().parent.parent / "config" / "vehicle_model.yaml"


DEFAUTS_VEHICULE = {
    "empattement": 2.800,
    "angle_braquage_max": 30.0,
    "hauteur_capteur": 0.26,
    "marge_courbure": 1.0,
}


def _depuis_autoware(doc):
    """Lit un vehicle_info.param.yaml Autoware ; max_steer_angle est en rad."""
    for cle, val in doc.items():
        if not isinstance(val, dict):
            continue
        prm = val.get("ros__parameters")
        if isinstance(prm, dict) and "wheel_base" in prm:
            v = {"empattement": float(prm["wheel_base"])}
            if "max_steer_angle" in prm:
                v["angle_braquage_max"] = math.degrees(float(prm["max_steer_angle"]))
            return v
    return None


def charger_vehicule(chemin, nom_modele, jrn):
    """Renvoie (parametres, etiquette). Format natif du projet ou format
    Autoware vehicle_info.param.yaml."""
    params = dict(DEFAUTS_VEHICULE)
    if chemin is not None and not Path(chemin).is_file() and chemin != VEHICULE_DEFAUT:
        indice = ""
        if VEHICULE_DEFAUT.is_file() and yaml is not None:
            doc = yaml.safe_load(VEHICULE_DEFAUT.read_text(encoding="utf-8")) or {}
            if str(chemin) in (doc.get("vehicules") or {}):
                indice = ("\n'%s' est un nom de modele : le mettre dans "
                          "vehicule.modele de config.yaml" % chemin)
        raise SystemExit("Fichier vehicule introuvable : %s%s" % (chemin, indice))
    if chemin is None or not Path(chemin).is_file():
        if nom_modele:
            raise SystemExit("vehicule.modele %s demande mais aucun fichier vehicule "
                             "lisible (%s)" % (nom_modele, chemin))
        return params, "valeurs par defaut integrees"
    if yaml is None:
        raise SystemExit("PyYAML est requis pour lire %s (pip install pyyaml)"
                         % chemin)

    doc = yaml.safe_load(Path(chemin).read_text(encoding="utf-8")) or {}

    aw = _depuis_autoware(doc)
    if aw is not None:
        if nom_modele:
            jrn("  (vehicule.modele ignore : %s est un fichier Autoware mono-vehicule)"
                % Path(chemin).name)
        params.update(aw)
        return params, "%s [format Autoware]" % Path(chemin).name

    vehicules = doc.get("vehicules")
    if not isinstance(vehicules, dict) or not vehicules:
        raise SystemExit("%s : ni cle 'vehicules', ni format Autoware reconnu."
                         % chemin)

    nom = nom_modele or doc.get("modele_par_defaut") or next(iter(vehicules))
    if nom not in vehicules:
        raise SystemExit("Modele '%s' absent de %s. Disponibles : %s"
                         % (nom, chemin, ", ".join(sorted(vehicules))))
    fiche = vehicules[nom] or {}
    for cle in DEFAUTS_VEHICULE:
        if cle in fiche:
            params[cle] = float(fiche[cle])
    etiquette = nom
    if fiche.get("description"):
        etiquette += " (%s)" % fiche["description"]
    return params, etiquette


def lister_vehicules(chemin):
    if yaml is None or chemin is None or not Path(chemin).is_file():
        print("Aucun fichier vehicule lisible : %s" % chemin)
        return 1
    doc = yaml.safe_load(Path(chemin).read_text(encoding="utf-8")) or {}
    veh = doc.get("vehicules")
    if not isinstance(veh, dict):
        aw = _depuis_autoware(doc)
        if aw:
            print("%s : fichier Autoware mono-vehicule" % chemin)
            print("  empattement=%.3f m  braquage max=%.2f deg"
                  % (aw["empattement"], aw.get("angle_braquage_max", float("nan"))))
            return 0
        print("%s : aucun vehicule declare." % chemin)
        return 1
    defaut = doc.get("modele_par_defaut")
    print("Vehicules declares dans %s :\n" % chemin)
    print("  %-22s %-11s %-11s %-10s %s"
          % ("modele", "empatt.(m)", "braquage", "R_min(m)", "description"))
    for nom in sorted(veh):
        f = veh[nom] or {}
        emp = float(f.get("empattement", DEFAUTS_VEHICULE["empattement"]))
        ang = float(f.get("angle_braquage_max",
                          DEFAUTS_VEHICULE["angle_braquage_max"]))
        r = emp / math.tan(math.radians(ang))
        print("  %-22s %-11.3f %-11s %-10.2f %s%s"
              % (nom, emp, "%.2f deg" % ang, r, f.get("description", ""),
                 "   <- defaut" if nom == defaut else ""))
    return 0


def courbure_max_vehicule(angle_braquage_deg, empattement, marge=1.0):
    """Modele bicyclette : R_min = L / tan(delta_max), k_max = 1 / R_min.
    `marge` < 1 durcit la contrainte (garde cinematique)."""
    delta = math.radians(angle_braquage_deg)
    if not (0.0 < delta < math.pi / 2):
        raise ValueError("angle de braquage hors domaine (0, 90) degres")
    return marge * math.tan(delta) / empattement
