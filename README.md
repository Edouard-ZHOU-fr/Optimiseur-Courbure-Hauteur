# Optimiseur d'altitude Lanelet2 avec correction de courbure

Recale l'altitude (`ele`) des noeuds d'une carte Lanelet2 a partir d'une
trajectoire reellement parcourue, et verifie que le trace reste compatible
avec la cinematique du vehicule (rayon de braquage minimal).

Un second mode, `aoi: true`, se passe de trajectoire : il corrige directement la
courbure des lanelets situes dans une zone AOI dessinee sur la carte et lisse
leur altitude (cf. [Mode AOI](#mode-aoi)).

### Contenu du dossier

| Fichier | Role |
|---|---|
| `optimiseur_hauteur_courbure.py` | Script principal : arguments, enchainement des etapes (`traiter_trajectoire`, `traiter_aoi`, `main`), constante `AOI`. |
| `utils/` | Fonctions, par theme (cf. ci-dessous). |
| `config/config.yaml` | **Tous les parametres** d'execution, relus a chaque lancement. |
| `config/vehicle_model.yaml` | Parametres des vehicules du parc. |
| `traj.csv`, `exemples/Test.osm` | Jeu de donnees d'exemple du mode trajectoire. |
| `log/` | Rapports horodates. |

Le dossier est autonome : aucune dependance vers un autre repertoire du depot.

| Module `utils/` | Contenu |
|---|---|
| `journal.py` | `Journal` : affichage console et rapport texte. |
| `vehicule.py` | Lecture du YAML projet ou Autoware, `courbure_max_vehicule`. |
| `geometrie.py` | Courbure de Menger, cap, tests polygone, echantillonnage de polylignes. |
| `trajectoire.py` | Mode trajectoire : decimation, rebroussements, reechantillonnage, lissage Z. |
| `courbure.py` | Correction de courbure sous contrainte d'ecart, bissection, residus. |
| `aoi.py` | Mode AOI : topologie Lanelet2, `Centrelignes`, lissage lateral et correction, lissage Z des bornes, insertion de noeuds, jacobien lat/lon. |
| `entrees_sorties.py` | Lecture / ecriture CSV et OSM, appariement trajectoire <-> carte. |
| `figures.py` | Figures de controle PNG. |

### Rapport avec l'ancien outil

Cet outil **remplace** `../optimisseur_KNN/optimiseur_hauteur_KNN.py`, conserve
a part pour reference historique et a ne plus utiliser : il produit une
altitude decalee d'un echantillon et ne verifie pas la courbure
(cf. [Corrections](#corrections-par-rapport-a-lancienne-version)).

```
optimisation_laneletMap_axeZ/
├── optimiseur_courbure_hauteur/   <- cet outil
└── optimisseur_KNN/               <- ancienne version, figee
```

---

## Sommaire

- [Prerequis](#prerequis)
- [Prise en main](#prise-en-main)
- [Modele de vehicule](#modele-de-vehicule)
- [Chaine de traitement](#chaine-de-traitement)
- [Mode AOI](#mode-aoi)
- [Parametres](#parametres)
- [Fichiers produits](#fichiers-produits)
- [Lire le rapport](#lire-le-rapport)
- [Points de vigilance](#points-de-vigilance)
- [Compatibilite VectorMap](#compatibilite-vectormap)
- [Corrections par rapport a l'ancienne version](#corrections-par-rapport-a-lancienne-version)

---

## Prerequis

Python 3.10+ et :

```bash
pip install numpy scipy pyyaml matplotlib
```

Versions validees : Python 3.10.12, numpy 1.26.4, scipy 1.10.1, PyYAML 5.4.1,
matplotlib 3.5.1. `matplotlib` n'est necessaire qu'avec `entrees_sorties.graphique`.

### Entrees attendues

**Trajectoire** — CSV avec en-tete, colonnes `x,y,z` (une 4e colonne `yaw`
est toleree mais ignoree, elle est recalculee) :

```csv
x,y,z,yaw
48571.61386891395,12782.825104412777,62.79046812548029,0
```

Les coordonnees sont dans le meme repere local que la carte (`local_x`,
`local_y`). En pratique on exporte le topic ROS `/tf` (et `/tf_static`).

**Carte** — Lanelet2 `.osm`, chaque `node` portant les tags `local_x`,
`local_y` et `ele`.

---

## Prise en main

Le script **ne prend aucun argument** : tous les parametres sont dans
[`config/config.yaml`](config/config.yaml), relu a chaque execution.

```bash
python3 optimiseur_hauteur_courbure.py
```

Exemples de reglages dans `config/config.yaml` :

```yaml
aoi: false                      # mode trajectoire
entrees_sorties:
  carte: exemples/Test.osm
  graphique: true               # ajouter les figures de controle
mode_trajectoire:
  traj: traj.csv
vehicule:
  modele: milla_bus_l3h2        # choisir un vehicule
  lister_vehicules: true        # lister les vehicules disponibles puis quitter
```

Toutes les cles sont **obligatoires** ; une cle manquante, inconnue (faute de
frappe) ou de mauvais type est refusee, avec la liste complete des erreurs.
`null` vaut valeur par defaut (indiquee en commentaire dans le fichier). Les
chemins relatifs partent du dossier du script, quel que soit le dossier
courant. Lancer le script avec un argument est refuse, pour ne pas laisser
croire qu'il est pris en compte.

Duree typique : **~1.6 s** pour 9 100 points de trajectoire et 28 700 noeuds
de carte (mode trajectoire).

### Codes de retour

| Code | Signification |
|---|---|
| `0` | Traitement complet, aucune zone hors limite cinematique. |
| `2` | Traitement complet, **mais** des virages restent infranchissables (detail dans le rapport). |
| `1` | Erreur (fichier introuvable, format inattendu, modele inconnu...). |

Le code `2` est concu pour etre exploite en CI : la carte est ecrite, mais
elle merite une relecture avant mise en service.

---

## Modele de vehicule

La contrainte de courbure vient du modele bicyclette :

```
R_min = empattement / tan(angle_braquage_max)
k_max = marge_courbure / R_min
```

### Trois sources, par ordre de priorite

**1. Surcharges de `config.yaml`** (prioritaires sur le fichier vehicule)

```yaml
vehicule:
  empattement: 3.682
  angle_braquage_max: 34.38
```

**2. `config/vehicle_model.yaml`** — format du projet, plusieurs vehicules :

```yaml
modele_par_defaut: lohr_cristal

vehicules:
  lohr_cristal:
    description: "Lohr Cristal"
    empattement: 2.670          # [m]
    angle_braquage_max: 30.00   # [deg]
    hauteur_capteur: 0.26       # [m]
    marge_courbure: 1.0
```

**3. Fichier Autoware natif** — pointer directement le package de description
evite de dupliquer les valeurs :

```yaml
vehicule:
  fichier: .../lohr_cristal_description/config/vehicle_info.param.yaml
```

Le script y lit `wheel_base` et `max_steer_angle`. Autoware exprime l'angle en
**radians**, la conversion est automatique (0.60 rad -> 34.38 deg).

`vehicule.fichier` doit etre un chemin : y mettre un nom de modele est refuse
avec l'indication de le placer dans `vehicule.modele`.

### Parc actuellement decrit

| Modele | Empattement | Braquage max | R_min |
|---|---|---|---|
| `lohr_cristal` | 2.670 m | 30.00 deg | 4.62 m |
| `milla_car_l2h2` | 3.682 m | 34.38 deg | 5.38 m |
| `renault_master_l2h2` | 3.682 m | 34.38 deg | 5.38 m |
| `milla_car_l3h2` | 4.332 m | 34.38 deg | 6.33 m |
| `milla_bus_l3h2` | 4.332 m | 34.38 deg | 6.33 m |
| `renault_master_l3h2` | 4.332 m | 34.38 deg | 6.33 m |

`marge_courbure` durcit la contrainte sans toucher aux valeurs physiques :
`0.8` reserve 20 % de debattement de braquage au controleur.

---

## Chaine de traitement

```
CSV trajectoire ─┐
                 ├─> 2. nettoyage ─> 3. reechantillonnage ─> 4. courbure ─┐
carte .osm ──────┘                                                        │
                                          5. cap  6. lissage Z <──────────┘
                                                      │
                          7. appariement KD-tree <────┘
                                                      │
                          8. ecriture carte + trajectoire + rapport
```

### 2. Nettoyage — l'etape qui compte le plus

Deux defauts rendent la courbure brute inexploitable :

- **Points quasi immobiles.** Vehicule a l'arret, pas descendant jusqu'a
  0.4 mm. La courbure discrete y explose (division par ~0). `mode_trajectoire.pas_min`
  (0.30 m) les elimine.
- **Rebroussements.** Le cap s'inverse de 180 deg en moins d'un metre :
  marche arriere, manoeuvre en plusieurs points, ou saut de localisation.
  La courbure n'a aucun sens sur un point de rebroussement. `mode_trajectoire.seuil_rebroussement`
  (90 deg) y **coupe la trajectoire en troncons** traites separement.

Sur le jeu de donnees de reference, cette seule etape fait passer la courbure
maximale de **0.936 (R = 1.07 m) a 0.19 (R = 5.19 m)**. Les rebroussements
detectes sont listes avec leurs coordonnees dans le rapport : ils meritent
une verification, un saut de localisation signale un probleme d'enregistrement.

### 3. Reechantillonnage

Le pas brut varie de 0.4 mm a 23 m (mediane 8.65 m) : deux courbures mesurees
a des pas differents ne sont pas comparables, et un virage a R = 5 m a un arc
plus court qu'un seul intervalle d'echantillonnage. Chaque troncon est donc
rejoue par spline cubique parametree par l'abscisse curviligne, a pas constant
(`courbure.pas_reech`, 1.0 m).

Effet secondaire utile : l'appariement avec la carte devient beaucoup plus fin
(distance mediane ~1.3 m au lieu de ~8.6 m).

### 4. Correction de courbure

Pour chaque zone ou `k > k_max`, on resout

```
min  Σ λᵢ ‖pᵢ₋₁ - 2pᵢ + pᵢ₊₁‖²  +  Σ ‖pᵢ - refᵢ‖²
```

systeme bande resolu directement, `λ` cherche par **bissection** sous la
contrainte `courbure.ecart_max`.

Deux pieges evites, tous deux verifies experimentalement :

- **Un lissage laplacien pur ne convient pas.** C'est un flot de
  raccourcissement de courbe : sur un arc soutenu il retrecit le rayon, donc
  *augmente* la courbure. Le terme d'ancrage sur `ref` neutralise cet effet.
- **Aucun ecretage a posteriori.** Ramener de force les points dans le budget
  d'ecart apres resolution recree des angles vifs — et donc des violations
  supplementaires. La contrainte est dans la recherche du `λ`.

Si le budget d'ecart ne suffit pas, on retient parmi les lissages admissibles
celui qui **minimise** la courbure maximale, en partant de la zone non
modifiee : la correction ne peut jamais degrader le trace.

### 6. Lissage de l'altitude

Moyenne glissante centree a phase nulle, fenetre exprimee **en metres**
(`altitude.fenetre_z`, 25 m) et convertie en echantillons selon `courbure.pas_reech`.

### 7. Appariement

KD-tree (`scipy.spatial.cKDTree`) : chaque noeud de la carte prend l'altitude
du point de trajectoire le plus proche, moins `hauteur_capteur`. Un noeud est
modifie si son plus proche voisin est a moins de `mode_trajectoire.distance_max` **et** que
l'ecart d'altitude atteint `mode_trajectoire.ecart_min_z`.

---

## Mode AOI

```yaml
aoi: true
entrees_sorties:
  carte: lanelet2_map.osm
```

Aucune trajectoire n'est lue : toute la section `mode_trajectoire` est
ignoree. La carte
elle-meme est corrigee : **position (`local_x`, `local_y`, `lat`, `lon`) et
altitude (`ele`)** des noeuds de bornes.

### Zone AOI

Un way ferme portant `subtype=AOI`, tel que le produit VectorMap Builder :

```xml
<way id="22383">
  <nd ref="22378"/> <nd ref="22382"/> <nd ref="22381"/> <nd ref="22379"/>
  <tag k="type" v="detection_area"/>
  <tag k="subtype" v="AOI"/>
  <tag k="area" v="yes"/>
</way>
```

Plusieurs AOI peuvent coexister. Seuls les lanelets **entierement couverts**
par une AOI sont traites : tous les noeuds de leurs deux bornes a l'interieur,
et aucun cote de borne qui traverse le bord de l'AOI (AOI non convexe). Les
lanelets qui la traversent ou y debordent sont comptes dans le rapport
(`partiellement (ignores)`) et laisses intacts. Les lanelets `crosswalk`,
`walkway` et `stairs` sont ignores.

### Noeuds mobiles et noeuds figes

Un noeud de borne n'est deplace que s'il n'appartient a **aucun** lanelet
roulable hors selection. Les noeuds partages avec le reste de la carte
(jonction avec un successeur hors AOI, borne commune avec une voie voisine
non selectionnee) restent figes et servent de conditions aux limites : la
carte reste continue. Les predecesseurs et successeurs directs sont integres
au calcul, figes, pour que la continuite aux jonctions soit prise en compte
(cf. [Raccords](#raccords-avec-les-voies-hors-aoi)). Pour qu'une borne
partagee entre deux voies soit corrigee, l'AOI doit couvrir les deux voies.

### Insertion de noeuds

Les bornes sont des polylignes a sommets espaces (~3.9 m) : en ne deplacant
que les sommets existants, un virage reste une suite d'angles vifs.
`mode_aoi.noeuds_ajoutes_max: N` (0 = aucune insertion) autorise
jusqu'a **N noeuds inseres par borne** avant le lissage. Chaque noeud est
place au milieu du segment de plus fort score

```
longueur × (deviation au sommet de debut + deviation au sommet de fin)
```

les deviations aux bouts de borne etant mesurees avec la voie precedente /
suivante : on densifie les virages et les raccords, jamais les lignes
droites. Un segment n'est coupe que s'il mesure au moins `2 × courbure.pas_reech`,
ce qui plafonne naturellement le nombre de noeuds utiles (sur les zones de
test, l'effet sature vers N = 5).

Seules les bornes utilisees **exclusivement** par des lanelets traites sont
densifiees : une borne partagee avec un lanelet hors AOI changerait aussi
celui-ci. Les nouveaux noeuds recoivent de nouveaux ids (apres le plus grand
id de la carte), `lat/lon`, `local_x/local_y` et `ele` interpoles, puis sont
lisses comme les autres. Ils sont inseres apres le dernier `node` du fichier
et dans la liste `nd` du way.

### A2. Lissage lateral et correction de courbure

La centreligne de chaque lanelet est echantillonnee au pas `courbure.pas_reech`, par
fractions egales de la longueur des deux bornes (meme construction que
`generateFineCenterline` dans Autoware). Chaque echantillon est une
combinaison lineaire fixe des noeuds des bornes (`C = M·P`), et `L = M_l·P`
est le vecteur gauche - droite. Les variables sont les **noeuds** `X`, pas les
echantillons. Deux passes, sous un **meme budget** `courbure.ecart_max` mesure depuis
la geometrie d'origine :

**1. Lissage lateral global** — applique a tous les lanelets selectionnes,
qu'ils depassent ou non `k_max`, pour une courbure a **variation continue** :

```
min  λ Σ ‖Δ³C‖²  +  λ μ₁ Σ ‖L - L⁰‖²  +  Σ aᵢ ‖Xᵢ - X⁰ᵢ‖²,     λ = (fenetre / 2π·pas)⁶
```

`Δ³C` (difference troisieme) est la derivee de la courbure : on penalise les
sauts de courbure, pas la courbure elle-meme. `mode_aoi.fenetre_laterale` (10 m) est
la longueur d'onde de coupure du filtre. Le budget d'ecart est tenu **noeud
par noeud** : l'ancrage `aᵢ` des seuls noeuds qui le depassent est renforce
iterativement, pour ne pas affaiblir le lissage sur toute l'AOI a cause d'un
seul noeud contraint. Le rapport donne `|dk/ds|` (max et p95) avant/apres.

Sur trois zones de test, 10 m est le meilleur compromis : p95 de `|dk/ds|`
reduit de 25 a 30 % pour un deplacement moyen de ~6 cm ; a 20 m, la moitie
des noeuds butent sur le budget pour un gain inferieur.

**2. Correction locale des depassements** — la ou `k > k_max` subsiste,
meme formulation que pour la trajectoire (difference seconde, poids en
cloche autour des violations, `λ` par bissection) :

```
min  λ Σ wᵢ ‖Δ²C‖²  +  λ μ₂ Σ wᵢ ‖L - L⁰‖²  +  Σ ‖X - X_dep‖²
```

Le terme de largeur (`μ₁ = POIDS_LARGEUR_LISSAGE = 1.0`,
`μ₂ = POIDS_LARGEUR = 0.3`) evite qu'une borne commune a deux voies opposees
se deforme pour servir les deux centrelignes : la largeur derive de 0.18 m
sans lui, de 0.04 m avec.

Le lissage lateral ne doit pas degrader la courbure : s'il ajoute un
depassement de `k_max` ou releve la courbure maximale au-dela de
`max(k_max, k d'origine)`, `λ` est reduit par bissection (`attenue pour ne pas
degrader la courbure max` dans le rapport).

`mode_aoi.sans_lissage_lateral` desactive la passe 1, `courbure.sans_courbure` la passe 2.

### Raccords avec les voies hors AOI

Les extremites d'un lanelet traite sont partagees avec ses predecesseurs /
successeurs hors AOI : elles restent en place (continuite de position), et
ces voisins figes entrent dans le calcul comme conditions aux limites. Les
differences troisiemes **a cheval sur un raccord** sont ponderees par
`POIDS_RACCORD = 10` : la centreligne corrigee doit prolonger la tangente et
la courbure de la voie voisine. Le rapport donne, aux raccords, l'angle de
cassure de la centreligne et `|dk/ds|` avant/apres :

```
  raccords aux voisins figes  : 5, cassure max 23.50 -> 9.18 deg, |dk/ds| max 0.1502 -> 0.1274 1/m2
```

Cassure max aux raccords sur trois zones de test :

| Zone | Origine | Sans insertion | `noeuds_ajoutes_max: 5` |
|---|---|---|---|
| F | 23.5 deg | 21.1 deg | 9.2 deg |
| G | 26.2 deg | 19.0 deg | 13.4 deg |
| H | 17.8 deg | 13.2 deg | 6.3 deg |

La cassure residuelle vient des voisins figes eux-memes : si la voie hors AOI
arrive deja avec un angle au noeud de jonction, seule une AOI plus grande
permet de le reprendre.

**Mesure de la courbure.** Les bornes sont des polylignes (pas median entre
noeuds ~3.9 m, angles aux sommets jusqu'a 30 deg) : la courbure de Menger a
1 m de base verrait un virage serre a chaque sommet. Elle est donc mesuree sur
des triplets espaces de `± mode_aoi.base_courbure`, par defaut **l'empattement** du
vehicule. Sur `lanelet2_map.osm` le nombre de lanelets `road` en exces passe
de 94 (base 1 m) a 46 (base 3 m) ; la base a donc un effet direct sur le
resultat.

### A3. Lissage de l'altitude

Pas de recalage : l'altitude existante est lissee **le long des bornes**
(jamais en travers, le devers est conserve), en resolvant

```
min  Σ mᵢ (zᵢ - z⁰ᵢ)²  +  ℓ² Σ₍ᵢⱼ₎ (zᵢ - zⱼ)² / dᵢⱼ ,   ℓ = fenetre / √24
```

discretisation d'un lissage a noyau exponentiel de meme variance qu'une
moyenne glissante de largeur `altitude.fenetre_z`. Les noeuds figes gardent leur
altitude : pas de marche a la frontiere de la zone. `mode_trajectoire.offset_capteur`,
`mode_trajectoire.distance_max` et `mode_trajectoire.ecart_min_z` n'interviennent pas.

### lat / lon

La projection de la carte n'est pas connue du script. Le deplacement en
`lat/lon` est obtenu par un jacobien local `d(lat,lon)/d(x,y)` ajuste sur
les noeuds des bornes traitees (residu typique ~1 mm), les noeuds dont
`lat/lon` ne concordent pas avec `local_x/local_y` etant ecartes.

### Limite connue : virages de carrefour

Sur les lanelets de tourne-a-gauche / tourne-a-droite d'un carrefour, la
courbure est souvent **non corrigeable quel que soit `courbure.ecart_max`** : leurs
extremites sont les points de divergence partages avec les lanelets tout
droit, qui fixent a la fois la position et la tangente de depart. Augmenter le
rayon demanderait de reculer ces points de jonction. Ces zones sont listees
dans le rapport (colonne `lanelets`) et donnent le code de retour 2 : c'est
le trace du carrefour qu'il faut revoir.

---

## Parametres

Tous dans [`config/config.yaml`](config/config.yaml), valeurs livrees
ci-dessous. Le fichier separe les parametres **communs** aux deux modes, ceux
du **mode AOI** (`mode_aoi`) et ceux du **mode trajectoire**
(`mode_trajectoire`, sans effet en mode AOI).

### `aoi`

| Cle | Valeur | Role |
|---|---|---|
| `aoi` | `true` | `true` : mode AOI, correction directe des lanelets entierement couverts par une AOI. `false` : mode trajectoire. |

### Communs : `entrees_sorties`

| Cle | Valeur | Role |
|---|---|---|
| `carte` | `lanelet2_map.osm` | Carte Lanelet2. |
| `sortie_carte` | `null` : `<carte>_OUT.osm` | Carte produite. |
| `rapport` | `null` : `log/<carte>_rapport_<AAAAMMJJ_HHMMSS>.txt` | Rapport texte, toujours horodate (un chemin donne recoit aussi l'horodatage). |
| `graphique` | `false` | Figures PNG. |
| `silencieux` | `false` | Pas de sortie console (le rapport est ecrit). |

### Communs : `vehicule`

| Cle | Valeur | Role |
|---|---|---|
| `fichier` | `config/vehicle_model.yaml` | Fichier de modeles, ou `vehicle_info.param.yaml` Autoware. |
| `modele` | `null` : `modele_par_defaut` | Modele a utiliser. |
| `lister_vehicules` | `false` | Affiche les modeles et quitte. |
| `empattement` | `null` : fichier vehicule | Empattement (m). |
| `angle_braquage_max` | `null` : fichier vehicule | Braquage max des roues avant (deg). |
| `marge_courbure` | `null` : fichier vehicule | Facteur de garde (`<1` = plus severe). |

### Communs : `courbure`

| Cle | Valeur | Role |
|---|---|---|
| `ecart_max` | 0.30 m | Deplacement lateral maximal d'un point (trajectoire) ou d'un noeud (AOI). |
| `pas_reech` | 1.0 m | Pas de reechantillonnage de la trajectoire / des centrelignes. |
| `sans_courbure` | `false` | Desactive l'etape 4 (passe 2 de A2 en mode AOI). |

### Communs : `altitude`

| Cle | Valeur | Role |
|---|---|---|
| `fenetre_z` | 25.0 m | Fenetre de lissage. |
| `sans_lissage_z` | `false` | Desactive le lissage. |
| `precision` | 4 | Decimales ecrites dans `ele`, `local_x`, `local_y`. |

### Mode AOI : `mode_aoi`

| Cle | Valeur | Role |
|---|---|---|
| `base_courbure` | `null` : empattement | Demi-base de mesure de la courbure des centrelignes (m). |
| `fenetre_laterale` | 10.0 m | Longueur d'onde de coupure du lissage lateral. |
| `sans_lissage_lateral` | `false` | Desactive le lissage lateral (passe 1 de A2). |
| `noeuds_ajoutes_max` | 1 | Nombre maximal de noeuds inseres par borne (0 = aucun). |

### Mode trajectoire : `mode_trajectoire`

Lus et verifies dans les deux modes, mais **sans effet en mode AOI**.

| Cle | Valeur | Role |
|---|---|---|
| `traj` | `traj.csv` | Trajectoire CSV. |
| `sortie_traj` | `null` : `<traj>_OUT.csv` | Trajectoire corrigee. |
| `pas_min` | 0.30 m | Distance min entre points bruts conserves. |
| `seuil_rebroussement` | 90 deg | Variation de cap au-dela de laquelle on coupe. |
| `distance_max` | 12.0 m | Rayon max d'appariement. |
| `ecart_min_z` | 0.20 m | En deca, le noeud n'est pas modifie. |
| `offset_capteur` | `null` : `hauteur_capteur` du modele | Hauteur soustraite a Z (m). |
| `ele_nul_mediane` | `false` | Altitude mediane pour les noeuds `ele=0` sans voisin. |

Les constantes internes des algorithmes (`LAMBDA_*`, `POIDS_LARGEUR*`,
`POIDS_RACCORD`, `N_ANCRAGE`...) restent dans `utils/courbure.py` et
`utils/aoi.py`.

---

## Fichiers produits

| Fichier | Contenu |
|---|---|
| `<carte>_OUT.osm` | Carte Lanelet2. **Seul le tag `ele` est modifie** ; geometrie, ways, relations et autres tags sont conserves a l'identique. |
| `<traj>_OUT.csv` | Trajectoire nettoyee, reechantillonnee, corrigee, `yaw` recalcule. |
| `log/<carte>_rapport_<AAAAMMJJ_HHMMSS>.txt` | Journal complet (identique a la sortie console), horodate : un lancement n'ecrase jamais le rapport precedent. |
| `<carte>_OUT.png` | Avec `entrees_sorties.graphique` : altitude brute, altitude lissee, courbure vs limite. En mode AOI : vue en plan avant/apres, courbure des centrelignes, altitude des noeuds mobiles. |

En mode AOI, pas de `<traj>_OUT.csv`, et `<carte>_OUT.osm` porte aussi les
nouvelles positions (`local_x`, `local_y`, `lat`, `lon`) des noeuds deplaces.

---

## Lire le rapport

Le bloc a surveiller est celui de l'etape 4 :

```
  points en exces de courbure : 42 -> 32

  /!\ 3 zone(s) NON corrigeables dans le budget d'ecart :
      x              y              R reel(m)  k(1/m)     braquage requis
      56271.40       15584.99       5.42       0.1846     38.7 deg
```

Lecture : a cet endroit le trace impose un rayon de 5.42 m, soit **38.7 deg de
braquage**, alors que le vehicule choisi n'en offre que 34.38. Trois issues :

1. **Augmenter `courbure.ecart_max`** si la voie le permet. Attention : si le
   deplacement lateral effectivement utilise reste bien en dessous du budget
   (ligne `deplacement lateral`), la contrainte n'est pas le budget mais la
   geometrie elle-meme — augmenter le budget ne changera rien.
2. **Revoir le trace de la carte** a ces coordonnees.
3. **Affecter un vehicule plus court** a cette ligne.

---

## Points de vigilance

### `hauteur_capteur` n'est pas verifiee

La valeur `0.26` m est reprise telle quelle de l'ancien script et appliquee a
tous les modeles. **Elle n'a ete confirmee pour aucun vehicule.** Autoware ne
la fournit pas dans `vehicle_info.param.yaml`.

Il s'agit de la hauteur, au-dessus du sol, du repere dans lequel la
trajectoire a ete enregistree. Si c'est `base_link` et que `base_link` est au
niveau du sol, la valeur correcte est `0.0` — sinon toute la carte est
decalee de 26 cm en altitude. **A mesurer avant exploitation.**

### Couverture partielle

Seuls les noeuds situes a moins de `mode_trajectoire.distance_max` d'un point de trajectoire
sont recales (50 % sur le jeu de reference : la trajectoire ne couvre que les
voies parcourues). Les autres gardent leur altitude d'origine. La ligne
`noeuds avec voisin` du rapport donne le taux de couverture.

### `mode_trajectoire.ele_nul_mediane`

Affecte l'altitude mediane de la trajectoire aux noeuds restes a `ele = 0`
sans voisin. C'est un depannage grossier : sur une carte a relief marque la
mediane globale n'a pas de sens local. Desactive par defaut.

---

## Compatibilite VectorMap

Les tags sont lus **par cle**, pas par position, donc les variantes de format
et les tags optionnels (`color`, `mgrs_code`...) sont sans effet.

Format depuis oct. 2024 :

```xml
<node id="1" lat="46.12982337722" lon="-1.0794278177">
  <tag k="local_x" v="48365.1015"/>
  <tag k="local_y" v="10264.8907"/>
  <tag k="ele" v="53.9878"/>
</node>
```

Format aout 2024 (`mgrs_code` en tete) :

```xml
<node id="48" lat="48.72244562283628" lon="2.257069031480778">
  <tag k="mgrs_code" v="31UDP453968"/>
  <tag k="local_x" v="45359.0525"/>
  <tag k="local_y" v="96868.3611"/>
  <tag k="ele" v="129.5001"/>
</node>
```

Les deux sont acceptes. Les noeuds depourvus de `local_x` / `local_y` / `ele`
sont ignores et comptes dans le rapport.

---

## Corrections par rapport a l'ancienne version

`../optimisseur_KNN/optimiseur_hauteur_KNN.py` presentait les defauts
suivants, tous corriges :

| Defaut | Consequence | Correction |
|---|---|---|
| Decoupage errone apres `np.convolve(..., 'full')` | Altitude lissee **decalee d'un echantillon** sur tout le signal | Remplissage par les bords + mode `valid` ; ecart interieur verifie a 3e-14 |
| Recherche du plus proche voisin en boucle Python `O(N*M)` | ~48 s pour 2.6e8 distances | KD-tree : chaine complete en ~1.6 s |
| Tags lus par position (`node[0]`, `node[1]`, `node[2]`) | Lecture silencieusement fausse sur le format aout 2024 | Lecture par cle |
| Chemins absolus codes en dur (`/home/admin01/...`) | Script inutilisable en l'etat | `config/config.yaml`, chemins relatifs au script |
| `str(float)` a l'ecriture | `ele` a 17 decimales | `altitude.precision` (4 par defaut) |
| `plt.show()` bloquant | Inutilisable sans affichage | Backend `Agg`, PNG optionnel |
| Titres/labels poses avant `plt.subplot` | Appliques a un axe fantome, invisibles | Figure construite par `subplots` |
| `z_median` calcule seulement si le filtrage etait actif | Noeuds `ele=0` ecrits a 0 | Mediane toujours calculee, option explicite |
| Courbure jamais verifiee | Rien ne garantissait un trace roulable | Etapes 2 a 4 |
| Constantes en dur (`-0.26`, `max_distance`, `min_dega`) | Non reglables | Modele vehicule + options |
