# -*- coding: utf-8 -*-
"""Figures de controle (matplotlib, backend Agg, jamais bloquant)."""

import numpy as np


def tracer(chemin, brut, corrige, k, k_max, jrn):
    try:
        import matplotlib
        matplotlib.use("Agg")           # jamais d'affichage bloquant
        import matplotlib.pyplot as plt
    except ImportError:
        jrn("  matplotlib absent, figure ignoree.")
        return
    fig, ax = plt.subplots(3, 1, figsize=(12, 11))
    ax[0].plot(brut, color="tab:red", lw=0.7)
    ax[0].set(title="Altitude brute", xlabel="index", ylabel="Z (m)")
    ax[1].plot(corrige, color="tab:green", lw=0.7)
    ax[1].set(title="Altitude lissee", xlabel="index", ylabel="Z (m)")
    ax[2].plot(k, color="tab:blue", lw=0.5)
    ax[2].axhline(k_max, color="tab:red", ls="--",
                  label="limite vehicule %.4f 1/m" % k_max)
    ax[2].set(title="Courbure apres correction", xlabel="index",
              ylabel="k (1/m)", yscale="log")
    ax[2].legend()
    fig.tight_layout()
    fig.savefig(chemin, dpi=110)
    plt.close(fig)
    jrn("  figure ecrite : %s" % chemin)


def tracer_aoi(chemin, zones, bornes, selection, xy_de, cl, P, k_av, k_ap,
               k_max, z_av, z_ap, residus, jrn):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        jrn("  matplotlib absent, figure ignoree.")
        return
    xy_ap = {int(i): p for i, p in zip(cl.ids, P)}
    fig, ax = plt.subplots(3, 1, figsize=(12, 15),
                           gridspec_kw={"height_ratios": [2, 1, 1]})
    for wid, poly in zones:
        f = np.vstack([poly, poly[:1]])
        ax[0].fill(f[:, 0], f[:, 1], color="tab:orange", alpha=0.15)
        ax[0].plot(f[:, 0], f[:, 1], color="tab:orange", lw=1.2)
    for lid in selection:
        for b in bornes[lid]:
            av = np.array([xy_de[int(i)] for i in b])
            ap = np.array([xy_ap.get(int(i), xy_de[int(i)]) for i in b])
            ax[0].plot(av[:, 0], av[:, 1], color="0.6", lw=0.8, ls="--")
            ax[0].plot(ap[:, 0], ap[:, 1], color="tab:blue", lw=1.0, marker=".",
                       ms=3)
    C = cl.M @ P
    ax[0].plot(C[:, 0], C[:, 1], ".", color="tab:green", ms=1.5)
    for r in residus:
        ax[0].plot(r["x"], r["y"], "x", color="tab:red", ms=10, mew=2)
    ax[0].set(title="Lanelets de l'AOI : avant (gris) / apres (bleu), "
                    "centrelignes (vert), residus (x rouge)",
              xlabel="local_x (m)", ylabel="local_y (m)", aspect="equal")
    # cadrage sur les AOI : les voisins figes peuvent etre tres longs
    coins = np.vstack([poly for _, poly in zones])
    ax[0].set_xlim(coins[:, 0].min() - 10.0, coins[:, 0].max() + 10.0)
    ax[0].set_ylim(coins[:, 1].min() - 10.0, coins[:, 1].max() + 10.0)
    if k_av is not None:
        ax[1].plot(k_av, ".", color="tab:red", ms=2, label="avant")
        ax[1].plot(k_ap, ".", color="tab:green", ms=2, label="apres")
        ax[1].axhline(k_max, color="tab:red", ls="--",
                      label="limite vehicule %.4f 1/m" % k_max)
        ax[1].legend()
    ax[1].set(title="Courbure des centrelignes", xlabel="triplet",
              ylabel="k (1/m)")
    ax[2].plot(z_av, ".", color="tab:red", ms=3, label="avant")
    ax[2].plot(z_ap, ".", color="tab:green", ms=3, label="apres")
    ax[2].set(title="Altitude des noeuds mobiles", xlabel="noeud",
              ylabel="ele (m)")
    ax[2].legend()
    fig.tight_layout()
    fig.savefig(chemin, dpi=110)
    plt.close(fig)
    jrn("  figure ecrite : %s" % chemin)
