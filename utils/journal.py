# -*- coding: utf-8 -*-
"""Journal de traitement : affichage console et rapport texte."""

from pathlib import Path


class Journal:
    def __init__(self, verbeux=True):
        self.verbeux = verbeux
        self.lignes = []

    def __call__(self, texte=""):
        self.lignes.append(texte)
        if self.verbeux:
            print(texte)

    def titre(self, texte):
        self("")
        self("=" * 74)
        self("  " + texte)
        self("=" * 74)

    def ecrire(self, chemin):
        Path(chemin).write_text("\n".join(self.lignes) + "\n", encoding="utf-8")
