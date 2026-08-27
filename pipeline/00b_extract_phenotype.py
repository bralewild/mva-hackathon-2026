#!/usr/bin/env python3
"""
==============================================================================
00b_extract_phenotype.py — Extrae los terminos HPO del documento clinico

ENTRADA : $RAW/Challenge_Clinical_Phenotype_1.docx
SALIDA  : $RAW/patient_hpo.tsv        (hpo_id, label, feature)

PROPOSITO / PRIVACIDAD
----------------------
El .docx dice "Confidential - Do not redistribute". Los terminos HPO SON
informacion clinica del paciente, asi que NO pueden vivir en el repositorio
publico. Este script los extrae del documento y los deja junto al resto de los
datos del paciente en WSL, donde el .gitignore los bloquea.

Consecuencia: cualquiera que clone el repo debe correr este paso con su propia
copia autorizada del dataset. Eso es lo correcto — y ademas hace el pipeline
reproducible sin filtrar nada.
==============================================================================
"""
import os, re, sys, zipfile, html

RAW  = os.path.expanduser("~/mva/data/raw")
DOCX = f"{RAW}/Challenge_Clinical_Phenotype_1.docx"
OUT  = f"{RAW}/patient_hpo.tsv"

def docx_text(path):
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8", "ignore")
    xml = xml.replace("</w:p>", "\n").replace("</w:tr>", "\n").replace("</w:tc>", "\t")
    return html.unescape(re.sub(r"<[^>]+>", "", xml))

def main():
    if not os.path.exists(DOCX):
        sys.exit(f"falta {DOCX}")
    txt = docx_text(DOCX)

    seen, rows = set(), []
    for line in txt.split("\n"):
        m = re.search(r"(HP:\d{7})", line)
        if not m:
            continue
        hpo = m.group(1)
        if hpo in seen:
            continue
        seen.add(hpo)
        # la fila viene como: Feature <tab> Termino HPO <tab> HP:xxxxxxx <tab> Notas
        cells = [c.strip(" |") for c in line.split("\t") if c.strip(" |")]
        feature = cells[0] if cells else ""
        label   = cells[1] if len(cells) > 1 else ""
        rows.append((hpo, label, feature))

    if not rows:
        sys.exit("no se encontro ningun termino HP:xxxxxxx en el documento")

    with open(OUT, "w") as f:
        f.write("hpo_id\tlabel\tclinical_feature\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    print(f"{len(rows)} terminos HPO extraidos -> {OUT}")
    for hpo, label, feat in rows:
        print(f"  {hpo}  {label}")

if __name__ == "__main__":
    main()
