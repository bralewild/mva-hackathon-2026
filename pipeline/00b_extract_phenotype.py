#!/usr/bin/env python3
"""
==============================================================================
00b_extract_phenotype.py - Extract HPO terms from the clinical document

INPUT  : $RAW/Challenge_Clinical_Phenotype_1.docx
OUTPUT : $RAW/patient_hpo.tsv        (hpo_id, label, clinical_feature)

PURPOSE / PRIVACY
-----------------
The .docx is marked "Confidential - Do not redistribute". The HPO terms ARE the
patient's clinical information, so they cannot live in a public repository. This
script extracts them and writes them alongside the rest of the patient data in
WSL, where .gitignore blocks them.

Consequence: anyone cloning the repository must run this step with their own
authorised copy of the dataset. That is the correct behaviour - and it keeps the
pipeline reproducible without leaking anything.
==============================================================================
"""
import html
import os
import re
import sys
import zipfile

from _paths import RAW  # honours MVA_BASE; see pipeline/_paths.py
DOCX = RAW + "/Challenge_Clinical_Phenotype_1.docx"
OUT = RAW + "/patient_hpo.tsv"


def docx_text(path):
    """Flatten a .docx into plain text, preserving paragraph and cell breaks."""
    xml = zipfile.ZipFile(path).read("word/document.xml").decode("utf-8", "ignore")
    xml = xml.replace("</w:p>", "\n").replace("</w:tr>", "\n").replace("</w:tc>", "\t")
    return html.unescape(re.sub(r"<[^>]+>", "", xml))


def main():
    if not os.path.exists(DOCX):
        sys.exit("missing " + DOCX)
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
        # rows look like: Feature <tab> HPO term <tab> HP:xxxxxxx <tab> Notes
        cells = [c.strip(" |") for c in line.split("\t") if c.strip(" |")]
        feature = cells[0] if cells else ""
        label = cells[1] if len(cells) > 1 else ""
        rows.append((hpo, label, feature))

    if not rows:
        sys.exit("no HP:xxxxxxx term found in the document")

    with open(OUT, "w") as f:
        f.write("hpo_id\tlabel\tclinical_feature\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    print("{} HPO terms extracted -> {}".format(len(rows), OUT))
    for hpo, label, _feat in rows:
        print("  {}  {}".format(hpo, label))


if __name__ == "__main__":
    main()
