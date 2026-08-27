#!/usr/bin/env python3
"""
==============================================================================
04_frequency_clinical.py — Frecuencia poblacional y evidencia clinica

ENTRADA : $WORK/03_candidates.tsv
SALIDA  : $WORK/04_annotated_candidates.tsv   (todas, con anotacion)
          $WORK/04_rare_candidates.tsv        (solo las que pasan frecuencia)
          $RESULTS/04_frequency_summary.txt

PROPOSITO
---------
Una enfermedad ultra-rara (<50 personas en el mundo) no puede estar causada por
una variante frecuente en la poblacion general. Este paso aplica ese principio.

Umbrales (docs/02_metodologia.md):
  AR_COMPOUND_HET / AR_HOMOZYGOUS : gnomAD AF < 0.01
  Ausencia total de gnomAD        : evidencia PM2 de ACMG

FUENTE: Ensembl VEP REST (POST /vep/human/region), lotes de 200 variantes.
No se almacena ninguna base de datos local — todo por API, reproducible.
Devuelve ademas: transcrito MANE, HGVSc/HGVSp, CADD, rsID y ClinVar.

NOTA sobre limites: la API permite ~15 req/s y 55.000 req/h. Usamos 3 req/s con
reintento exponencial ante HTTP 429.
==============================================================================
"""
import csv, json, os, sys, time, urllib.request, urllib.error, collections

BASE    = os.path.expanduser("~/mva")
IN_TSV  = f"{BASE}/work/03_candidates.tsv"
OUT_ALL = f"{BASE}/work/04_annotated_candidates.tsv"
OUT_RARE= f"{BASE}/work/04_rare_candidates.tsv"
OUT_SUM = f"{BASE}/results/04_frequency_summary.txt"

VEP_URL = ("https://rest.ensembl.org/vep/human/region"
           "?hgvs=1&canonical=1&mane=1&numbers=1&protein=1"
           "&check_existing=1&variant_class=1&CADD=1")
BATCH, SLEEP, MAX_RETRY = 200, 0.35, 5
AF_MAX_RECESSIVE = 0.01


def vep_batch(variants):
    """POST un lote a VEP REST. Devuelve lista de records."""
    body = json.dumps({"variants": variants}).encode()
    for attempt in range(MAX_RETRY):
        req = urllib.request.Request(VEP_URL, data=body, headers={
            "Content-Type": "application/json", "Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (429, 503):
                wait = 2 ** attempt
                print(f"    HTTP {e.code}, reintento en {wait}s", file=sys.stderr)
                time.sleep(wait); continue
            print(f"    HTTP {e.code}: {e.read()[:200]}", file=sys.stderr); return []
        except Exception as e:
            wait = 2 ** attempt
            print(f"    {type(e).__name__}: {e} — reintento en {wait}s", file=sys.stderr)
            time.sleep(wait)
    return []


def extract(rec):
    """Saca frecuencia gnomAD, rsID, ClinVar y el transcrito MANE."""
    out = dict(gnomad_af=None, rsid=".", clinvar=".", mane=".",
               hgvsc=".", hgvsp=".", cadd=None, consequence=".", impact=".", gene=".")
    for cv in rec.get("colocated_variants", []) or []:
        if cv.get("id", "").startswith("rs") and out["rsid"] == ".":
            out["rsid"] = cv["id"]
        cs = cv.get("clin_sig")
        if cs:
            out["clinvar"] = ",".join(cs)
        for _alt, d in (cv.get("frequencies") or {}).items():
            for k in ("gnomade", "gnomadg", "gnomad"):
                if d.get(k) is not None:
                    v = float(d[k])
                    out["gnomad_af"] = v if out["gnomad_af"] is None else max(out["gnomad_af"], v)
    best = None
    for tc in rec.get("transcript_consequences", []) or []:
        if tc.get("mane_select"):
            best = tc; break
        if best is None and tc.get("canonical"):
            best = tc
    if best is None and rec.get("transcript_consequences"):
        best = rec["transcript_consequences"][0]
    if best:
        out.update(
            mane=best.get("mane_select", "."),
            hgvsc=best.get("hgvsc", "."), hgvsp=best.get("hgvsp", "."),
            cadd=best.get("cadd_phred"), gene=best.get("gene_symbol", "."),
            consequence=",".join(best.get("consequence_terms", [])),
            impact=best.get("impact", "."))
    return out


def main():
    if not os.path.exists(IN_TSV):
        sys.exit(f"falta {IN_TSV} — corre antes 03_inheritance_models.py")

    rows = list(csv.DictReader(open(IN_TSV), delimiter="\t"))
    print(f"candidatas de entrada: {len(rows):,}")

    keyed, order = {}, []
    for r in rows:
        k = f"{r['chrom']} {r['pos']} . {r['ref']} {r['alt']} . . ."
        if k not in keyed:
            keyed[k] = r; order.append(k)

    ann, nb = {}, (len(order) + BATCH - 1) // BATCH
    for i in range(0, len(order), BATCH):
        chunk = order[i:i + BATCH]
        print(f"  lote {i//BATCH+1}/{nb} ({len(chunk)} variantes)", file=sys.stderr)
        for rec in vep_batch(chunk):
            inp = rec.get("input", "").split()
            if len(inp) >= 5:
                ann[f"{inp[0]}_{inp[1]}_{inp[3]}_{inp[4]}"] = extract(rec)
        time.sleep(SLEEP)

    cols = list(rows[0].keys()) + ["gnomad_af", "rsid", "clinvar", "mane",
                                   "vep_hgvsc", "vep_hgvsp", "cadd",
                                   "vep_consequence", "vep_impact"]
    stats = collections.Counter()
    with open(OUT_ALL, "w", newline="") as fa, open(OUT_RARE, "w", newline="") as fr:
        wa = csv.DictWriter(fa, fieldnames=cols, delimiter="\t"); wa.writeheader()
        wr = csv.DictWriter(fr, fieldnames=cols, delimiter="\t"); wr.writeheader()
        for r in rows:
            a = ann.get(f"{r['chrom']}_{r['pos']}_{r['ref']}_{r['alt']}", {})
            af = a.get("gnomad_af")
            r = dict(r, gnomad_af=("" if af is None else af), rsid=a.get("rsid", "."),
                     clinvar=a.get("clinvar", "."), mane=a.get("mane", "."),
                     vep_hgvsc=a.get("hgvsc", "."), vep_hgvsp=a.get("hgvsp", "."),
                     cadd=(a.get("cadd") if a.get("cadd") is not None else ""),
                     vep_consequence=a.get("consequence", "."),
                     vep_impact=a.get("impact", "."))
            wa.writerow(r); stats["total"] += 1
            if af is None:
                stats["ausente_gnomad"] += 1
            if af is None or af < AF_MAX_RECESSIVE:
                wr.writerow(r); stats["raras"] += 1
            else:
                stats["comunes_descartadas"] += 1

    txt = "\n".join([
        "=" * 66, " 04 — FRECUENCIA POBLACIONAL Y EVIDENCIA CLINICA", "=" * 66, "",
        f"  variantes evaluadas        {stats['total']:>10,}",
        f"  ausentes de gnomAD (PM2)   {stats['ausente_gnomad']:>10,}",
        f"  raras (AF < {AF_MAX_RECESSIVE})        {stats['raras']:>10,}",
        f"  comunes descartadas        {stats['comunes_descartadas']:>10,}", "",
        f"  -> {OUT_RARE}"])
    os.makedirs(os.path.dirname(OUT_SUM), exist_ok=True)
    open(OUT_SUM, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
