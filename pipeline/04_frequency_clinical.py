#!/usr/bin/env python3
"""
==============================================================================
04_frequency_clinical.py - Frecuencia poblacional y evidencia clinica

ENTRADA : $WORK/03_candidates.tsv
SALIDA  : $WORK/04_annotated_candidates.tsv   (todas, con anotacion)
          $WORK/04_rare_candidates.tsv        (las que pasan el filtro)
          $WORK/04_vep_cache.jsonl            (cache incremental, permite reanudar)
          $RESULTS/04_frequency_summary.txt

PROPOSITO
---------
Una enfermedad ultra-rara (<50 personas en el mundo) no puede estar causada por
una variante frecuente en la poblacion general. Este paso aplica ese principio.

Umbral: gnomAD AF < 0.01 para modelos recesivos.
La ausencia total de gnomAD se registra como evidencia PM2 de ACMG.

FUENTE: Ensembl VEP REST (POST /vep/human/region), lotes de 200 variantes.
No se guarda ninguna base de datos local: la anotacion queda siempre al dia y
evitamos descargar decenas de GB de cache para unos pocos miles de variantes.
Devuelve ademas transcrito MANE, HGVSc/HGVSp, CADD, rsID y ClinVar.

REANUDACION (leccion aprendida)
-------------------------------
La primera version escribia resultados solo al final: 15 minutos sin ninguna
senal, y si se cortaba la red habia que empezar de cero. Un pipeline que hay
que rearrancar entero por un timeout no es reproducible en la practica.

Ahora cada lote se persiste en $WORK/04_vep_cache.jsonl apenas llega. Si el
script se corta, al volver a correrlo saltea todo lo ya cacheado. El progreso
se imprime lote a lote con flush inmediato.

LIMITES DE LA API: ~15 req/s y 55.000 req/h. Usamos ~3 req/s con reintento
exponencial ante HTTP 429/503.
==============================================================================
"""
import collections
import csv
import json
import os
import sys
import time
import urllib.error
import urllib.request

BASE = os.path.expanduser("~/mva")
IN_TSV = BASE + "/work/03_candidates.tsv"
OUT_ALL = BASE + "/work/04_annotated_candidates.tsv"
OUT_RARE = BASE + "/work/04_rare_candidates.tsv"
CACHE = BASE + "/work/04_vep_cache.jsonl"
OUT_SUM = BASE + "/results/04_frequency_summary.txt"

VEP_URL = ("https://rest.ensembl.org/vep/human/region"
           "?hgvs=1&canonical=1&mane=1&numbers=1&protein=1"
           "&check_existing=1&variant_class=1&CADD=1")
BATCH = 200
SLEEP = 0.35
MAX_RETRY = 5
AF_MAX_RECESSIVE = 0.01


def log(msg):
    """Progreso con flush inmediato: sin esto no se ve nada hasta el final."""
    print(msg, file=sys.stderr, flush=True)


def vep_batch(variants):
    """POST de un lote a VEP REST con reintento exponencial."""
    body = json.dumps({"variants": variants}).encode()
    for attempt in range(MAX_RETRY):
        req = urllib.request.Request(VEP_URL, data=body, headers={
            "Content-Type": "application/json",
            "Accept": "application/json",
        })
        try:
            with urllib.request.urlopen(req, timeout=180) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (429, 503):
                wait = 2 ** attempt
                log("    HTTP {} - reintento en {}s".format(e.code, wait))
                time.sleep(wait)
                continue
            log("    HTTP {}: {}".format(e.code, e.read()[:200]))
            return []
        except Exception as e:
            wait = 2 ** attempt
            log("    {}: {} - reintento en {}s".format(type(e).__name__, e, wait))
            time.sleep(wait)
    log("    lote abandonado tras {} intentos".format(MAX_RETRY))
    return []


def extract(rec):
    """Frecuencia gnomAD, rsID, ClinVar y el transcrito MANE de un registro VEP."""
    out = {
        "gnomad_af": None, "rsid": ".", "clinvar": ".", "mane": ".",
        "hgvsc": ".", "hgvsp": ".", "cadd": None,
        "consequence": ".", "impact": ".", "gene": ".",
    }
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
            best = tc
            break
        if best is None and tc.get("canonical"):
            best = tc
    if best is None and rec.get("transcript_consequences"):
        best = rec["transcript_consequences"][0]
    if best:
        out.update({
            "mane": best.get("mane_select", "."),
            "hgvsc": best.get("hgvsc", "."),
            "hgvsp": best.get("hgvsp", "."),
            "cadd": best.get("cadd_phred"),
            "gene": best.get("gene_symbol", "."),
            "consequence": ",".join(best.get("consequence_terms", [])),
            "impact": best.get("impact", "."),
        })
    return out


def key_of(rec_input):
    p = rec_input.split()
    return "{}_{}_{}_{}".format(p[0], p[1], p[3], p[4]) if len(p) >= 5 else None


def load_cache():
    """Lee el cache incremental para poder reanudar tras un corte."""
    ann = {}
    if not os.path.exists(CACHE):
        return ann
    with open(CACHE, encoding="utf-8") as f:
        for line in f:
            try:
                d = json.loads(line)
                ann[d["key"]] = d["ann"]
            except Exception:
                continue
    return ann


def main():
    if not os.path.exists(IN_TSV):
        sys.exit("falta " + IN_TSV + " - corre antes 03_inheritance_models.py")

    rows = list(csv.DictReader(open(IN_TSV), delimiter="\t"))
    log("candidatas de entrada: {:,}".format(len(rows)))

    order, seen = [], set()
    for r in rows:
        k = "{} {} . {} {} . . .".format(r["chrom"], r["pos"], r["ref"], r["alt"])
        kk = key_of(k)
        if kk not in seen:
            seen.add(kk)
            order.append(k)
    log("variantes unicas: {:,}".format(len(order)))

    ann = load_cache()
    if ann:
        log("cache: {:,} variantes ya anotadas, se saltean".format(len(ann)))

    pending = [k for k in order if key_of(k) not in ann]
    nb = (len(pending) + BATCH - 1) // BATCH
    log("lotes pendientes: {}".format(nb))

    t0 = time.time()
    with open(CACHE, "a", encoding="utf-8") as cf:
        for i in range(0, len(pending), BATCH):
            chunk = pending[i:i + BATCH]
            n = i // BATCH + 1
            for rec in vep_batch(chunk):
                k = key_of(rec.get("input", ""))
                if k:
                    a = extract(rec)
                    ann[k] = a
                    cf.write(json.dumps({"key": k, "ann": a}) + "\n")
            cf.flush()
            el = time.time() - t0
            eta = (el / n) * (nb - n) if n else 0
            log("  lote {}/{}  ({:,} anotadas)  ETA {:.0f} min".format(
                n, nb, len(ann), eta / 60))
            time.sleep(SLEEP)

    cols = list(rows[0].keys()) + [
        "gnomad_af", "rsid", "clinvar", "mane",
        "vep_hgvsc", "vep_hgvsp", "cadd", "vep_consequence", "vep_impact"]
    stats = collections.Counter()

    with open(OUT_ALL, "w", newline="") as fa, open(OUT_RARE, "w", newline="") as fr:
        wa = csv.DictWriter(fa, fieldnames=cols, delimiter="\t")
        wa.writeheader()
        wr = csv.DictWriter(fr, fieldnames=cols, delimiter="\t")
        wr.writeheader()
        for r in rows:
            a = ann.get("{}_{}_{}_{}".format(r["chrom"], r["pos"], r["ref"], r["alt"]), {})
            af = a.get("gnomad_af")
            out = dict(r)
            out.update({
                "gnomad_af": "" if af is None else af,
                "rsid": a.get("rsid", "."),
                "clinvar": a.get("clinvar", "."),
                "mane": a.get("mane", "."),
                "vep_hgvsc": a.get("hgvsc", "."),
                "vep_hgvsp": a.get("hgvsp", "."),
                "cadd": a.get("cadd") if a.get("cadd") is not None else "",
                "vep_consequence": a.get("consequence", "."),
                "vep_impact": a.get("impact", "."),
            })
            wa.writerow(out)
            stats["total"] += 1
            if not a:
                stats["sin_respuesta_vep"] += 1
            if af is None:
                stats["ausente_gnomad"] += 1
            if af is None or af < AF_MAX_RECESSIVE:
                wr.writerow(out)
                stats["raras"] += 1
            else:
                stats["comunes_descartadas"] += 1

    txt = "\n".join([
        "=" * 66,
        " 04 - FRECUENCIA POBLACIONAL Y EVIDENCIA CLINICA",
        "=" * 66,
        "",
        "  variantes evaluadas        {:>10,}".format(stats["total"]),
        "  sin respuesta de VEP       {:>10,}".format(stats["sin_respuesta_vep"]),
        "  ausentes de gnomAD (PM2)   {:>10,}".format(stats["ausente_gnomad"]),
        "  raras (AF < {})        {:>10,}".format(AF_MAX_RECESSIVE, stats["raras"]),
        "  comunes descartadas        {:>10,}".format(stats["comunes_descartadas"]),
        "",
        "  -> " + OUT_RARE,
    ])
    os.makedirs(os.path.dirname(OUT_SUM), exist_ok=True)
    open(OUT_SUM, "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
