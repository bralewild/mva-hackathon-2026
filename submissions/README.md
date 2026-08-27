# Submissions

Los archivos `*.csv` de predicciones **no estan versionados en este repo**
mientras la competencia siga abierta (cierre: 24 oct 2026).

**Motivo:** el challenge exige un repositorio publico. Publicar las coordenadas
causales antes del cierre equivaldria a entregarle la respuesta a los demas
equipos. Lo que se juzga aqui es el **metodo** — reproducibilidad, innovacion y
escalabilidad — no el resultado masticado.

El CSV se sube directamente al formulario del challenge, que es privado hasta
la evaluacion. Se agregara a este repo despues del cierre para dejar el
registro completo.

## Convencion de nombres (exigida por el challenge)

```
<usuario-hf>_<nombre-del-modelo>.csv     ->  bralewild_blind-wgs-triage.csv
<usuario-hf>_track1_report.md            ->  bralewild_track1_report.md
```

## Formato (verificado contra tabs/submit_track1.py y evaluation.py)

| Campo | Tipo | Nota |
|---|---|---|
| `proband_id` | string | **debe ser `PROBAND01`** |
| `chrom_1` / `chrom_2` | string | **con prefijo `chr`** (ej. `chr15`) |
| `pos_1` / `pos_2` | int | GRCh38 |
| `ref_*` / `alt_*` | string | |
| `epcr` | float | rango `(0, 1]`, filas ordenadas descendente |
| `finding_type` | string | `primary` o `secondary` |
| `notes` | string | opcional, justificacion |

Maximo 10 filas. Maximo **6 submissions** por participante; solo la de mayor
puntaje aparece en el leaderboard.

`evaluation.py` **no normaliza** los nombres de cromosoma: compara tuplas
exactas `(chrom, int(pos), ref.upper(), alt.upper())`. El VCF de origen usa
naming Ensembl (`15`), asi que el submission debe anteponer `chr`.
