# Rare Disease, Real Kid — MVA Hackathon 2026

Pipeline de diagnostico genomico para un caso real de **Mosaic Variegated
Aneuploidy (MVA)**. Track 1 — identificacion de variantes causales.

* Challenge: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026
* Dataset (gated): https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data
* Ventana de envio: 24 ago 2026 – 24 oct 2026 · **maximo 6 submissions**
* Evaluacion cualitativa: 24 oct – 24 nov 2026

---

## Arquitectura: por que los datos no estan en esta carpeta

Son **dos sistemas de archivos distintos**:

| | Ubicacion | Contenido | Motivo |
|---|---|---|---|
| **Codigo** | `C:\Users\user\Documents\real-kid-mva-hackathon\` (NTFS) | scripts, docs, submissions, reportes | versionable con git |
| **Datos** | `/home/user/mva/` en WSL (ext4) | VCF, intermedios, salidas pesadas | I/O 5-10x mas rapido que `/mnt/c` |

Ademas: si el VCF viviera en este repo, un `git add .` distraido publicaria el
genoma de un menor. Que el archivo **no este aca** es defensa en profundidad;
el `.gitignore` es la segunda linea, no la primera.

### Como ver los datos desde Windows

1. **`results/`** — los reportes livianos se espejan aca con `pipeline/sync_results.sh`
2. **`DATOS-EN-WSL.lnk`** — doble clic, abre el Explorador en los datos
3. **Ruta UNC** — `\wsl.localhost\Ubuntu-24.04\home\user\mva`

---

## Estructura

```
pipeline/
  00_config.sh                rutas, umbrales, constantes compartidas
  01_qc_baseline.sh           caracterizacion del VCF crudo
  02_annotate_genomewide.sh   snpEff sobre las 5.012.204 variantes (CIEGO)
  03_inheritance_models.py    calidad + impacto + compound het / homocigosis
  04_frequency_clinical.py    gnomAD + ClinVar + CADD via Ensembl VEP REST
  05_phenotype_rank.py        ranking por los terminos HPO del paciente
  99_data_inventory.sh        inventario de fuentes (cumplimiento del DUA)
  sync_results.sh             espeja reportes livianos a results/
results/                      reportes espejados desde WSL
submissions/                  CSV para el evaluador
docs/                         metodologia, hallazgos, razonamiento ACMG
```

## Como ejecutar

Los scripts corren dentro de WSL. **Siempre con shell de login** (`bash -lc`),
porque un shell no-login no lee `/etc/profile` y el PATH del entorno `bio`
no se carga:

```bash
wsl -d Ubuntu-24.04 -- bash -lc "bash /mnt/c/Users/user/Documents/real-kid-mva-hackathon/pipeline/01_qc_baseline.sh"
```

## Entorno

WSL2 Ubuntu 24.04 + micromamba env `bio`:
bcftools/samtools/htslib 1.24 · bwa · minimap2 · whatshap · bedtools · seqkit ·
snpEff 5.4c · nextflow · bbmap · python 3.12 · OpenJDK 21+

## Datos del paciente y Data Use Agreement

Correr `pipeline/99_data_inventory.sh` para el inventario completo.

Al terminar el hackathon:

```bash
rm -rf ~/mva/data/raw ~/mva/work ~/mva/results
```

y notificar a **RarediseaserealkidMVAhackathon2026@synapse.org (y copia a MVAHackathon2026@synapse.org)**.

## Formato del submission (verificado contra el codigo del evaluador)

```
proband_id,chrom_1,pos_1,ref_1,alt_1,chrom_2,pos_2,ref_2,alt_2,epcr,finding_type,notes
```

* `proband_id` **debe ser `PROBAND01`** (hardcodeado en `tabs/submit_track1.py`)
* `chrom` con prefijo **`chr`** — el VCF usa naming Ensembl (`15`), el evaluador
  espera `chr15`. **`evaluation.py` no normaliza**: compara tuplas exactas.
* `epcr` en `(0, 1]`, filas ordenadas por EPCR descendente, maximo 10

---

## Cumplimiento normativo

Ver [docs/02_compliance.md](docs/02_compliance.md) para el checklist operativo
completo derivado de las Reglas Oficiales.

Puntos clave:

* **Borrado obligatorio** dentro de los 30 dias del cierre, incluyendo repos
  privados y todo dataset intermedio o derivado. Correr
  `pipeline/99_data_inventory.sh` y notificar a **ambas** direcciones oficiales.
* **Las submissions pueden ser reejecutadas** por los organizadores: la
  reproducibilidad es un requisito funcional, no estetico.
* **El codigo se puede publicar en cualquier momento** — el embargo solo aplica
  a manuscritos para revision por pares.
* Protocolo aprobado por **WCG IRB #20252010**. Submissions bajo **CC BY 4.0**.

## Acknowledgement

Toda publicacion derivada debe incluir textualmente:

> "This work was made possible through the Hackathon, organized by Sage
> Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON
> (The Benchmarking, Evaluation, and Assessment Consortium for Science), with
> prize sponsorship from AWS and Anthropic. We are deeply grateful to the child
> and their family who generously contributed their data and their story to
> advance research into this rare disease. We acknowledge their trust in making
> this Hackathon possible."
