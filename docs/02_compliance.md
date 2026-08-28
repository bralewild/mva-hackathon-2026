# Cumplimiento de las Reglas Oficiales

Checklist operativo derivado de las *Official Rules* del MVA Hackathon 2026.
Este documento existe para que el cumplimiento sea **verificable**, no una
promesa: cada obligación tiene un mecanismo concreto asociado.

---

## 1. Privacidad y restricciones de recontacto

| Obligación | Cómo se cumple |
|---|---|
| No recontactar al paciente, su familia ni a la MVA Society | Ningún paso del pipeline consulta identidad, contacto ni redes. Solo APIs científicas públicas (Ensembl, NCBI, HPO). |
| No liberar ni dar acceso a los datos a nadie | Los datos viven únicamente en `~/mva/` dentro de WSL. Nunca en el repo, nunca en un servicio externo. |
| Salvaguardas contra uso no autorizado | Tres barreras, ver §3. |
| Reportar cualquier divulgación no autorizada | Sage Help Center. Contacto documentado abajo. |

**No se subió ningún dato del paciente a ningún servicio.** Las consultas a
Ensembl VEP envían coordenadas genómicas (`15 40209701 T G`) sin identificador
alguno del sujeto — es exactamente el uso previsto de una API de anotación.

---

## 2. Borrado de datos — 30 días desde el cierre

Las Reglas Oficiales son más amplias que el formulario del dataset:

> *"All data must be deleted within 30 days of Hackathon close from all
> environments (local machines, cloud instances, notebooks, **private repos**,
> and any **intermediate or derived datasets**)."*

### Qué hay que borrar

Correr `pipeline/99_data_inventory.sh` para el inventario exacto. Resumen:

| Categoría | Ruta | ¿Se borra? |
|---|---|---|
| VCF, índice, documento clínico | `~/mva/data/raw/` | **SÍ** |
| Términos HPO extraídos | `~/mva/data/raw/patient_hpo.tsv` | **SÍ** (derivado) |
| VCF anotado, candidatas, cache VEP | `~/mva/work/` | **SÍ** (derivados) |
| Reportes con coordenadas del paciente | `~/mva/results/` | **SÍ** (derivados) |
| Reportes espejados en el proyecto | `results/` | **SÍ** (derivados) |
| Referencia genómica, BD snpEff, HPO | `~/micromamba`, `~/mva/data/annot` | No — recursos públicos |
| Código del pipeline | repo git | No — sin datos del paciente |

### Comando

```bash
rm -rf ~/mva/data/raw ~/mva/work ~/mva/results
rm -rf "$PROJECT/results"
```

### Notificación

Hay **dos direcciones** en las fuentes oficiales, así que se notifica a ambas:

| Fuente | Dirección |
|---|---|
| Official Rules | `RarediseaserealkidMVAhackathon2026@synapse.org` |
| Formulario del dataset | `MVAHackathon2026@synapse.org` |

> *"If you do not contact us, we may contact you directly in 30 days."*

---

## 3. Salvaguardas implementadas

Tres barreras independientes para que el genoma del chico no llegue a GitHub:

1. **Separación física.** Los datos viven en ext4 dentro de WSL; el repo está en
   NTFS. No es que estén ignorados: **no están en la carpeta del repo.**
2. **`.gitignore`.** Bloquea `*.vcf*`, `*.bam`, `*.cram`, `*.fastq*`, `*.docx`,
   `patient_hpo.tsv`, `results/*.csv`, `results/*.tsv`, `results/*.html`,
   `data/`, `work/`.
3. **Auditoría.** `99_data_inventory.sh` lista dónde está cada byte y clasifica
   qué es dato del paciente y qué es recurso público.

Verificación rápida antes de cada push:

```bash
git ls-files | grep -E '\.(vcf|bam|cram|fastq|docx)$|patient_hpo'
# no debe devolver nada
```

---

## 4. Reproducibilidad — es un requisito funcional

> *"By registering, participants acknowledge that submissions **may be rerun**
> by the Hackathon organizers."*

Los organizadores pueden **ejecutar este código**. No alcanza con que se vea
prolijo:

- `pipeline/run_all.sh` corre el pipeline de punta a punta
- Cada etapa es **idempotente**: si su salida existe y es válida, se saltea
- La etapa 04 tiene **cache incremental y reanudación** ante cortes de red
- El entorno está declarado en el README (WSL2 + micromamba, versiones exactas)
- `.gitattributes` fuerza `LF` para que los scripts corran en Linux

**Requisito para quien reejecute:** debe tener su propio acceso autorizado al
dataset. El pipeline no incluye datos del paciente por diseño; `00b` extrae los
términos HPO del `.docx` original, que cada participante debe obtener por la vía
gated de Hugging Face.

---

## 5. Embargo y publicación

| Acción | ¿Permitida? |
|---|---|
| Publicar **código, modelos y salidas derivadas** | **Sí, en cualquier momento** |
| Manuscritos para revisión por pares durante el embargo | No |
| Abstracts o pósters en congresos | Sí, con **aprobación escrita previa** |

El embargo empieza al cierre y termina cuando los organizadores publiquen su
reporte resumen o preprint.

**Consecuencia práctica:** mantener el repo privado hasta el cierre es una
decisión **estratégica** (que otros equipos no copien el método mientras la
competencia está abierta), no una obligación legal.

---

## 6. Acknowledgement obligatorio

Debe incluirse **textual** en toda publicación, preprint, abstract o
comunicación pública derivada del hackathon:

> "This work was made possible through the Hackathon, organized by Sage
> Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON
> (The Benchmarking, Evaluation, and Assessment Consortium for Science), with
> prize sponsorship from AWS and Anthropic. We are deeply grateful to the child
> and their family who generously contributed their data and their story to
> advance research into this rare disease. We acknowledge their trust in making
> this Hackathon possible."

Además:

- **Privacidad del sujeto:** ninguna publicación puede incluir información que
  permita reidentificar al paciente o su familia más allá de lo ya público en
  los blogs y comunicaciones de la propia familia.
- **Cita del dataset:** usar la referencia indicada en la página de Synapse del
  hackathon al momento de publicar.
- **Aprobación ética:** protocolo aprobado por **WCG IRB #20252010**.
- **Licencia:** las submissions se liberan bajo **CC BY 4.0**, con atribución
  nominal del participante.

---

## 7. Entregables por track

| | Track 1 | Track 2 |
|---|---|---|
| Envíos permitidos | 6 | **1** |
| CSV de predicciones | Sí | — |
| Reporte escrito | Sí | Sí |
| Repo público de GitHub | Sí | Sí |
| Video de pitch de 3 min | Ver nota | Sí |
| Evaluación | Automática (rank points + F-max) contra la variante **validada por el NHS** | Panel experto: rigor 35%, impacto 25%, innovación 25%, escalabilidad 15% |

> **Nota sobre el video:** las Reglas Oficiales dicen *"Each team's submission
> includes a written report, a GitHub repository, and a 3-minute recorded pitch
> video"* de forma general, mientras que la pestaña de envío de Track 1 solo
> pide CSV + reporte + repo. Ante la ambigüedad, conviene tener el video listo.

---

## 8. Otros requisitos

- Ser mayor de 18 años.
- Cada integrante de un equipo debe registrarse **individualmente** y aceptar
  las reglas. *(Participación individual como `bralewild` — no aplica.)*
- El cómputo corre por cuenta del participante. *(Todo local: WSL2 sobre
  i9-14900HX, sin costos de nube.)*
- Los datos no pueden recompartirse por ningún canal.
