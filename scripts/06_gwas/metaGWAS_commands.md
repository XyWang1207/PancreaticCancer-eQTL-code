# GWAS meta-analysis commands (META v1.7)

This file records the command-line workflow used to perform GWAS meta-analysis using **META v1.7**.

> Notes  
> - Replace all input/output filenames with your local files.  
> - `./meta` refers to the META v1.7 executable in the current working directory.  
> - Cohort input files should be formatted according to META v1.7 requirements (summary statistics per cohort).

---

## Check available options

```bash
./meta
````

---

## Basic fixed-effects meta-analysis (example)

```bash
./meta \
  --method 1 \
  --threshold 0.3 \
  --cohort EUR_input.txt CHN_input.txt \
  --output EUR_CHN_meta_output.txt
```

**Arguments**

* `--method 1`: meta-analysis method implemented by META v1.7 (as used in this study)
* `--threshold 0.3`: cohort-level filter threshold 
* `--cohort`: input summary statistics files for each cohort
* `--output`: output meta-analysis result file

---

## Meta-analysis with cohort-specific genomic control lambdas (example)

If genomic control (GC) lambdas are pre-estimated for each cohort, they can be provided explicitly:

```bash
./meta \
  --method 1 \
  --lambda 1.05 1.08 \
  --cohort EUR_input.txt CHN_input.txt \
  --output EUR_CHN_meta_GC_output.txt
```

* `--lambda 1.05 1.08` corresponds to GC lambdas for the cohorts listed in `--cohort` (same order).

---

## Meta-analysis of three cohorts (example)

```bash
./meta \
  --method 1 \
  --cohort input_EUR.txt input_CHN.txt input_rep.txt \
  --output EUR_CHN_rep_meta_output.txt
```


