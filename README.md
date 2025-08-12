# Ontology Work Toolkit

*A lightweight Python package distilled from the Jupyter notebook **`Ontology_Work (2).ipynb`** for loading, enriching, and analysing a phytotherapy-focused OWL/Turtle ontology.*

---

## Table of Contents

1. [Background](#background)
2. [Repository Layout](#repository-layout)
3. [Quick Start](#quick-start)
4. [Command-line Reference](#command-line-reference)
5. [Workflow Details](#workflow-details)
6. [Extending the Toolkit](#extending-the-toolkit)
7. [Contributing](#contributing)
8. [Citation](#citation)

---

## Background

Researchers in natural-product drug discovery often work with sprawling ontologies that link **plants → chemicals → targets → therapeutic roles**. The original notebook captured a proof-of-concept pipeline for:

* Parsing large OWL/Turtle files
* Adding computed fingerprints (ECFP) from a CSV
* Running SPARQL metrics queries
* Exporting summary statistics for publication figures

This repository converts that exploratory notebook into **clean, script-first code** that can be version-controlled, tested, and deployed.

---

## Repository Layout

```
.
├── src/
│   └── ontology_work.py      # main script (refactored notebook code)
├── data/                     # ontologies / CSVs here 
├── requirements.txt          # pinned dependencies for SRC.py
├── CTAPred/
│   └── predict_targets.py    # CTAPred pipeline drivers script
│   └── SharedFunc.py         # CTApred helper functions
│   └── requirements.txt      # dependencies for CTAPred pipeline
│   └── README.md             # instructions for using the CTAPred pipeline
└── README.md                 # you are here
```

| File / Dir             | Purpose                                                     |
| ---------------------- | ----------------------------------------------------------- |
| `SRC.py`               | Entry-point that orchestrates parsing, enrichment & metrics |
| `data/`                | Placeholder for ontology TTL/RDF files and CSV inputs       |
| `requirements.txt`     | Packages required to run ontology workflow                  |
| `CTAPred/`             | Modified CTAPred pipeline for predicting protein-phytochemical interactions |


---

## Quick Start

1. **Clone & install**

   ```bash
   git clone https://github.com/<your-org>/ontology-work.git
   cd ontology-work
   python -m venv .venv          # optional but recommended
   source .venv/bin/activate     # Windows: .venv\Scripts\activate
   pip install -r requirements.txt
   ```

2. **Run the script**

   ```bash
   python -m src.ontology_work \
         --ontology data/phytotherapies.ttl \
         --ecfp-csv data/phytotherapy_ecfp.csv \
         --out metrics.json
   ```

   Results:

   * `metrics.json` – counts of classes, individuals, triples, etc.
   * `phytotherapies_enriched.ttl` – copy of the ontology with `hasECFP` data properties inserted.

---

## Command-line Reference

```text
usage: ontology_work.py [-h] --ontology PATH [--ecfp-csv PATH]
                        [--out PATH] [--log-level {DEBUG,INFO,WARNING,ERROR}]

optional arguments:
  --ontology PATH      input OWL/Turtle/RDF file (required)
  --ecfp-csv PATH      CSV with columns `smiles,ecfp` for enrichment
  --out PATH           where to save metrics JSON (default: stdout)
  --log-level LEVEL    set log verbosity (default: INFO)
```

---

## Workflow Details

1. **Parse ontology** using **rdflib** – handles TTL, RDF/XML, OWL.
2. **Compute / insert fingerprints** – matches individuals by `hasSmiles` and adds `hasECFP`.
3. **SPARQL metrics** – counts classes, properties, individuals, & custom queries.
4. **Export artifacts** – enriched TTL plus JSON metrics for figures or dashboards.

---

## Extending the Toolkit

* **Modularise** – split into smaller modules and expose a proper Python API.
* **Unit tests** – add `pytest` and GitHub Actions for continuous integration.
* **Visualisation** – integrate network diagrams or heatmaps (e.g., via `matplotlib`).
* **Docker** – package the pipeline for reproducible runs on any machine.

---

## Contributing

1. Fork → feature branch → commit + push → Pull Request.
2. Follow [PEP 8](https://peps.python.org/pep-0008/) and add/ update tests.
3. Be descriptive in PR titles & commit messages.

---

## Citation

If you use this code in your research, please cite:

```bibtex
@software{hewryk_ontology_work_2025,
  author       = {Hewryk, Oresta S. I. and Pan, Ian Tong and Romano, Joseph D.},
  title        = {{Ontology Work}: A Python toolkit for phytotherapy ontology enrichment},
  year         = {2025},
  publisher    = {GitHub},
  url          = {https://github.com/RomanoLab/ontology-work}
}
```

