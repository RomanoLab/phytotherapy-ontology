#!/usr/bin/env python3
"""
phyto_ontology_pipeline.py

Refactored, publication-quality pipeline for PSB repository:
- Validates & inserts plant instances into an ontology from COCONUT DB
- Adds chemicals & links (CMAUP + NPASS + COCONUT) with data props (SMILES, MW, formula, names)
- Classifies phytochemicals by SMILES heuristics
- Enriches with DrugCentral (targets/action types) and ATC codes
- Enriches with ChEMBL mechanisms (targets, action types)
- Computes fingerprints (MACCS) and adds ECFP from CSV
- Exports helpful CSVs and ontology stats
- Builds example heatmaps from CTA predictions

Usage (examples):
  python phyto_ontology_pipeline.py plants-from-coconut --in phyth.rdf --out augmented.rdf
  python phyto_ontology_pipeline.py add-chemicals --in augmented.rdf --out augmented_with_chems.rdf --cma_dir data/CMAUP
  python phyto_ontology_pipeline.py add-maccs --in augmented_with_chems.rdf --out with_maccs.rdf
  python phyto_ontology_pipeline.py drugcentral-enrich --in with_maccs.rdf --out enriched_dc.rdf
  python phyto_ontology_pipeline.py add-atc --in enriched_dc.rdf --out enriched_dc_atc.rdf
  python phyto_ontology_pipeline.py chembl-mechanisms --in enriched_dc_atc.rdf --chembl_csv data/ChEMBL_drug_mechanisms.csv --out final_chembl.rdf
  python phyto_ontology_pipeline.py export-links --in final_chembl.rdf --csv Chemical_SMILES_isDerivedFrom_PlantLabel.csv
  python phyto_ontology_pipeline.py stats --in final_with_GeneticEntityOnly.ttl
  python phyto_ontology_pipeline.py add-ecfp --in final_with_GeneticEntityOnly.ttl --ecfp_csv phytotherapy_ecfp.csv --out final_with_ECFP.ttl
  python phyto_ontology_pipeline.py heatmap --pred cta_predictions_all_phytotherapies.csv --out fig_heatmap_top50x50.png
"""

# ========== DEPENDENCIES ==========
# pip install rdflib biopython pandas rdkit-pypi requests tqdm psycopg2-binary matplotlib
# (Optional) seaborn for nicer heatmaps: pip install seaborn
# ==================================

from __future__ import annotations
import argparse
import csv
import dataclasses
import logging
import os
import re
import sys
import time
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Set, Tuple

import pandas as pd
import requests
from tqdm import tqdm

from rdflib import Graph, Namespace, RDF, RDFS, OWL, Literal, URIRef
from rdflib.namespace import XSD

# RDKit (ensure rdkit-pypi is installed)
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import MACCSkeys

# Bio.Entrez for NCBI taxonomy checks
from Bio import Entrez

# Database (use env vars; do NOT hardcode credentials)
import psycopg2
import psycopg2.extras


# ===========================
# ENV VARS (set these safely)
# ===========================
# Entrez email (required by NCBI)
ENTREZ_EMAIL = os.getenv("ENTREZ_EMAIL", "your_email@domain.com")

# COCONUT / OpenData Postgres (example)
PG_HOST = os.getenv("PG_HOST", "")
PG_PORT = int(os.getenv("PG_PORT", "5432"))
PG_DB   = os.getenv("PG_DB", "")
PG_USER = os.getenv("PG_USER", "")
PG_PASS = os.getenv("PG_PASS", "")
PG_SSLMODE = os.getenv("PG_SSLMODE", "prefer")
PG_GSSENC = os.getenv("PG_GSSENC", "disable")

# DrugCentral public DB
DC_HOST = os.getenv("DC_HOST", "unmtid-dbs.net")
DC_PORT = int(os.getenv("DC_PORT", "5433"))
DC_DB   = os.getenv("DC_DB", "drugcentral")
DC_USER = os.getenv("DC_USER", "drugman")
DC_PASS = os.getenv("DC_PASS", "dosage")

# Namespaces
PHYT = Namespace("http://www.semanticweb.org/orestah/ontologies/2024/9/phytotherapies#")


# ============
# LOGGING SETUP
# ============
logging.basicConfig(
    level=os.getenv("LOGLEVEL", "INFO"),
    format="%(asctime)s | %(levelname)s | %(message)s",
)
logger = logging.getLogger("phyto-ontology")


# ============
# UTILITIES
# ============
def uri_safe(s: str) -> str:
    s = re.sub(r"\W+", "_", str(s))
    s = re.sub(r"_+", "_", s)
    return s.strip("_")


def sanitize_for_uri(text: str) -> str:
    # Keep aligned with previous behavior
    return re.sub(r"\W|^(?=\d)", "_", text or "")


def normalize_scientific_name(name: str) -> str:
    name = (name or "").strip()
    parts = name.split()
    if len(parts) >= 2:
        parts[0] = parts[0].capitalize()
        parts[1] = parts[1].lower()
        return f"{parts[0]} {parts[1]}"
    elif len(parts) == 1:
        return parts[0].capitalize()
    return ""


def requests_get_json(url: str, timeout: float = 10.0) -> Optional[dict]:
    try:
        r = requests.get(url, timeout=timeout)
        if r.status_code == 200:
            return r.json()
    except requests.RequestException as e:
        logger.debug(f"Request failed: {e}")
    return None


def connect_pg(
    host: str, port: int, dbname: str, user: str, password: str, sslmode: str, gssencmode: str
) -> psycopg2.extensions.connection:
    return psycopg2.connect(
        host=host,
        port=port,
        dbname=dbname,
        user=user,
        password=password,
        sslmode=sslmode,
        gssencmode=gssencmode,
    )


# ============
# VALIDATION
# ============
def is_powo_plant(scientific_name: str) -> bool:
    query = scientific_name.strip()
    if not query:
        return False
    url = f"https://powo.science.kew.org/api/2/search?q={query.replace(' ', '%20')}"
    data = requests_get_json(url, timeout=10)
    if not data:
        return False
    results = data.get("results") or []
    if not results:
        return False
    # Conservative: consider a match if returned
    for result in results:
        if str(result.get("name", "")).lower() == query.lower():
            return True
    return True


def is_plantae_ncbi(scientific_name: str) -> bool:
    Entrez.email = ENTREZ_EMAIL
    try:
        search = Entrez.esearch(db="taxonomy", term=scientific_name)
        res = Entrez.read(search)
        search.close()
        if not res.get("IdList"):
            return False
        tax_id = res["IdList"][0]
        fetch = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        records = Entrez.read(fetch)
        fetch.close()
        rec = records[0]
        for rank in rec.get("LineageEx", []):
            if rank.get("Rank", "").lower() == "kingdom":
                sci = rank.get("ScientificName", "").lower()
                if "plantae" in sci or "viridiplantae" in sci:
                    return True
        lineage = rec.get("Lineage", "").lower()
        return ("plantae" in lineage) or ("viridiplantae" in lineage)
    except Exception as e:
        logger.debug(f"NCBI taxonomy error: {e}")
        return False


# ===========================
# PHENOTYPE CLASSIFICATION
# ===========================
PHENOTYPE_CLASSES = {
    "Carotenoid": PHYT.Carotenoid,
    "DietaryFiber": PHYT.DietaryFiber,
    "Isoprenoid": PHYT.Isoprenoid,
    "Phytosterol": PHYT.Phytosterol,
    "Polyphenol": PHYT.Polyphenol,
    "Polysaccharide": PHYT.Polysaccharide,
    "Saponin": PHYT.Saponin,
    "Unknown": PHYT.Unknown,
}


def classify_smiles(smiles: str) -> List[str]:
    mol = Chem.MolFromSmiles(smiles or "")
    if mol is None:
        return ["Unknown"]
    classes: List[str] = []
    mw = Descriptors.MolWt(mol)
    num_rings = Descriptors.RingCount(mol)
    num_oh = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
    num_c = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
    num_o = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "O")
    num_double = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))

    if num_rings >= 2 and num_oh >= 2:
        classes.append("Polyphenol")
    if num_c >= 30 and num_double >= 5:
        classes.append("Carotenoid")
    if num_rings >= 4 and num_oh >= 1:
        classes.append("Phytosterol")
    if num_c >= 40 and num_o >= 5:
        classes.append("Saponin")
    if mw > 500 and num_o > 10:
        classes.append("Polysaccharide")
    if num_c > 20 and num_o > 15:
        classes.append("DietaryFiber")
    if "C=C" in smiles and (num_c % 5 == 0):
        classes.append("Isoprenoid")

    return classes or ["Unknown"]


# ===========================
# ONTOLOGY HELPERS
# ===========================
def ensure_class(g: Graph, cls: URIRef, parent: Optional[URIRef] = None, label: Optional[str] = None) -> None:
    if (cls, RDF.type, OWL.Class) not in g:
        g.add((cls, RDF.type, OWL.Class))
        if parent:
            g.add((cls, RDFS.subClassOf, parent))
        if label:
            g.add((cls, RDFS.label, Literal(label)))


def add_literal(g: Graph, s: URIRef, p: URIRef, v: Optional[object]) -> None:
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return
    g.add((s, p, Literal(v)))


# ===========================
# SUBCOMMAND IMPLEMENTATIONS
# ===========================
def cmd_plants_from_coconut(args: argparse.Namespace) -> None:
    """
    Pull organisms from COCONUT DB, validate as plants (NCBI -> POWO fallback),
    and add as instances of PHYT.Plant in the ontology.
    """
    # Connect COCONUT
    with connect_pg(PG_HOST, PG_PORT, PG_DB, PG_USER, PG_PASS, PG_SSLMODE, PG_GSSENC) as conn:
        cur = conn.cursor()
        cur.execute("SELECT id, name FROM coconut.organisms;")
        organisms = cur.fetchall()

    g = Graph()
    g.parse(args.in_file, format=args.in_format)

    ensure_class(g, PHYT.Plant, parent=PHYT.PlantTaxonomy if (PHYT.PlantTaxonomy, RDF.type, OWL.Class) in g else None, label="Plant")

    cache_path = args.cache or "plant_validation_cache.pkl"
    try:
        validation_cache: Dict[str, bool] = pd.read_pickle(cache_path)  # type: ignore
    except Exception:
        validation_cache = {}

    added = 0
    for org_id, name in tqdm(organisms, desc="Validating organisms"):
        if not name:
            continue
        clean = normalize_scientific_name(name)
        if not clean:
            continue
        if clean in validation_cache:
            ok = validation_cache[clean]
        else:
            ok = is_plantae_ncbi(clean) or is_powo_plant(clean)
            validation_cache[clean] = ok
            time.sleep(0.2)  # API courtesy

        if ok:
            org_uri = PHYT[f"Organism_{org_id}"]
            g.add((org_uri, RDF.type, OWL.NamedIndividual))
            g.add((org_uri, RDF.type, PHYT.Plant))
            g.add((org_uri, RDFS.label, Literal(name)))
            added += 1

    logger.info("Plants added: %d", added)
    pd.to_pickle(validation_cache, cache_path)
    g.serialize(args.out_file, format=args.out_format)


def cmd_add_chemicals(args: argparse.Namespace) -> None:
    """
    Adds chemicals + data from COCONUT + CMAUP/NPASS merges; links chemicals to plants via isDerivedFrom.
    Requires CSVs from CMAUP directory (see README in repo).
    """
    g = Graph()
    g.parse(args.in_file, format=args.in_format)

    # Load CMAUP/NPASS merged CSV (expects columns like Plant_Name, SMILES, MW, Molecular_Formula, iupac_name, pref_name, Phytochemical_Class)
    df = pd.read_csv(args.cma_merged_csv)
    df = df.dropna(subset=["SMILES", "Phytochemical_Class", "Plant_Name"]).copy()
    df["SMILES"] = df["SMILES"].astype(str).str.strip()
    df["Plant_Name"] = df["Plant_Name"].astype(str).str.strip()

    hasMW = PHYT.hasMolecularWeight
    hasMF = PHYT.hasMolecularFormula
    hasCN = PHYT.hasCommonName
    hasIUPAC = PHYT.hasIUPACName
    hasSMI = PHYT.hasSMILES
    isDerivedFrom = PHYT.isDerivedFrom

    # Ensure Plant instances exist for Plant_Name (create if missing)
    existing_labels_lower = {
        str(o).strip().lower()
        for s, p, o in g.triples((None, RDF.type, PHYT.Plant))
        for _, _, lab in g.triples((s, RDFS.label, None))
        if (o == PHYT.Plant)
    }

    for plant in df["Plant_Name"].dropna().drop_duplicates():
        if plant.strip().lower() in existing_labels_lower:
            continue
        plant_uri = PHYT[f"Plant_{uri_safe(plant)}"]
        g.add((plant_uri, RDF.type, OWL.NamedIndividual))
        g.add((plant_uri, RDF.type, PHYT.Plant))
        g.add((plant_uri, RDFS.label, Literal(plant)))

    existing_smiles = {str(o) for s, p, o in g.triples((None, hasSMI, None))}
    added_chems = 0
    links = 0

    for _, row in df.iterrows():
        smi = row["SMILES"]
        chem_uri = PHYT[f"Chemical_{abs(hash(smi))}"]
        if smi not in existing_smiles:
            g.add((chem_uri, RDF.type, OWL.NamedIndividual))
            for cls in str(row["Phytochemical_Class"]).split(";"):
                cls_uri = PHENOTYPE_CLASSES.get(cls.strip(), PHYT.Unknown)
                g.add((chem_uri, RDF.type, cls_uri))
            g.add((chem_uri, hasSMI, Literal(smi)))
            g.add((chem_uri, RDFS.label, Literal(smi)))
            add_literal(g, chem_uri, hasMW, pd.to_numeric(row.get("MW"), errors="coerce"))
            add_literal(g, chem_uri, hasMF, row.get("Molecular_Formula"))
            add_literal(g, chem_uri, hasIUPAC, row.get("iupac_name"))
            add_literal(g, chem_uri, hasCN, row.get("pref_name"))
            existing_smiles.add(smi)
            added_chems += 1

        plant_uri = PHYT[f"Plant_{uri_safe(row['Plant_Name'])}"]
        g.add((chem_uri, isDerivedFrom, plant_uri))
        links += 1

    logger.info("Chemicals added: %d | Links added: %d", added_chems, links)
    g.serialize(args.out_file, format=args.out_format)


def cmd_add_maccs(args: argparse.Namespace) -> None:
    g = Graph()
    g.parse(args.in_file, format=args.in_format)

    hasSMI = PHYT.hasSMILES
    hasMACCs = PHYT.hasMACCs
    chemical_classes = set(PHENOTYPE_CLASSES.values())

    added = 0
    for s, p, o in g.triples((None, RDF.type, None)):
        if o not in chemical_classes:
            continue
        # Skip if MACCS already present
        if (s, hasMACCs, None) in g:
            continue
        smi = None
        for _, _, val in g.triples((s, hasSMI, None)):
            smi = str(val).strip()
            break
        if not smi:
            continue
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        fp = MACCSkeys.GenMACCSKeys(mol).ToBitString()
        g.add((s, hasMACCs, Literal(fp, datatype=XSD.string)))
        added += 1

    logger.info("MACCS keys added to %d chemical instances.", added)
    g.serialize(args.out_file, format=args.out_format)


def _drugcentral_query(inchikeys: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Returns (action_type_df, target_df, atc_df)
    """
    if not inchikeys:
        return (pd.DataFrame(), pd.DataFrame(), pd.DataFrame())

    with psycopg2.connect(host=DC_HOST, port=DC_PORT, dbname=DC_DB, user=DC_USER, password=DC_PASS) as conn:
        action_type_query = """
            SELECT s.inchikey, atf.action_type
            FROM structures s
            JOIN act_table_full atf ON s.id = atf.struct_id
            WHERE s.inchikey = ANY(%s) AND atf.action_type IS NOT NULL
        """
        target_query = """
            SELECT s.inchikey, td.name AS target
            FROM structures s
            JOIN act_table_full atf ON s.id = atf.struct_id
            JOIN target_dictionary td ON atf.target_id = td.id
            WHERE s.inchikey = ANY(%s)
        """
        atc_query = """
            SELECT s.inchikey, a.code AS atc_code
            FROM structures s
            JOIN struct2atc s2a ON s.id = s2a.struct_id
            JOIN atc a ON s2a.atc_code = a.code
            WHERE s.inchikey = ANY(%s)
        """
        action_type_df = pd.read_sql(action_type_query, conn, params=(inchikeys,))
        target_df = pd.read_sql(target_query, conn, params=(inchikeys,))
        atc_df = pd.read_sql(atc_query, conn, params=(inchikeys,))
    return action_type_df, target_df, atc_df


def cmd_drugcentral_enrich(args: argparse.Namespace) -> None:
    g = Graph()
    g.parse(args.in_file, format=args.in_format)

    hasSMI = PHYT.hasSMILES
    hasInChIKey = PHYT.hasInChIKey
    TherapeuticEffect = PHYT.TherapeuticEffect
    TargetedPathway = PHYT.TargetedPathway
    hasTherapeuticEffect = PHYT.hasTherapeuticEffect
    targetsPathway = PHYT.targetsPathway
    hasValue = PHYT.hasValue

    # Build SMILES → subject and InChIKey mappings
    smiles_to_subject = {}
    for subj, _, obj in g.triples((None, hasSMI, None)):
        smiles_to_subject[str(obj)] = subj

    def to_inchikey(smiles: str) -> Optional[str]:
        m = Chem.MolFromSmiles(smiles)
        return Chem.MolToInchiKey(m) if m else None

    smiles_inchikey = {s: to_inchikey(s) for s in tqdm(smiles_to_subject.keys(), desc="SMILES→InChIKey")}
    # store inchikeys back to ontology
    for smi, subj in smiles_to_subject.items():
        ik = smiles_inchikey.get(smi)
        if ik:
            g.add((subj, hasInChIKey, Literal(ik)))

    inchikeys = [ik for ik in smiles_inchikey.values() if ik]
    action_df, target_df, atc_df = _drugcentral_query(inchikeys)

    # Enrich ontology
    effect_instances: Dict[str, URIRef] = {}
    target_instances: Dict[str, URIRef] = {}

    # action types → TherapeuticEffect
    for _, row in action_df.iterrows():
        ik = row["inchikey"]
        effect = str(row["action_type"]).strip()
        if not effect:
            continue
        # find subject
        subj = None
        for s, _, o in g.triples((None, hasInChIKey, Literal(ik))):
            subj = s
            break
        if not subj:
            continue
        eff_uri = effect_instances.get(effect)
        if not eff_uri:
            eff_uri = PHYT[uri_safe(effect)]
            g.add((eff_uri, RDF.type, TherapeuticEffect))
            g.add((eff_uri, hasValue, Literal(effect)))
            effect_instances[effect] = eff_uri
        g.add((subj, hasTherapeuticEffect, eff_uri))

    # targets → TargetedPathway
    for _, row in target_df.iterrows():
        ik = row["inchikey"]
        target = str(row["target"]).strip()
        if not target:
            continue
        subj = None
        for s, _, o in g.triples((None, hasInChIKey, Literal(ik))):
            subj = s
            break
        if not subj:
            continue
        t_uri = target_instances.get(target)
        if not t_uri:
            t_uri = PHYT[uri_safe(target)]
            g.add((t_uri, RDF.type, TargetedPathway))
            g.add((t_uri, hasValue, Literal(target)))
            target_instances[target] = t_uri
        g.add((subj, targetsPathway, t_uri))

    # Optionally attach ATC to enriched subjects
    hasATC = PHYT.hasATC
    # Only add ATC to subjects that have either effect or target defined
    enriched_subjects = {
        s for s, _, _ in g.triples((None, hasTherapeuticEffect, None))
    } | {s for s, _, _ in g.triples((None, targetsPathway, None))}

    # inchikey → subject map (enriched only)
    ik_to_subj: Dict[str, URIRef] = {}
    for s in enriched_subjects:
        for _, _, o in g.triples((s, hasInChIKey, None)):
            ik_to_subj[str(o)] = s

    added_atc = 0
    for _, row in atc_df.iterrows():
        ik = row["inchikey"]
        atc = row["atc_code"]
        subj = ik_to_subj.get(ik)
        if subj and atc:
            g.add((subj, hasATC, Literal(atc)))
            added_atc += 1

    logger.info("DrugCentral enrichment complete. Effects: %d, Targets: %d, ATC added: %d",
                len(action_df), len(target_df), added_atc)
    g.serialize(args.out_file, format=args.out_format)


def cmd_chembl_mechanisms(args: argparse.Namespace) -> None:
    """
    Adds TargetedPathway + TherapeuticEffect from a ChEMBL mechanisms CSV
    (columns required: 'Smiles', 'Target Name', 'Action Type').
    """
    g = Graph()
    g.parse(args.in_file, format=args.in_format)

    df = pd.read_csv(args.chembl_csv, sep=";") if args.chembl_csv.endswith(".csv") else pd.read_csv(args.chembl_csv)
    df["Smiles"] = df["Smiles"].astype(str).str.strip()

    hasSMI = PHYT.hasSMILES
    TherapeuticEffect = PHYT.TherapeuticEffect
    TargetedPathway = PHYT.TargetedPathway
    hasTherapeuticEffect = PHYT.hasTherapeuticEffect
    targetsPathway = PHYT.targetsPathway
    hasValue = PHYT.hasValue

    # Canonical mapping for ontology SMILES
    smi_to_subj: Dict[str, URIRef] = {}
    for subj, _, obj in g.triples((None, hasSMI, None)):
        smi_to_subj[str(obj)] = subj

    def canonical(s: str) -> Optional[str]:
        m = Chem.MolFromSmiles(s)
        return Chem.MolToSmiles(m, canonical=True) if m else None

    can_to_subj = {}
    for smi, subj in smi_to_subj.items():
        c = canonical(smi)
        if c:
            can_to_subj[c] = subj

    added_target = 0
    added_effect = 0
    pathway_instances: Dict[str, URIRef] = {}
    effect_instances: Dict[str, URIRef] = {}

    df["Smiles_canonical"] = df["Smiles"].apply(canonical)
    for _, row in df.iterrows():
        c = row["Smiles_canonical"]
        if not c or c not in can_to_subj:
            continue
        subj = can_to_subj[c]
        tname = str(row.get("Target Name", "")).strip()
        ename = str(row.get("Action Type", "")).strip()

        if tname:
            tu = pathway_instances.get(tname)
            if not tu:
                tu = PHYT[uri_safe(tname)]
                g.add((tu, RDF.type, TargetedPathway))
                g.add((tu, hasValue, Literal(tname)))
                pathway_instances[tname] = tu
            g.add((subj, targetsPathway, tu))
            added_target += 1

        if ename:
            eu = effect_instances.get(ename)
            if not eu:
                eu = PHYT[uri_safe(ename)]
                g.add((eu, RDF.type, TherapeuticEffect))
                g.add((eu, hasValue, Literal(ename)))
                effect_instances[ename] = eu
            g.add((subj, hasTherapeuticEffect, eu))
            added_effect += 1

    logger.info("ChEMBL mechanisms added. targets=%d, effects=%d", added_target, added_effect)
    g.serialize(args.out_file, format=args.out_format)


def cmd_export_links(args: argparse.Namespace) -> None:
    g = Graph()
    g.parse(args.in_file, format=args.in_format)

    isDerivedFrom = PHYT.isDerivedFrom
    hasSMILES = PHYT.hasSMILES

    def get_smiles(s: URIRef) -> str:
        for _, _, v in g.triples((s, hasSMILES, None)):
            if isinstance(v, Literal):
                return str(v)
        return "N/A"

    def plant_label(puri: URIRef) -> str:
        for _, _, lab in g.triples((puri, RDFS.label, None)):
            if isinstance(lab, Literal):
                return str(lab)
        return str(puri).split("#")[-1]

    rows: List[Tuple[str, str]] = []
    for chem, _, plant in g.triples((None, isDerivedFrom, None)):
        rows.append((get_smiles(chem), plant_label(plant)))

    out = args.csv
    with open(out, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Chemical (SMILES)", "Plant (Label)"])
        writer.writerows(rows)
    logger.info("Exported %d rows to %s", len(rows), out)


def cmd_stats(args: argparse.Namespace) -> None:
    g = Graph()
    g.parse(args.in_file, format=args.in_format)

    total_triples = len(g)
    classes: Set[URIRef] = set()
    individuals: Set[URIRef] = set()
    object_props: Set[URIRef] = set()
    data_props: Set[URIRef] = set()

    instances_by_class: Dict[URIRef, int] = defaultdict(int)

    for s, p, o in g:
        if p == RDF.type:
            if o in (OWL.Class, RDFS.Class):
                classes.add(s)
            elif o not in (OWL.NamedIndividual, RDF.Property):
                instances_by_class[o] += 1
                individuals.add(s)

    # Count properties by usage
    for s, p, o in g:
        if isinstance(o, Literal):
            if p != RDF.type:
                data_props.add(p)
        else:
            if p != RDF.type:
                object_props.add(p)

    def label(u: URIRef) -> str:
        for _, _, lab in g.triples((u, RDFS.label, None)):
            return str(lab)
        return str(u).split("#")[-1]

    logger.info("Total RDF Triples: %s", total_triples)
    logger.info("Total Classes: %s", len(classes))
    logger.info("Total Object Properties: %s", len(object_props))
    logger.info("Total Data Properties: %s", len(data_props))
    logger.info("Total Individuals: %s", len(individuals))
    logger.info("Instance counts by class:")
    for cls, cnt in sorted(instances_by_class.items(), key=lambda x: x[1], reverse=True)[:30]:
        logger.info("  %s: %d", label(cls), cnt)


def cmd_add_ecfp(args: argparse.Namespace) -> None:
    """
    Add ECFP strings from CSV to individuals matched by exact SMILES via hasSmiles.
    CSV must have columns: smiles, ecfp
    """
    g = Graph()
    g.parse(args.in_file, format=args.in_format)

    df = pd.read_csv(args.ecfp_csv)
    if not {"smiles", "ecfp"}.issubset(df.columns.str.lower()):
        raise ValueError("ECFP CSV must include 'smiles' and 'ecfp' columns.")

    # Normalize column names
    col_map = {c.lower(): c for c in df.columns}
    smi_col = col_map["smiles"]
    ecfp_col = col_map["ecfp"]

    hasSmiles = PHYT.hasSmiles  # Note: your ontology may use hasSMILES vs hasSmiles (align to your schema)
    hasECFP = PHYT.hasECFP

    added = 0
    for _, row in df.iterrows():
        smi = str(row[smi_col])
        ecfp = str(row[ecfp_col])
        for subj in g.subjects(predicate=hasSmiles, object=Literal(smi, datatype=XSD.string)):
            g.add((subj, hasECFP, Literal(ecfp, datatype=XSD.string)))
            added += 1

    logger.info("ECFP properties added: %d", added)
    g.serialize(args.out_file, format=args.out_format)


def cmd_heatmap(args: argparse.Namespace) -> None:
    """
    Build a top-50x50 heatmap from CTA predictions CSV with columns:
      ['np_id','uniprot','mean_score']
    Saves a PNG figure.
    """
    import matplotlib.pyplot as plt
    try:
        import seaborn as sns  # optional
        use_seaborn = True
    except ImportError:
        use_seaborn = False

    df = pd.read_csv(args.pred_csv)
    if not {"np_id", "uniprot", "mean_score"}.issubset(df.columns):
        raise ValueError("Prediction CSV must include columns: np_id, uniprot, mean_score")

    top_chems = df["np_id"].value_counts().head(50).index
    top_prots = df["uniprot"].value_counts().head(50).index
    filt = df[df["np_id"].isin(top_chems) & df["uniprot"].isin(top_prots)]
    mat = filt.pivot_table(index="np_id", columns="uniprot", values="mean_score", aggfunc="mean").fillna(0)

    plt.figure(figsize=(16, 12))
    if use_seaborn:
        sns.heatmap(mat, cmap="viridis", linewidths=0.1, linecolor="gray")
    else:
        # Minimal matplotlib fallback
        plt.imshow(mat.values, aspect="auto")
        plt.colorbar(label="Mean Score")
        plt.xticks(range(len(mat.columns)), mat.columns, rotation=90, fontsize=6)
        plt.yticks(range(len(mat.index)), mat.index, fontsize=6)
    plt.title("Phytochemical vs Protein Target Predictions (Top 50×50)")
    plt.xlabel("Protein Target (UniProt)")
    plt.ylabel("Phytochemical (InChIKey)")
    plt.tight_layout()
    plt.savefig(args.out_png, dpi=300)
    logger.info("Saved heatmap: %s", args.out_png)


# ============
# CLI
# ============
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Phytotherapy Ontology Pipeline (PSB-ready)")

    sub = p.add_subparsers(dest="cmd", required=True)

    # plants-from-coconut
    sp = sub.add_parser("plants-from-coconut", help="Validate organisms as plants and add to ontology")
    sp.add_argument("--in", dest="in_file", required=True, help="Input ontology file")
    sp.add_argument("--in-format", dest="in_format", default="xml", help="rdflib format (xml,turtle,...)")
    sp.add_argument("--out", dest="out_file", required=True, help="Output ontology file")
    sp.add_argument("--out-format", dest="out_format", default="xml")
    sp.add_argument("--cache", dest="cache", default="plant_validation_cache.pkl")
    sp.set_defaults(func=cmd_plants_from_coconut)

    # add-chemicals
    sp = sub.add_parser("add-chemicals", help="Add chemicals & link to plants from CMAUP/NPASS merges")
    sp.add_argument("--in", dest="in_file", required=True)
    sp.add_argument("--in-format", dest="in_format", default="xml")
    sp.add_argument("--cma-merged-csv", dest="cma_merged_csv", required=True,
                    help="CSV with Plant_Name, SMILES, MW, Molecular_Formula, iupac_name, pref_name, Phytochemical_Class")
    sp.add_argument("--out", dest="out_file", required=True)
    sp.add_argument("--out-format", dest="out_format", default="xml")
    sp.set_defaults(func=cmd_add_chemicals)

    # add-maccs
    sp = sub.add_parser("add-maccs", help="Compute and add MACCS keys to chemical instances")
    sp.add_argument("--in", dest="in_file", required=True)
    sp.add_argument("--in-format", dest="in_format", default="xml")
    sp.add_argument("--out", dest="out_file", required=True)
    sp.add_argument("--out-format", dest="out_format", default="xml")
    sp.set_defaults(func=cmd_add_maccs)

    # drugcentral-enrich
    sp = sub.add_parser("drugcentral-enrich", help="Enrich ontology with DrugCentral targets, action types, and ATC")
    sp.add_argument("--in", dest="in_file", required=True)
    sp.add_argument("--in-format", dest="in_format", default="xml")
    sp.add_argument("--out", dest="out_file", required=True)
    sp.add_argument("--out-format", dest="out_format", default="xml")
    sp.set_defaults(func=cmd_drugcentral_enrich)

    # chembl-mechanisms
    sp = sub.add_parser("chembl-mechanisms", help="Add targets/effects from ChEMBL mechanisms CSV")
    sp.add_argument("--in", dest="in_file", required=True)
    sp.add_argument("--in-format", dest="in_format", default="xml")
    sp.add_argument("--chembl-csv", dest="chembl_csv", required=True)
    sp.add_argument("--out", dest="out_file", required=True)
    sp.add_argument("--out-format", dest="out_format", default="xml")
    sp.set_defaults(func=cmd_chembl_mechanisms)

    # export-links
    sp = sub.add_parser("export-links", help="Export SMILES → Plant label mappings to CSV")
    sp.add_argument("--in", dest="in_file", required=True)
    sp.add_argument("--in-format", dest="in_format", default="xml")
    sp.add_argument("--csv", dest="csv", required=True)
    sp.set_defaults(func=cmd_export_links)

    # stats
    sp = sub.add_parser("stats", help="Compute ontology statistics (triples, classes, instances)")
    sp.add_argument("--in", dest="in_file", required=True)
    sp.add_argument("--in-format", dest="in_format", default="turtle")
    sp.set_defaults(func=cmd_stats)

    # add-ecfp
    sp = sub.add_parser("add-ecfp", help="Attach ECFP values from CSV to instances by matching SMILES")
    sp.add_argument("--in", dest="in_file", required=True)
    sp.add_argument("--in-format", dest="in_format", default="turtle")
    sp.add_argument("--ecfp-csv", dest="ecfp_csv", required=True)
    sp.add_argument("--out", dest="out_file", required=True)
    sp.add_argument("--out-format", dest="out_format", default="turtle")
    sp.set_defaults(func=cmd_add_ecfp)

    # heatmap
    sp = sub.add_parser("heatmap", help="Render a top-50×50 heatmap from CTA predictions")
    sp.add_argument("--pred", dest="pred_csv", required=True, help="cta_predictions_all_phytotherapies.csv")
    sp.add_argument("--out", dest="out_png", required=True, help="Output PNG path")
    sp.set_defaults(func=cmd_heatmap)

    return p


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    # sanity check critical envs
    if not ENTREZ_EMAIL:
        logger.warning("ENTREZ_EMAIL is not set; NCBI validation may fail.")
    try:
        args.func(args)
        return 0
    except Exception as e:
        logger.exception("Error: %s", e)
        return 1


if __name__ == "__main__":
    sys.exit(main())
