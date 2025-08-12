# Modified CTAPred Pipeline

CTAPred pipeline adapted for parallelization by splitting input dataset and processing each chunk on a separate core which massively reduces runtime for large data inputs. Also replaces some Pandas functionality with Polars operations for greater efficiency. 

## Usage
```
python predict_targets.py --input=InputFolderPath --output=OutputFolderPath [additional params]
```

`--input` should be the full path to directory containing input CSV file with name "QueryList1_smiles.csv". The input file should have a "smiles" column which contains SMILES strings and a "np_id" column for your choice of compound identifier.

### Advanced Usage

```
python predict_targets.py \
    --input my_compounds.csv \
    --output results/ \
    --fingerprint ecfp \ # fingerprint choice [avalon, ecfp, fcfp, maccs]
    --nBits 1024 \ # number of bits for ecfp/fcfp fingerprint
    --radius 3 \ # radius for ecfp fingerprint
    --k 1 5 10 \ # top-k reference compounds
    --n_jobs 8 \ # Number of cores on your machine
    --verbose
```

## References

This work built upon the original CTAPred pipeline:

```bibtex
@software{CTAPred,
  author       = {Alhasbary, Abeer Abdulhakeem Mansour},
  title        = {CTAPred},
  year         = {2024},
  publisher    = {GitHub},
  url          = {https://github.com/Alhasbary/CTAPred}
}
```

```bibtex
@article{ABDULHAKEEMMANSOURALHASBARY2025109351,
title = {Exploring natural products potential: A similarity-based target prediction tool for natural products},
journal = {Computers in Biology and Medicine},
volume = {184},
pages = {109351},
year = {2025},
issn = {0010-4825},
doi = {https://doi.org/10.1016/j.compbiomed.2024.109351},
url = {https://www.sciencedirect.com/science/article/pii/S0010482524014367},
author = {Abeer {Abdulhakeem Mansour Alhasbary} and Nurul {Hashimah Ahamed Hassain Malim} and Siti {Zuraidah Mohamad Zobir}},
keywords = {Similarity searching, Drug discovery, Compound-target interaction, Natural products, ChEMBL database},
abstract = {Natural products are invaluable resources in drug discovery due to their substantial structural diversity. However, predicting their interactions with druggable protein targets remains a challenge, primarily due to the limited availability of bioactivity data. This study introduces CTAPred (Compound-Target Activity Prediction), an open-source command-line tool designed to predict potential protein targets for natural products. CTAPred employs a two-stage approach, combining fingerprinting and similarity-based search techniques to identify likely drug targets for these bioactive compounds. Despite its simplicity, the tool's performance is comparable to that of more complex methods, demonstrating proficiency in target retrieval for natural product compounds. Furthermore, this study explores the optimal number of reference compounds most similar to the query compound, aiming to refine target prediction accuracy. The findings demonstrated the superior performance of considering only the most similar reference compounds for target prediction. CTAPred is freely available at https://github.com/Alhasbary/CTAPred, offering a valuable resource for deciphering natural product-target associations and advancing drug discovery.}
}
```