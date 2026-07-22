# Virtual Screen Tools
This repository contains virtual screen tools meant to be used with [SPRINT](https://github.com/abhinadduri/panspecies-dti) and [GNINA](https://github.com/gnina/gnina). This repository includes scripts for selecting unique drug-like compounds from Deepdock screens as well as scripts to visualize drug-target interaction predictions generated with SPRINT.

## Forward Screen
For screening virtual chemical libraries against a specific target.

**deepdock_helper_template.slurm**
This script runs deepdock on high performance computing clusters using the slurm scheduler. Make sure to edit values to including partition, email, and GNINA and SPRINT checkpoints to your needs.
```bash
sbatch deepdock_helper_template.slurm all library.parquet pdhe1/receptor.pdb pdhe1/ligand.pdb pdhe1/sequence.pdb pdhe1_screen
```

**drug_like_filter.py** filters Deepdock output predictions (topn.sdf) for drug-like properties including maximum molecular weight, predicted affinity (from GNINA), and predicted solubility (from [SolTranNet](https://github.com/gnina/SolTranNet)).
```bash
python drug_like_fitler.py --topn pdhe1/topn.sdf --outPath pdhe1/filtered_compounds.csv
```

**select_unique.py** selects the chemically unique compounds out of a list of SMILES using RDKIT. Meant to be used on drug_like_filter.py output.
```bash
python drug_like_fitler.py --drugFile pdhe1/filtered_compounds.csv --outPath pdhe1/selected_compounds.csv
```

## Reverse Screen
For screening known inhibitors against a proteome of interest.

**dti_analysis.py** generates a principal component analysis (PCA) of SPRINT compound embeddings, then underlays proteome SPRINT embeddings onto the graph to visualize drug-target interaction (DTI) predictions. It also creates a UMAP using the same methodology. It labels compound outliers from the PCA across both graphs, and can generate and save outlier interaction predictions.
```bash
python dti_analysis.py all  \
    --drugSmiles drugs.csv           \
    --drugEmbeds drug_embeds.npy           \
    --proteome proteome.csv             \
    --receptorEmbeds proteome_embeds.npy       \
    --interestEmbeds forward_screen_targets.npy \
    --interestLabel "Forward Screen Targets" \
    --plotPath dti_analysis.png \
    --outlierDir outlier_target_predictions
```

![dti analysis plot][dti_analysis]

**outlier_chemical_analysis** generates the same PCA of SPRINT compound embeddings as dti_analysis.py but also creates a PCA of the molecular fingerprints of each compound. Outliers from the SPRINT PCA are labeled in the fingerprint PCA, creating an easy visualization that shows whether SPRINT embedding outliers are also chemical outliers. Refer to dti_analysis.py for usage.
![chemical analysis plot][chemical_analysis]

[dti_analysis]: https://github.com/CFerrance/Virtual-Screen-Tools/tree/main/Images/dti_analysis_example.png "DTI Analysis"
[chemical_analysis]: https://github.com/CFerrance/Virtual-Screen-Tools/tree/main/Images/chemical_analysis_example.png "Chemical Analysis"