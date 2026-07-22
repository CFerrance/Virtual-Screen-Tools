import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


def is_not_similar(compundFingerprint, selectedFingerprints, cutoff) -> bool:
    for fp in selectedFingerprints:
        if DataStructs.TanimotoSimilarity(compundFingerprint, fp) > cutoff:
            return False
    return True


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Deepdock Hit Selection")
    parser.add_argument("--drugFile",  help=".csv with ID, SMILES, and CNN_VS", required=True)
    parser.add_argument("--similarityCutoff", help="Maximum allowed similarity", default=0.7)
    parser.add_argument("--outPath",  help="File output path (.txt?)", required=True)
    
    args = parser.parse_args()

    df = pd.read_csv(args.drugFile, header = 0)
    df = df.sort_values(by='CNN_VS', ascending=False)

    selections = []
    fingerprints = []
    fpgen = AllChem.GetMorganGenerator(radius=2, fpSize = 2048)

    for entry in df.itertuples(index=False):
        mol = Chem.MolFromSmiles(entry.SMILES)
        fp = fpgen.GetFingerprint(mol)
        if is_not_similar(fp, fingerprints, args.similarityCutoff):
            selections.append(entry)
            fingerprints.append(fp)
    
    selectDF = pd.DataFrame(selections)

    print(str(len(selectDF)) + " of " + str(len(df)) + " molecules selected.")

    selectDF.to_csv(args.outPath, index= False)