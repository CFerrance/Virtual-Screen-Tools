import argparse
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import soltrannet as stn


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filter for ligands based on molecular weight, affinity, and solubility")
    parser.add_argument("--topn",  help=".sdf file containing ligands", required=True)
    parser.add_argument("--maxWeight", help="Maximum weight", default=500)
    parser.add_argument("--minAffinity",  help="Minimum affinity", default=-7)
    parser.add_argument("--minSolubility",  help="Minimum solubility", default=-4)
    parser.add_argument("--outPath",  help=".csv out path", required=True)
    
    args = parser.parse_args()

    structures = Chem.SDMolSupplier(args.topn)
    molecules = []

    for mol in structures:
        #get weight
        try: 
            weight = Descriptors.MolWt(mol)
        except:
            print("Skipping Molecule Due To Mol Weight Error")
            continue
        
        #get affinity
        affinity = float(mol.GetProp("minimizedAffinity"))
        
        #get id and smiles from name
        name = mol.GetProp("_Name")
        spaceIndex = name.find(" ")
        assert spaceIndex != -1
        id = name[:spaceIndex]
        smiles = name[spaceIndex + 1:]
        
        #add CNN_VS and store molecule
        cnn = mol.GetProp("CNN_VS")
        molecules.append({"MoleculeID" : id, "SMILES" : smiles, "MolWeight" : weight, "CNN_VS" : cnn, "MinimizedAffinity" : affinity})

    filtered = pd.DataFrame(molecules) 

    #filter by weight and affinity
    filtered = filtered[filtered["MolWeight"] < args.maxWeight]
    filtered = filtered[filtered["MinimizedAffinity"] < args.minAffinity]

    #SolTranNet Predictions and Merge
    predictions = pd.DataFrame(list(stn.predict(filtered["SMILES"])), columns=["Solubility", "SMILES", "Log"])
    filtered = pd.merge(filtered, predictions, how="inner", on="SMILES")

    #filter by solubility
    filtered = filtered[filtered["Solubility"] > args.minSolubility]

    print(str(len(filtered)) + " of " + str(len(molecules)) + " selected.")

    #save
    filtered.to_csv(args.outPath, index=False)
