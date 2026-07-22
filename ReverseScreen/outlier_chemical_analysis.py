import argparse
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


def add_to_scatter(plot, data, color, desc, edge = False, order = 0):
    edgeColor = "none"
    if edge:
        edgeColor = 'k'
        
    plot.scatter(data[:, 0], data[:, 1], c = color, edgecolor = edgeColor, label = desc, zorder = order)


def ids_to_coordinates(ids, dict):
    coords = []
    for id in ids:
        coords.append(dict[id])
    return coords


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="DTI Embeddings vs Molecular Fingerprints")
    parser.add_argument("--drugSmiles",  help="Input file (.csv with id and SMILES columns)", required=True)
    parser.add_argument("--drugEmbeds",  help="Input file (SPRINT .npy file)", required=True)
    parser.add_argument("--proteome",  help="Input file (.csv with IDs and Gene Names)")
    parser.add_argument("--receptorEmbeds",  help="Input file (SPRINT .npy file)", required=True)
    parser.add_argument("--interestEmbeds",  help="Input file (SPRINT .npy file)")
    parser.add_argument("--interestLabel", default="Receptors of Interest",  help="Input file (SPRINT .npy file)")
    parser.add_argument("--xCutoff", default=5, help="PCA coordinate x cutoff")
    parser.add_argument("--yCutoff", default=5, help="PCA coordinate y cutoff")
    parser.add_argument("--plotPath",  help="Plot Output Path")
    args = parser.parse_args()

    drugEmbeddings = np.load(args.drugEmbeds)
    receptorEmbeddings = np.load(args.receptorEmbeds)

    interestEmbeddings = None
    if args.interestEmbeds != None:
        interestEmbeddings = np.load(args.interestEmbeds)

    #map drugs to embeddings
    smilesList = pd.read_csv(args.drugSmiles, header=0)["SMILES"].to_list()
    drugsDF = pd.read_csv(args.drugSmiles, header=0)
    indexToID = dict(zip(range(len(drugEmbeddings)), drugsDF["id"].to_list()))
    idToSmiles = dict(zip(drugsDF["id"].to_list(), drugsDF["SMILES"].to_list()))
    idsToEmbeddings = dict(zip(drugsDF["id"].to_list(), drugEmbeddings))

    #Initialize PCA and UMAP, fit to drugs
    pca = PCA(n_components=2)
    pca.fit(drugEmbeddings)

    #transform data
    pca_drugEmbeddings = pca.transform(drugEmbeddings)
    pca_receptorEmbeddings = pca.transform(receptorEmbeddings)

    #get PCA groups
    xCutoff = args.xCutoff
    yCutoff = args.yCutoff
    inliers = []
    xOutliers = []
    yOutliers = []
    xyOutliers = []
    idsToPCA = {}
    for i in range(len(pca_drugEmbeddings)):
        pcaEmbed = pca_drugEmbeddings[i]
        id = indexToID[i]
        idsToPCA[id] = pcaEmbed
        if pcaEmbed[0] < xCutoff and pcaEmbed[1] < yCutoff:
            inliers.append(id)
        elif pcaEmbed[0] >= xCutoff and pcaEmbed[1] < yCutoff:
            xOutliers.append(id)
        elif pcaEmbed[0] < xCutoff and pcaEmbed[1] >= yCutoff:
            yOutliers.append(id)
        elif pcaEmbed[0] >= xCutoff and pcaEmbed[1] >= yCutoff:
            xyOutliers.append(id)

    #create plots
    fig = plt.figure(figsize=(20,10))
    plotA = fig.add_subplot(121)
    plotB = fig.add_subplot(122)

    #create plotA: PCA of both Drug and Receptor Embeddings
    add_to_scatter(plotA, np.array(ids_to_coordinates(inliers, idsToPCA)), 'b', "Inliers", True, 3)
    add_to_scatter(plotA, np.array(ids_to_coordinates(xOutliers, idsToPCA)), 'c', "X Outliers", True, 3)
    add_to_scatter(plotA, np.array(ids_to_coordinates(yOutliers, idsToPCA)), 'y', "Y Outliers", True, 3)
    if len(xyOutliers) > 0:
        add_to_scatter(plotA, np.array(ids_to_coordinates(xyOutliers, idsToPCA)), "blueviolet", True, 3)
    
    add_to_scatter(plotA, pca_receptorEmbeddings, "pink", str(len(receptorEmbeddings)) + " Proteome Proteins", False, 1)

    if interestEmbeddings is not None:
        pca_interestEmbeddings = pca.transform(interestEmbeddings)
        add_to_scatter(plotA, pca_interestEmbeddings, 'r', args.interestLabel, False, 2)
    
    plotA.set_xlabel("PCA1: Exp. Var = " + str(round(pca.explained_variance_ratio_[0], 4)), fontsize = 14)
    plotA.set_ylabel("PCA2: Exp. Var = " + str(round(pca.explained_variance_ratio_[1], 4)), fontsize = 14)
    plotA.set_title("PCA of Inhibitor and Receptor SPRINT Embeddings", fontsize = 16)
    plotA.legend()

    #plotB, chemical analysis
    fpgen = AllChem.GetMorganGenerator(radius=2, fpSize = 2048)

    idToFingerprint = {}
    fingerprints = []
    for id, smiles in idToSmiles.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol == None:
            print("Fingerprint failed for " + smiles)
        fp = fpgen.GetFingerprint(mol)
        arr = np.zeros((2048,), dtype=bool)
        DataStructs.ConvertToNumpyArray(fp, arr)
        idToFingerprint[id] = fp
        fingerprints.append(fp)

    #pca of molecular fingerprints colored by PCA outlier status
    pca_fps = PCA(n_components=2)
    pca_fps.fit(fingerprints)

    def EmbedGroupToFingerprintPCA(group, color, description):
        fps = []
        for id in group:
            fps.append(idToFingerprint[id])
        
        pca_group = pca_fps.transform(fps)
        plotB.scatter(pca_group[:, 0], pca_group[:, 1], c = color, edgecolor = 'k', label = description)
        
    EmbedGroupToFingerprintPCA(inliers, 'b', "Embed Inliers")
    EmbedGroupToFingerprintPCA(xOutliers, 'c', "Embed X Outliers")
    EmbedGroupToFingerprintPCA(yOutliers, 'y', "Embed Y Outliers")
    if len(xyOutliers) > 0:
        EmbedGroupToFingerprintPCA(xyOutliers, "blueviolet", "Embed XY Outliers")

    #pca labels
    plotB.set_xlabel("PCA1: Exp. Var = " + str(round(pca_fps.explained_variance_ratio_[0], 4)), fontsize = 14)
    plotB.set_ylabel("PCA2: Exp. Var = " + str(round(pca_fps.explained_variance_ratio_[1], 4)), fontsize = 14)
    plotB.set_title("PCA of Molecular Fingerprints Labeled with SPRINT Embedding Outliers", fontsize = 16)
    plotB.legend()

    plt.savefig(args.plotPath)