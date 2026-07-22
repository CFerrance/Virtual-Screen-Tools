import os
import argparse
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import umap.umap_ as umap
import matplotlib.pyplot as plt


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


def write_outlier_essentials(outliers, idToSmiles, outPath):
    smiles = [idToSmiles[id] for id in outliers]
    df = pd.DataFrame({'IDs': outliers, 'SMILES': smiles})
    df.to_csv(outPath, index = False)


def write_outlier_neighbors(outliers, indexToProtein, outPath):
    distances, indices = nearest.kneighbors(outliers)
    outlierUniprotIds = set()
    for sub in indices:
        for index in sub:
            outlierUniprotIds.add(indexToProtein[index])

    with open(outPath, "w") as file:
        for uniprotID in outlierUniprotIds:
            file.write(uniprotID + "\n") 


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="DTI Analysis for Sprint Embeddings")
    parser.add_argument( "command", choices=[
            "all",
            "generate_plot",
            "get_outlier_neighbors" ],
        help="Command to execute",
    )
    parser.add_argument("--drugSmiles",  help="Input file (.csv with id and SMILES columns)", required=True)
    parser.add_argument("--drugEmbeds",  help="Input file (SPRINT .npy file)", required=True)
    parser.add_argument("--proteome",  help="Input file (.csv with IDs and Gene Names)")
    parser.add_argument("--receptorEmbeds",  help="Input file (SPRINT .npy file)", required=True)
    parser.add_argument("--interestEmbeds",  help="Input file (SPRINT .npy file)")
    parser.add_argument("--interestLabel", default="Receptors of Interest",  help="Input file (SPRINT .npy file)")
    parser.add_argument("--xCutoff", default=5, help="PCA coordinate x cutoff")
    parser.add_argument("--yCutoff", default=5, help="PCA coordinate y cutoff")
    parser.add_argument("--plotPath",  help="Plot Output Path")
    parser.add_argument("--outlierDir",  help="Outlier Output Directory")
    args = parser.parse_args()

    drugEmbeddings = np.load(args.drugEmbeds)
    receptorEmbeddings = np.load(args.receptorEmbeds)

    interestEmbeddings = None
    if args.interestEmbeds != None:
        interestEmbeddings = np.load(args.interestEmbeds)

    #map drugs to embeddings
    drugsDF = pd.read_csv(args.drugSmiles, header=0)
    indexToID = dict(zip(range(len(drugEmbeddings)), drugsDF["id"].to_list()))
    idToSmiles = dict(zip(drugsDF["id"].to_list(), drugsDF["SMILES"].to_list()))
    idsToEmbeddings = dict(zip(drugsDF["id"].to_list(), drugEmbeddings))

    #Initialize PCA and UMAP, fit to drugs
    pca = PCA(n_components=2)
    reducer = umap.UMAP(random_state=0, metric='cosine')
    pca.fit(drugEmbeddings)
    reducer.fit(drugEmbeddings)

    #transform data
    pca_drugEmbeddings = pca.transform(drugEmbeddings)
    pca_receptorEmbeddings = pca.transform(receptorEmbeddings)
    umap_drugEmbeddings = reducer.transform(drugEmbeddings)
    umap_receptorEmbeddings = reducer.transform(receptorEmbeddings)

    #get PCA groups
    xCutoff = args.xCutoff
    yCutoff = args.yCutoff
    inliers = []
    xOutliers = []
    yOutliers = []
    xyOutliers = []
    idsToPCA = {}
    idsToUMAP = {}
    for i in range(len(pca_drugEmbeddings)):
        pcaEmbed = pca_drugEmbeddings[i]
        id = indexToID[i]
        idsToPCA[id] = pcaEmbed
        idsToUMAP[id] = umap_drugEmbeddings[i]
        if pcaEmbed[0] < xCutoff and pcaEmbed[1] < yCutoff:
            inliers.append(id)
        elif pcaEmbed[0] >= xCutoff and pcaEmbed[1] < yCutoff:
            xOutliers.append(id)
        elif pcaEmbed[0] < xCutoff and pcaEmbed[1] >= yCutoff:
            yOutliers.append(id)
        elif pcaEmbed[0] >= xCutoff and pcaEmbed[1] >= yCutoff:
            xyOutliers.append(id)

    if args.command == "generate_plot" or args.command == "all":
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
        
        plotA.set_xlabel("PCA1: Exp. Var = " + str(round(pca.explained_variance_ratio_[0], 4)))
        plotA.set_ylabel("PCA2: Exp. Var = " + str(round(pca.explained_variance_ratio_[1], 4)))
        plotA.set_title("PCA of Drug and Receptor DTI Embeddings")
        plotA.legend()


        #generate plotB: Umap of both Drug and Receptor Embeddings with PCA outliers labelled
        add_to_scatter(plotB, np.array(ids_to_coordinates(inliers, idsToUMAP)), 'b', "Inliers", True, 3)
        add_to_scatter(plotB, np.array(ids_to_coordinates(xOutliers, idsToUMAP)), 'c', "X Outliers", True, 3)
        add_to_scatter(plotB, np.array(ids_to_coordinates(yOutliers, idsToUMAP)), 'y', "Y Outliers", True, 3)
        if len(xyOutliers) > 0:
            add_to_scatter(plotB, np.array(ids_to_coordinates(xyOutliers, idsToUMAP)), "blueviolet", True, 3)

        add_to_scatter(plotB, umap_receptorEmbeddings, "pink", str(len(receptorEmbeddings)) + " Proteome Proteins", False, 1)

        if interestEmbeddings is not None:
            umap_interestEmbeddings = reducer.transform(interestEmbeddings)
            add_to_scatter(plotB, umap_interestEmbeddings, 'r', args.interestLabel, False, 2)

        plotB.set_xlabel("UMAP 1")
        plotB.set_ylabel("UMAP 2")
        plotB.set_title("UMAP of Drug and Receptor DTI Embeddings with PCA Outliers Labelled")
        plotB.legend()

        plt.savefig(args.plotPath)
    

    if args.command == "get_outlier_neighbors" or args.command == "all":
        os.makedirs(args.outlierDir, exist_ok=True)
        #construct protein embedding dictionary
        proteinDF = pd.read_csv(args.proteome, header=0)
        indexToProtein = dict(zip(range(len(receptorEmbeddings)), proteinDF["uniprot_id"]))

        #write drug outlier ids to file
        write_outlier_essentials(xOutliers, idToSmiles, args.outlierDir + "/xOutliers.csv")
        write_outlier_essentials(yOutliers, idToSmiles, args.outlierDir + "/yOutliers.csv")
        if len(xyOutliers) > 0:
            write_outlier_essentials(xyOutliers, idToSmiles, args.outlierDir + "/xyOutliers.csv")

        #get nearest neighbors of outliers
        nearest = NearestNeighbors(n_neighbors= 5, algorithm='auto', metric='cosine').fit(receptorEmbeddings)
        write_outlier_neighbors([idsToEmbeddings[id] for id in xOutliers], indexToProtein, args.outlierDir + "/xOutlierNeighbors.txt")
        write_outlier_neighbors([idsToEmbeddings[id] for id in yOutliers], indexToProtein, args.outlierDir + "/yOutlierNeighbors.txt")
        if len(xyOutliers) > 0:
            write_outlier_neighbors(xyOutliers, indexToProtein, args.outlierDir + "/xyOutlierNeighbors.txt")
    