import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from matplotlib import pyplot as plt
from sklearn.mixture import GaussianMixture as GMM
from scipy.stats import norm

def getOverlapGenes(genes1, genes2, homology, species1, species2):
    gene2_121 = np.array(homology.loc[:,species2])
    for i, gene2 in enumerate(genes2):
        gene1 = " "
        if gene2 in gene2_121:
            gene1 = homology.loc[np.where(gene2_121 == gene2)[0][0],species1]
        genes2[i] = f"{gene1}/{gene2}"

    gene1_121 = np.array(homology.loc[:,species1])
    for i, gene1 in enumerate(genes1):
        gene2 = " "
        if gene1 in gene1_121:
            gene2 = homology.loc[np.where(gene1_121 == gene1)[0][0],species2]
        genes1[i] = f"{gene1}/{gene2}"
    
    return(genes1, genes2)

def getOverlapGenesMulti(homology, geneSpecDict):
    for data in geneSpecDict:
        genes = data["genes"]
        species = data["species"]
        gene_121 = np.array(homology.loc[:,species])
        for i, gene in enumerate(genes):
            allGenes = gene
            if gene in gene_121:
                allGenes = "/".join(homology.loc[np.where(gene_121==gene)[0][0],:])
            genes[i] = allGenes
    
    return(geneSpecDict)

def plotScores(scoreMat, figCols = 4, cutoff = 1.5):
    fig, axs = plt.subplots((len(scoreMat.columns)+1)//figCols, figCols, figsize=(10, 10))
    #sigs = []
    for i,scoreCol in enumerate(scoreMat.columns):
        scoreData = scoreMat[scoreCol]
        x, y = i//figCols, i%figCols
        axs[x, y].hist(((scoreData-np.mean(scoreData))/np.std(scoreData)), bins=100)
        axs[x, y].axvline(cutoff, color='k', linestyle='dashed', linewidth=1)
        axs[x, y].set_title(scoreCol)
        #sigs.append(((scoreData-np.mean(scoreData))/np.std(scoreData)) > cutoff)

    fig.show()
    #return(sigs)
    
def scoreGeneSig(adata, geneSig, translate = False, toGenes=None, fromGenes=None):
    for j,sigName in enumerate(geneSig.columns):
        clustGenes = geneSig.iloc[:,j].dropna()
        if translate:
            for i,mGene in enumerate(clustGenes):
                indexOver = np.where(fromGenes == mGene)[0]
                if indexOver.size > 0:
                    clustGenes[i] = toGenes[indexOver[0]]
        sc.tl.score_genes(adata, clustGenes, score_name=f"{sigName}Score")

def getNewLabels(adata, ogLabels, scoreNames,labelDict):
#   ogLabelScore = pd.DataFrame(np.zeros((len(ogLabels),len(scoreNames))),index=ogLabels, columns=scoreNames)
    ogLabelScoreMe = pd.DataFrame(np.zeros((len(ogLabels),len(scoreNames))),index=ogLabels, columns=scoreNames)
    newBClabel = list(ogLabels.copy())
    
    for score in ogLabelScoreMe.columns:
        scorMe = np.median(adata.obs[score])
        #print(f"{score}: {scorMe}\n")
        for i,leid in enumerate(ogLabelScoreMe.index):
            adataCat = adata[adata.obs.leiden==leid]
            #ogLabelScore.loc[leid,score] = scorMe
            #ogLabelScoreMe.loc[leid,score] = np.round(sum(adataCat.obs[score] > scorMe)/len(adataCat.obs[score]),decimals=4)
            ogLabelScoreMe.loc[leid,score] = np.round(sum(adataCat.obs[score] > scorMe)/len(adataCat.obs[score])*scorMe,decimals=4)

    #print(ogLabelScoreMe)
    
    for i,leid in enumerate(ogLabelScoreMe.index):
        if(np.max(ogLabelScoreMe.loc[leid,:])): #  > 0.5 and scorMe > 0
            newBClabel[i] = labelDict[ogLabelScoreMe.columns[np.argmax(ogLabelScoreMe.loc[leid,:])]]
        else:
            newBClabel[i] = "inter"
        
    adata.obs["cellState"] = [newBClabel[int(lei)] for lei in adata.obs.leiden]
    return(newBClabel, ogLabelScoreMe)


def scoreAndLabel(adata, sigGenes, labelOfSigGenes, ogLabel="leiden",translate = False, toGenes=None, fromGenes=None):
    scoreGeneSig(adata, sigGenes, translate = translate, toGenes=toGenes, fromGenes=fromGenes)
    ogLabels = adata.obs[ogLabel].cat.categories
    scoreNames = [f"{sigName}Score" for sigName in sigGenes.columns]
    labelDict = dict(zip(scoreNames,labelOfSigGenes))
    newBClabel, ogLabelScoreMe = getNewLabels(adata, ogLabels, scoreNames, labelDict)
    return(scoreNames, newBClabel, ogLabelScoreMe)

def addIndvLabel(adata, scoreNames, obsLabel="zsig", cutoff=1):
    scoreMat = adata.obs[scoreNames]
    if(obsLabel=="topPer"):
        out = topPercent(scoreMat, cutoff = cutoff)
    else:#(obsLabel=="zsig")
        out = zScores(scoreMat, cutoff = cutoff)
    adata.obs[obsLabel] = out
    return(adata)

def zScores(scoreMat, cutoff = 1.5):
    sigs = []
    for i,scoreCol in enumerate(scoreMat.columns):
        scoreData = scoreMat[scoreCol]
        zscore = ((scoreData-np.mean(scoreData))/np.std(scoreData))
        sigs.append(zscore)
    sigScore = pd.DataFrame(np.array(sigs)).T
    simple = []
    for i,cell in enumerate(sigScore.index):
        names = scoreMat.columns
        sigNames = names[np.array((sigScore.loc[cell]> cutoff).values)]
        sigName = "out"
        if len(sigNames) > 0:
            sigName = str(names[np.argmax(sigScore.loc[cell])])[:-5]
        simple.append(sigName)
    return(np.array(simple))

def topPercent(scoreMat, cutoff = .80):
    perc = pd.DataFrame(np.zeros(scoreMat.shape),columns=scoreMat.columns,index=scoreMat.index)
    for i,scoreCol in enumerate(scoreMat.columns):
        scoreData = scoreMat[scoreCol]
        for j, cell in enumerate(scoreData.index):
            perc.loc[cell,scoreCol] = sum(scoreData<scoreData[cell])/len(scoreData)
    simple=[]
    for i,cell in enumerate(perc.index):
        names = perc.columns
        sigNames = names[np.array((perc.loc[cell]> cutoff).values)]
        sigName = "out"
        if len(sigNames) > 0:
            sigName = str(names[np.argmax(perc.loc[cell])])[:-5]
        simple.append(sigName)
    return(np.array(simple))


from sklearn.mixture import BayesianGaussianMixture as GMM
#from sklearn.mixture import GaussianMixture as GMM

def gmmScoreGeneSig(scoreMat, plotLen = 3, show=False):
    scoreNames = scoreMat.columns
    numScores = len(scoreNames)
    if(show):
        fig, axs = plt.subplots((numScores//plotLen)+1,plotLen)
        plt.rcParams["figure.figsize"] = (15,5)

    dfScoreBoundry = pd.DataFrame(np.zeros(numScores),scoreNames, columns=["boundry"])
    gmm = GMM(n_components = 2, random_state=10)#, init_params="random_from_data")#, means_init=meansInit)
    #binEx = np.arange(0.5,10,10/200).reshape(-1,1)

    for i, scoreName in enumerate(scoreNames):
        scoreCount = np.array(scoreMat[scoreName]).reshape(-1, 1)
        fitGMM = gmm.fit(scoreCount)
        mean = fitGMM.means_  
        covs  = fitGMM.covariances_
        weights = fitGMM.weights_
        #binEx = np.arange(min(min(mean),max(mean)),max(scoreCount),0.01).reshape(-1,1)
        #binEx = np.arange(min(min(mean),max(mean))[0],max(scoreCount)[0],0.01).reshape(-1,1)
        #binEx = np.arange(min(mean),np.percentile(scoreCount,95),0.005).reshape(-1,1)
        binEx = np.arange(np.percentile(scoreCount,10),np.percentile(scoreCount,95),0.01).reshape(-1,1)
        
        #print(f"{min(mean)} {np.percentile(scoreCount,85)}")
        fitGmmBound = fitGMM.predict(binEx)
        furtherBound = fitGmmBound[-1]
        fitGmmBoundUniq = np.unique(fitGmmBound)

        #print(f"bound {fitGmmBound}")
        #print(furtherBound)
        #print(fitGmmBound)
        if (len(fitGmmBoundUniq) == 2):
            if(fitGmmBound[0] == fitGmmBound[-1]):
                furtherBound = fitGmmBoundUniq[fitGmmBoundUniq != furtherBound][0]
            scoreBoundry = binEx[np.where(fitGmmBound == furtherBound)[0][0]][0]
        else:
            scoreBoundry = max(scoreCount)
        #naiveBoundry = np.log10(int(scoreDisc.loc["90%",scoreName])+1)

        dfScoreBoundry.loc[scoreName] = scoreBoundry
        
        if(show):
            x_axis = np.arange(min(scoreCount), max(scoreCount), 0.01)
            y_axis0 = norm.pdf(x_axis, float(mean[0][0]), np.sqrt(float(covs[0][0][0])))*weights[0] # 1st gaussian
            y_axis1 = norm.pdf(x_axis, float(mean[1][0]), np.sqrt(float(covs[1][0][0])))*weights[1] # 2nd gaussian

            # Plot 2
            x,y = i//plotLen, i%plotLen
            axs[x,y].set_title(scoreName)
            #axs[x,y].axvline(naiveBoundry, c='C3', linestyle='dashed', linewidth=1) #red
            axs[x,y].axvline(scoreBoundry, c='C2', linestyle='dashed', linewidth=1)  #green
            axs[x,y].hist(scoreCount, density=True, color='black', bins=100)        
            axs[x,y].plot(x_axis, y_axis0, lw=3, c='C6')                            #pink
            axs[x,y].plot(x_axis, y_axis1, lw=3, c='C1')                            #orange
            axs[x,y].plot(x_axis, y_axis0+y_axis1, lw=3, c='C0', ls=':')            #dotted blue

    if(show):
        plt.tight_layout(pad=1.0)
        plt.show()
        #sc.set_figure_params(scanpy=True, dpi=100, dpi_save=150, fontsize=10, format='png')
    
    scoreIDs = scoreMat.copy()
    scoreZscores = scoreMat.apply(stats.zscore)
    scoreID = np.array(scoreNames)
    #pdb.set_trace()

    for scoreName in scoreNames:
        #print(scoreName)
        #print(dfScoreBoundry.loc[scoreName].values[0])
        #if(dfScoreBoundry.loc[scoreName].values[0] > 0):
        scoreIDs.loc[:,scoreName] = (scoreMat.loc[:,scoreName] > dfScoreBoundry.loc[scoreName].values[0]).astype("float")

    classification = np.empty(len(scoreMat), dtype="object")
    i = 0
    for cellBar, scoreBool in scoreIDs.iterrows():
        scoreBool = scoreBool.astype("bool")
        numscorees = sum(scoreBool)
        if (numscorees == 1):
            classif = scoreID[scoreBool.values][0][:-5]#.values
        elif (numscorees > 1):
            #pdb.set_trace()
            #classif = "Doublet"
            maxTrue = np.argmax(scoreZscores.loc[cellBar,scoreID[scoreBool]])   
            #pdb.set_trace()
            classif = scoreID[scoreBool][maxTrue][:-5]#.values
        else:
            classif = "Negative"
        classification[i] = classif
        i = i + 1
        
    return(classification)

from sklearn.mixture import BayesianGaussianMixture as GMM
#from sklearn.mixture import GaussianMixture as GMM

def gmmScoreGeneSig(scoreMat, plotLen = 3, show=False):
    scoreNames = scoreMat.columns
    numScores = len(scoreNames)
    if(show):
        fig, axs = plt.subplots((numScores//plotLen)+1,plotLen)
        plt.rcParams["figure.figsize"] = (15,5)

    dfScoreBoundry = pd.DataFrame(np.zeros(numScores),scoreNames, columns=["boundry"])
    gmm = GMM(n_components = 2, random_state=10)#, init_params="random_from_data")#, means_init=meansInit)
    #binEx = np.arange(0.5,10,10/200).reshape(-1,1)

    for i, scoreName in enumerate(scoreNames):
        scoreCount = np.array(scoreMat[scoreName]).reshape(-1, 1)
        fitGMM = gmm.fit(scoreCount)
        mean = fitGMM.means_  
        covs  = fitGMM.covariances_
        weights = fitGMM.weights_
        #binEx = np.arange(min(min(mean),max(mean)),max(scoreCount),0.01).reshape(-1,1)
        #binEx = np.arange(min(min(mean),max(mean))[0],max(scoreCount)[0],0.01).reshape(-1,1)
        binEx = np.arange(np.percentile(scoreCount,10),np.percentile(scoreCount,95),0.01).reshape(-1,1)
        
        #print(f"{min(mean)} {np.percentile(scoreCount,85)}")
        fitGmmBound = fitGMM.predict(binEx)
        furtherBound = fitGmmBound[-1]
        fitGmmBoundUniq = np.unique(fitGmmBound)

        #print(f"bound {fitGmmBound}")
        #print(furtherBound)
        #print(fitGmmBound)
        if (len(fitGmmBoundUniq) == 2):
            if(fitGmmBound[0] == fitGmmBound[-1]):
                furtherBound = fitGmmBoundUniq[fitGmmBoundUniq != furtherBound][0]
            scoreBoundry = binEx[np.where(fitGmmBound == furtherBound)[0][0]][0]
        else:
            scoreBoundry = max(scoreCount)
        #naiveBoundry = np.log10(int(scoreDisc.loc["90%",scoreName])+1)

        dfScoreBoundry.loc[scoreName] = scoreBoundry
        
        if(show):
            x_axis = np.arange(min(scoreCount), max(scoreCount), 0.01)
            y_axis0 = norm.pdf(x_axis, float(mean[0][0]), np.sqrt(float(covs[0][0][0])))*weights[0] # 1st gaussian
            y_axis1 = norm.pdf(x_axis, float(mean[1][0]), np.sqrt(float(covs[1][0][0])))*weights[1] # 2nd gaussian

            # Plot 2
            x,y = i//plotLen, i%plotLen
            axs[x,y].set_title(scoreName)
            #axs[x,y].axvline(naiveBoundry, c='C3', linestyle='dashed', linewidth=1) #red
            axs[x,y].axvline(scoreBoundry, c='C2', linestyle='dashed', linewidth=1)  #green
            axs[x,y].hist(scoreCount, density=True, color='black', bins=100)        
            axs[x,y].plot(x_axis, y_axis0, lw=3, c='C6')                            #pink
            axs[x,y].plot(x_axis, y_axis1, lw=3, c='C1')                            #orange
            axs[x,y].plot(x_axis, y_axis0+y_axis1, lw=3, c='C0', ls=':')            #dotted blue

    if(show):
        plt.tight_layout(pad=1.0)
        plt.show()
        #sc.set_figure_params(scanpy=True, dpi=100, dpi_save=150, fontsize=10, format='png')
    
    scoreIDs = scoreMat.copy()
    scoreID = np.array(scoreNames)
    #pdb.set_trace()

    for scoreName in scoreNames:
        #print(scoreName)
        #print(dfScoreBoundry.loc[scoreName].values[0])
        #if(dfScoreBoundry.loc[scoreName].values[0] > 0):
        scoreIDs.loc[:,scoreName] = (scoreMat.loc[:,scoreName] > dfScoreBoundry.loc[scoreName].values[0]).astype("float")

        
    classification = np.empty(len(scoreMat), dtype="object")
    i = 0
    for cellBar, scoreBool in scoreIDs.iterrows():
        scoreBool = scoreBool.astype("bool")
        numscorees = sum(scoreBool)
        if (numscorees == 1):
            classif = scoreID[scoreBool.values][0][:-5]#.values
        elif (numscorees > 1):
            #pdb.set_trace()
            #classif = "Doublet"
            maxTrue = np.argmax(scoreMat.loc[cellBar,scoreID[scoreBool]])   
            #pdb.set_trace()
            classif = scoreID[scoreBool][maxTrue][:-5]#.values
        else:
            classif = "Negative"
        classification[i] = classif
        i = i + 1
        
    return(classification)

def gmmScoreGeneSigOLD(scoreMat, meansInit=[[0],[0.5]],plotLen = 3, show=False):
    scoreNames = scoreMat.columns
    numScores = len(scoreNames)
    if(show):
        fig, axs = plt.subplots((numScores//plotLen)+1,plotLen)
        plt.rcParams["figure.figsize"] = (15,5)

    dfScoreBoundry = pd.DataFrame(np.zeros(numScores),scoreNames, columns=["boundry"])
    gmm = GMM(n_components = 2, random_state=10, covariance_type = 'full', n_init=2, means_init=meansInit)
    #binEx = np.arange(0.5,10,10/200).reshape(-1,1)

    for i, scoreName in enumerate(scoreNames):
        scoreCount = np.array(scoreMat[scoreName]).reshape(-1, 1)
        fitGMM = gmm.fit(scoreCount)
        mean = fitGMM.means_  
        covs  = fitGMM.covariances_
        weights = fitGMM.weights_
        #print(mean)
        binEx = np.arange(min(mean),max(mean),0.01).reshape(-1,1)
        fitGmmBound = fitGMM.predict(binEx)
        #pdb.set_trace()
        #print(fitGmmBound)
        try:
            scoreBoundry = binEx[np.where(fitGmmBound == 1)[0][0]][0]
        except:
            scoreBoundry = max(scoreCount)
        #naiveBoundry = np.log10(int(scoreDisc.loc["90%",scoreName])+1)

        dfScoreBoundry.loc[scoreName] = scoreBoundry
        
        if(show):
            x_axis = np.arange(min(scoreCount), max(scoreCount), 0.05)
            y_axis0 = norm.pdf(x_axis, float(mean[0][0]), np.sqrt(float(covs[0][0][0])))*weights[0] # 1st gaussian
            y_axis1 = norm.pdf(x_axis, float(mean[1][0]), np.sqrt(float(covs[1][0][0])))*weights[1] # 2nd gaussian

            # Plot 2
            x,y = i//plotLen, i%plotLen
            axs[x,y].set_title(scoreName)
            #axs[x,y].axvline(naiveBoundry, c='C3', linestyle='dashed', linewidth=1) #red
            axs[x,y].axvline(scoreBoundry, c='C2', linestyle='dashed', linewidth=1)  #green
            axs[x,y].hist(scoreCount, density=True, color='black', bins=100)        
            axs[x,y].plot(x_axis, y_axis0, lw=3, c='C6')                            #pink
            axs[x,y].plot(x_axis, y_axis1, lw=3, c='C1')                            #orange
            axs[x,y].plot(x_axis, y_axis0+y_axis1, lw=3, c='C0', ls=':')            #dotted blue

    if(show):
        plt.tight_layout(pad=1.0)
        plt.show()
        sc.set_figure_params(scanpy=True, dpi=100, dpi_save=150, fontsize=10, format='png')
    
    scoreIDs = scoreMat.copy()
    scoreID = np.array(scoreNames)
    for scoreName in scoreNames:
        #print(scoreName)
        #print(dfScoreBoundry.loc[scoreName].values[0])
        scoreIDs.loc[:,scoreName] = scoreMat.loc[:,scoreName] > dfScoreBoundry.loc[scoreName].values[0]
        
    classification = np.empty(len(scoreMat), dtype="object")
    i = 0
    for cellBar, scoreBool in scoreIDs.iterrows():
        numscorees = sum(scoreBool)
        if (numscorees == 1):
            classif = scoreID[scoreBool.values][0][:-5]#.values
        elif (numscorees > 1):
            classif = "Doublet"
        else:
            classif = "Negative"
        classification[i] = classif
        i = i + 1
        
    return(classification)


from sklearn.metrics.pairwise import cosine_similarity
from scipy import stats

lenPCs = 20

label1 = "mrtxScore"#'scBasalScore'#"S_score"
label2 = "vehScore"#'scClassicalScore'#"G2M_score"

#find pcs pearson correlated with given score
def findDiffPCs(adata, label1, label2=None, lenPCs=20, show=True):
    #intitalize to 0
    #cosSimDif = np.zeros(lenPCs)
    perCorrDif = np.zeros(lenPCs)

    #allCosB = np.zeros(lenPCs)
    #allCosC = np.zeros(lenPCs)

    allCorrB = np.zeros(lenPCs)
    allCorrC = np.zeros(lenPCs)
    
    #if pca hasnt been calcualted calulcate it 
    if("X_pca" not in adata.obsm): 
        sc.tl.pca(adata)
    
    #if given a categorical variable split it into two numbered scores
    if(label2==None):    #if(adata.obs[label1].dtype != float or adata.obs[label1].dtype != int or adata.obs[label2].dtype != float or adata.obs[label2].dtype != int):
        cat1 = adata.obs[label1].cat.categories[0]
        cat2 = adata.obs[label1].cat.categories[1]

        adata.obs[f"{cat1}Score"] = [1 if m==cat1 else 0 for m in adata.obs[label1]]
        adata.obs[f"{cat2}Score"] = [1 if m==cat2 else 0 for m in adata.obs[label1]]
        
        label1 = f"{cat1}Score"
        label2 = f"{cat2}Score"
    
    for i in range(lenPCs):
        #print(f"PC {i+1}")
        #cosB = cosine_similarity(adata.obs[f"{label1}"].values.reshape(1, -1), adata.obsm["X_pca"][:,i].reshape(1, -1))[0][0]
        #cosC = cosine_similarity(adata.obs[f"{label2}"].values.reshape(1, -1), adata.obsm["X_pca"][:,i].reshape(1, -1))[0][0]
        corrB, pval = stats.pearsonr(adata.obs[f"{label1}"].values, adata.obsm["X_pca"][:,i])
        corrC, pval = stats.pearsonr(adata.obs[f"{label2}"].values, adata.obsm["X_pca"][:,i])

        #allCosB[i] = cosB
        #allCosC[i] = cosC

        allCorrB[i] = corrB
        allCorrC[i] = corrC

        #cosSimDif[i] = np.abs(cosB - cosC)
        perCorrDif[i] = np.abs(corrB - corrC)

    #made with the help of chatGPT

    numPCs = len(allCorrB)
    differences = np.abs(allCorrB - allCorrC)
    indSort = np.argsort(differences)

    if(show):
        # Plotting
        plt.figure(figsize=(8, 6))
        bar_width = 0.35

        index = np.arange(1, numPCs+1)
        bars1 = plt.bar(index - bar_width/2, allCorrB, bar_width, color='b', label=label1)

        # Plot Grade 2
        bars2 =plt.bar(index + bar_width/2, allCorrC, bar_width, color='r', label=label2)

        for i in range(len(index)):
            plt.text((index[i] - bar_width/2)+0.2, max(0,max(allCorrB[i], allCorrC[i])) + 0.01, f'{differences[i]:.2f}', ha='center')

        # Highlight the student with the largest difference
        #plt.plot([max_diff_index + 1 - bar_width/2, max_diff_index + 1 + bar_width/2],[allCorrB[max_diff_index], allCorrC[max_diff_index]], 'g--')  

        #bars1[max_diff_index].set_edgecolor('yellow')
        #bars1[max_diff_index].set_linewidth(2)
        #bars2[max_diff_index].set_edgecolor('yellow')
        #bars2[max_diff_index].set_linewidth(2)

        #plt.plot(np.arange(1, numPCs+1), allCorrB, 'bo', label='Basal', markersize=10)  # Plot Grade 1
        #plt.plot(np.arange(1, numPCs+1), allCorrC, 'ro', label='Classical', markersize=10)  # Plot Grade 2
        #plt.plot([max_diff_index + 1] * 2, [allCorrB[max_diff_index], allCorrC[max_diff_index]], 'g--')  # Connect the grades for the student with the largest difference

        # Customizing plot
        plt.title(f'Correlation of PCs with {label1} and {label2}')
        plt.xlabel('PCs')
        plt.ylabel('Pearson Correlation')
        plt.xticks(np.arange(1, numPCs+1))  # Set ticks for each student
        plt.legend()

        plt.grid(True)
        plt.tight_layout()
        plt.show()
    
    return indSort[::-1]

def mannwhitScore(adata, obsLabel, cellstate1, cellState2, score, alt="greater"):

    return stats.mannwhitneyu(adata[adata.obs[obsLabel]==cellstate1].obs[score],
                              adata[adata.obs[obsLabel]==cellstate1].obs[score], alternative=alt)