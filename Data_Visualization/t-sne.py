import pandas as pd
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import numpy as np
import os
import pylab
from matplotlib.legend_handler import HandlerBase

def get_datasets():
    datasets = []
    # get all folders in a directory
    directory = 'C:/Users/Roman/Documents/Work/Depression_and_Immunology/Spring_Research/Data/'
    data_folders = [x[0] for x in os.walk(directory)][1:]

    print(data_folders)
    #print(data_folders)
    # go through each folder and get the files
    for data_folder in data_folders:
        subDir = data_folder + '/'
        data_files = [f for f in os.listdir(subDir)]
        # get gene counts
        print(data_files)
        gene_counts = subDir + [fileName for fileName in data_files if 'gene_counts' in fileName][0]
        # get meta data
        metadata = subDir + [fileName for fileName in data_files if 'metadata' in fileName][0]

        datasets.append((gene_counts, metadata, data_folder.split('/')[-1]))
    return datasets;

def get_geneCounts(gene_counts_path ,annot=None):
    # get the input files
    input_counts = gene_counts_path#"C:\\Users\\Roman\\Documents\\Work\\Depression_and_Immunology\\Data\\filteredData.csv"

    # get gene counts
    gene_counts = pd.read_csv(input_counts, index_col=0, header=0).T
    gene_counts.index = gene_counts.index.map(int)  # index is messed up for some reason

    if(not (annot == None)):
        # select only rows that appear in the annot file
        gene_counts = gene_counts.ix[annot.index.values]
        # we need to drop NaNs that were inserted from annot
    gene_counts = gene_counts.dropna()
    return gene_counts

def get_metadata(metadata_path):
    input_annot = metadata_path#"C:\\Users\\Roman\\Documents\\Work\\Depression_and_Immunology\\Scripts\\Data\\meta_data_with_paths_all_patients_final.csv"
    # get meta_data
    annot = pd.read_csv(input_annot).set_index('record_id')
    #annot = annot[annot['missing'] == False]

    return annot

def plot_tsne(df,annot, perplexity, plot_name = None):

    model = TSNE(n_components=2, perplexity=perplexity, n_iter = 5000)
    model = model.fit_transform(df)

    list_color = ['k', 'k', 'k', '#11c3be', '#89f5f2', '#d7fcfa', 'r']
    list_marker = ["o", "v", "^", "x", "x", "x", "x"]
    list_lab = ['Healthy', 'Depression', 'Bipolar', 'MDD B1', 'MDD B2', 'MDD B3', 'Control Dataset']

    tsne_df = pd.DataFrame()
    tsne_df['x'] = model[:, 0]
    tsne_df['y'] = model[:, 1]
    batchValues = annot['batch'].map({0.0: 'w', 1.0:'#11c3be', 2.0:'#89f5f2', 3.0:'#d7fcfa', 4.0:'r'})

    diagnosisValues = annot['scid_diagnosis'].map({0.0: '.', 1.0:'v', 2.0:'^', 3.0:'^', 4.0:'o'})

    class MarkerHandler(HandlerBase):
        def create_artists(self, legend, tup,xdescent, ydescent,
                            width, height, fontsize,trans):
            return [plt.Line2D([width/2], [height/2.],ls="",
                           marker=tup[1],color=tup[0], transform=trans)]

    for x,y,b,d in zip(tsne_df['x'], tsne_df['y'], batchValues, diagnosisValues):
        plt.scatter(x, y, c=b, marker=d)


    plt.legend(list(zip(list_color, list_marker)), list_lab, handler_map={tuple:MarkerHandler()})

    plt.title('T-SNE Perplexity: ' + str(perplexity) + ' ' + plot_name)
    plt.savefig('t-sne_perplex_' + str(perplexity) + '_' + plot_name + '.png')
    plt.clf()
    plt.cla()
    plt.close()

    #plt.show()

    #print(model)

#gene_counts = get_geneCounts(annot)
#print(gene_counts.shape)
## filter by mean
#gene_counts = gene_counts.loc[:, gene_counts.aggregate('mean') > 3]

#print(gene_counts)
#print(annot)

# filter by variance 1
#print(((gene_counts.var(axis=0)) > 1).sum())

#print(gene_counts.aggregate('mean'))
#plt.hist()
#print(gene_counts.var(axis=0).values)
#plt.hist(gene_counts.var(axis=0).values)
#gene_counts = gene_counts.loc[:, gene_counts.aggregate('') > 1]

#gene_counts = gene_counts.loc[:, gene_counts.aggregate('mean') > 3]
#print(gene_counts.shape)

#print(gene_counts.isnull().sum().sum())
datasets = get_datasets()

for dataset in datasets:
    for i in range(10,51,10):
        plot_tsne(get_geneCounts(dataset[0]), get_metadata(dataset[1]), i, dataset[2])
#plot_tsne(gene_counts, 20)
#plot_tsne(gene_counts, 30)


