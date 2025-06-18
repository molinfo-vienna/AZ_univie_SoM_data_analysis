import matplotlib.pyplot as plt

#plot_pca
def plot_pca(pcaDF,variance,databases,colors):
    plt.rcParams['font.family'] = 'STIXGeneral'
    # plt.rcParams['figure.figsize'] = 10.5, 7.5
    fig = plt.figure(dpi=150)
    axes = fig.add_subplot()
    axes.set_xlabel('PC1 ({:.2%})'.format(variance[0]), fontsize=18)
    axes.set_ylabel('PC2 ({:.2%})'.format(variance[1]), fontsize=18)

    for database, color in zip(databases, colors):
        indicesToKeep = pcaDF['db'] == database
        axes.scatter(pcaDF.loc[indicesToKeep, 'PC1'], pcaDF.loc[indicesToKeep, 'PC2'], c=color, s=3,alpha=0.8)
    axes.legend(databases, fontsize=16,loc='upper left')
    plt.tick_params(labelsize=16)
    #axes.set_xlim(left=-5.2, right=30.5)
    #axes.set_yticks([-10,0,10,20,30])
    #axes.set_ylim(bottom=-11.5, top=37)

    plt.tight_layout()
#     if len(set(pcaDF.db))==2:
#         plt.savefig('../../../20210901_analysis/pca.png')
    plt.show()


#plot_umap
def plot_umap_generic(umapDF, databases, colors, db_column='db'):
    """
    Generic UMAP plotting function that can handle different database column names
    
    Parameters:
    -----------
    umapDF : pandas DataFrame
        DataFrame containing UMAP coordinates and database labels
    databases : list
        List of database names to plot
    colors : list
        List of colors to use for plotting
    db_column : str, default='db'
        Name of the column containing database labels
    """
    plt.rcParams['font.family'] = 'STIXGeneral'
    fig = plt.figure(dpi=150)
    axes = fig.add_subplot()
    axes.set_xlabel('UMAP 1', fontsize=16)
    axes.set_ylabel('UMAP 2', fontsize=16)
    axes.set_aspect('equal', adjustable='box')
    
    colors = colors[:len(databases)]

    for database, color in zip(databases, colors):
        indicesToKeep = umapDF[db_column] == database
        axes.scatter(umapDF.loc[indicesToKeep, 'UMAP 1'], 
                    umapDF.loc[indicesToKeep, 'UMAP 2'], 
                    c=color, s=1, alpha=0.8)
    axes.legend(databases, fontsize=13, handletextpad=0.5, markerscale=2)
    plt.tick_params(labelsize=14)

    # Add small padding (e.g., 5%) to the limits
    x_padding = (umapDF['UMAP 1'].max() - umapDF['UMAP 1'].min()) * 0.05
    y_padding = (umapDF['UMAP 2'].max() - umapDF['UMAP 2'].min()) * 0.05
    plt.xlim(umapDF['UMAP 1'].min() - x_padding, umapDF['UMAP 1'].max() + x_padding)
    plt.ylim(umapDF['UMAP 2'].min() - y_padding, umapDF['UMAP 2'].max() + y_padding)

    plt.tight_layout()
    plt.show()


def plot_umap(umapDF, databases, colors):
    """Wrapper function maintaining original behavior"""
    plot_umap_generic(umapDF, databases, colors, db_column='db')


def plot_umap_subplots(umapDF, databases, colors):
    """Plot UMAP in current subplot using generic function"""
    # Call plot_umap_generic but skip figure creation
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.xlabel('UMAP 1', fontsize=16)
    plt.ylabel('UMAP 2', fontsize=16)
    plt.gca().set_aspect('equal', adjustable='box')
    
    colors = colors[:len(databases)]

    for database, color in zip(databases, colors):
        indicesToKeep = umapDF['db'] == database
        plt.scatter(umapDF.loc[indicesToKeep, 'UMAP 1'], 
                   umapDF.loc[indicesToKeep, 'UMAP 2'], 
                   c=color, s=1, alpha=0.8)
    plt.legend(databases, fontsize=13, handletextpad=0.5, markerscale=2)
    plt.tick_params(labelsize=14)

    # Add small padding (e.g., 5%) to the limits
    x_padding = (umapDF['UMAP 1'].max() - umapDF['UMAP 1'].min()) * 0.05
    y_padding = (umapDF['UMAP 2'].max() - umapDF['UMAP 2'].min()) * 0.05
    plt.xlim(umapDF['UMAP 1'].min() - x_padding, umapDF['UMAP 1'].max() + x_padding)
    plt.ylim(umapDF['UMAP 2'].min() - y_padding, umapDF['UMAP 2'].max() + y_padding)