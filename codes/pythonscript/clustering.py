import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib as mpl

# Set default figure size
mpl.rcParams['figure.figsize'] = (6, 6)

# Function to preprocess and filter data
def preprocess_data(file_path, sum_cols_no_atp, sum_cols_atp, ratio_cols):
    """
    Read and preprocess data from an Excel file.
    Filters rows based on the mean of specified columns and handles missing values.
    """
    # Read Excel file
    data = pd.read_excel(file_path)

    # Filter rows where the mean of columns exceeds 100
    data = data[(data[sum_cols_no_atp].mean(axis=1) > 100) | (data[sum_cols_atp].mean(axis=1) > 100)]

    # Handle missing values in ratio columns
    if data[ratio_cols].isnull().any().any():
        print("Warning: Missing values found in ratio columns. Filling with zeros.")
        data[ratio_cols] = data[ratio_cols].fillna(0)

    # Debugging: Print basic statistics
    print(data[ratio_cols].describe())
    
    return data

# Function to add noise if variability is low
def add_noise_if_low_variability(data, ratio_cols):
    ratio_data = data[ratio_cols].to_numpy()
    if np.var(ratio_data) < 1e-5:
        print("Adding slight noise to data to improve clustering.")
        ratio_data += np.random.normal(0, 1e-5, ratio_data.shape)
    return ratio_data

# Function to create AnnData object
def create_anndata(ratio_data, data, ratio_cols):
    adata = sc.AnnData(ratio_data)
    adata.obs['sample'] = data.index.astype(str)
    adata.obs['genes'] = data['genes'].tolist()
    adata.var['features'] = ratio_cols
    return adata

# Function to perform clustering and UMAP
def perform_clustering(adata, resolution=0.75):
    sc.tl.pca(adata, n_comps=min(adata.shape[0], adata.shape[1]) - 1)
    sc.pp.neighbors(adata, n_neighbors=10, use_rep='X_pca')
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", random_state=42)
    return adata

# Function to save UMAP plot with annotations
def save_umap_plot(adata, output_file1,output_file2):
    fig = sc.pl.umap(adata, color='leiden', title='UMAP with Leiden Clustering',
                     legend_loc='on data', return_fig=True)
    umap_coords = adata.obsm['X_umap']
    labels = adata.obs['genes'].tolist()
    ax = fig.axes[0]

    # for i in range(len(umap_coords)):
    #     x, y = umap_coords[i, 0], umap_coords[i, 1]
    #     label = labels[i]
    #     ax.text(x, y, str(label), fontsize=3, alpha=0.7)

    plt.xlabel('UMAP1', fontsize=12)
    plt.ylabel('UMAP2', fontsize=12)
    plt.savefig(output_file1, format='pdf')
    plt.close()
    print(f"UMAP plot saved to: {output_file1}")

    fig = sc.pl.umap(adata, color='leiden', title='UMAP with Leiden Clustering',
                     legend_loc='on data', return_fig=True)
    umap_coords = adata.obsm['X_umap']
    labels = adata.obs['genes'].tolist()
    ax = fig.axes[0]

    for i in range(len(umap_coords)):
        x, y = umap_coords[i, 0], umap_coords[i, 1]
        label = labels[i]
        ax.text(x, y, str(label), fontsize=3, alpha=0.7)

    plt.xlabel('UMAP1', fontsize=12)
    plt.ylabel('UMAP2', fontsize=12)
    plt.savefig(output_file2, format='pdf')
    plt.close()
    print(f"UMAP plot saved to: {output_file2}")


# Function to plot time series for each cluster
def plot_clusters_time_series(data, ratio_cols, cluster_col, output_file):
    unique_clusters = data[cluster_col].unique()
    clusters_per_page = 6  # For a 2x3 grid

    with PdfPages(output_file) as pdf:
        num_pages = int(np.ceil(len(unique_clusters) / clusters_per_page))
        for page in range(num_pages):
            fig, axes = plt.subplots(3, 2, figsize=(10, 8))
            axes = axes.flatten()

            clusters_on_page = unique_clusters[page * clusters_per_page:(page + 1) * clusters_per_page]
            for idx, cluster in enumerate(clusters_on_page):
                ax = axes[idx]
                cluster_data = data[data[cluster_col] == cluster][ratio_cols]
                if cluster_data.shape[0] == 0:
                    ax.set_visible(False)
                    continue

                cluster_mean = cluster_data.mean(axis=0)
                cluster_std = cluster_data.std(axis=0) * 0.5

                ax.plot(ratio_cols, cluster_mean, color='blue', linewidth=2)
                ax.fill_between(ratio_cols, cluster_mean - cluster_std, cluster_mean + cluster_std,
                                color='blue', alpha=0.2)
                ax.set_xticks([])

            for idx in range(len(clusters_on_page), len(axes)):
                axes[idx].set_visible(False)

            plt.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Time series plots saved to: {output_file}")

# Function to save ranked t180 plot
def plot_t180_ranking(data, output_file):
    ranked_data = data.sort_values(by='ratio_t180_x', ascending=False).reset_index(drop=True)
    ranked_data['Rank'] = ranked_data.index + 1

    def assign_color(order):
        if order > 1:
            return '#9e0142'
        elif order < -1:
            return '#5e4fa2'
        else:
            return 'lightgray'

    ranked_data['Color'] = ranked_data['ratio_t180_x'].apply(assign_color)
    plt.figure(figsize=(12, 8))
    plt.scatter(ranked_data['ratio_t180_x'], ranked_data['Rank'], c=ranked_data['Color'], s=20)

    for _, row in ranked_data.iterrows():
        if abs(row['ratio_t180_x']) > 1:
            plt.text(row['ratio_t180_x'], row['Rank'], row['genes'], fontsize=8, ha='center', rotation=90)

    plt.title("Entities Ranked by Ratio_t180", fontsize=16)
    plt.xlabel("Ratio_t180", fontsize=14)
    plt.ylabel("Rank", fontsize=14)
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    print(f"t180 ranking plot saved to: {output_file}")

# Main function
def main():
    file_path = 'filtered_table_with_reactions_m9_median_no_duplicates.xlsx'
    umap_output_file = "umap_plot_m9.pdf"
    umap_output_file2 = "umap_plot_m9_text.pdf"
    time_series_output_file = "time_series_clusters_m9.pdf"
    t180_ranking_output_file = "t180_ranking_plot_m9.pdf"
    output_excel_file = 'clustered_data_with_umap_m9.xlsx'

    sum_cols_no_atp = ['avg_t0', 'avg_t5', 'avg_t30', 'avg_t60', 'avg_t120', 'avg_t180']
    sum_cols_atp = ['avg_t0_ATP', 'avg_t5_ATP', 'avg_t30_ATP', 'avg_t60_ATP', 'avg_t120_ATP', 'avg_t180_ATP']
    ratio_cols = ['ratio_t0_x', 'ratio_t5_x', 'ratio_t30_x', 'ratio_t60_x', 'ratio_t120_x', 'ratio_t180_x']

    data = preprocess_data(file_path, sum_cols_no_atp, sum_cols_atp, ratio_cols)
    ratio_data = add_noise_if_low_variability(data, ratio_cols)
    adata = create_anndata(ratio_data, data, ratio_cols)
    adata = perform_clustering(adata)

    data['Leiden_Cluster'] = adata.obs['leiden'].tolist()

    save_umap_plot(adata, umap_output_file,umap_output_file2)
    plot_clusters_time_series(data, ratio_cols, 'Leiden_Cluster', time_series_output_file)
    # plot_t180_ranking(data, t180_ranking_output_file)

    data.to_excel(output_excel_file, index=False)
    print(f"Clustered data with UMAP and Leiden saved to: {output_excel_file}")

    ### for lb 

    file_path = 'filtered_table_with_reactions_lb_median_no_duplicates.xlsx'
    umap_output_file = "umap_plot_lb.pdf"
    umap_output_file2 = "umap_plot_lb_text.pdf"
    time_series_output_file = "time_series_clusters_lb.pdf"
    t180_ranking_output_file = "t180_ranking_plot_lb.pdf"
    output_excel_file = 'clustered_data_with_umap_lb.xlsx'

    sum_cols_no_atp = ['avg_t0', 'avg_t5', 'avg_t30', 'avg_t60', 'avg_t120', 'avg_t180']
    sum_cols_atp = ['avg_t0_ATP', 'avg_t5_ATP', 'avg_t30_ATP', 'avg_t60_ATP', 'avg_t120_ATP', 'avg_t180_ATP']
    ratio_cols = ['ratio_t0_x', 'ratio_t5_x', 'ratio_t30_x', 'ratio_t60_x', 'ratio_t120_x', 'ratio_t180_x']

    data = preprocess_data(file_path, sum_cols_no_atp, sum_cols_atp, ratio_cols)
    ratio_data = add_noise_if_low_variability(data, ratio_cols)
    adata = create_anndata(ratio_data, data, ratio_cols)
    adata = perform_clustering(adata)

    data['Leiden_Cluster'] = adata.obs['leiden'].tolist()

    save_umap_plot(adata, umap_output_file,umap_output_file2)
    plot_clusters_time_series(data, ratio_cols, 'Leiden_Cluster', time_series_output_file)
    # plot_t180_ranking(data, t180_ranking_output_file)

    data.to_excel(output_excel_file, index=False)
    print(f"Clustered data with UMAP and Leiden saved to: {output_excel_file}")



if __name__ == "__main__":
    main()
