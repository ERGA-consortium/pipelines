import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def parse_arguments():
    parser = argparse.ArgumentParser(description='Create protein length dot plot from BLASTP results')
    parser.add_argument('-i', '--input', required=True, help='Input file with BLASTP results')
    parser.add_argument('-o', '--output', default='Protein_length', help='Output base name')
    parser.add_argument('-r', '--ratio', type=float, default=1.5,
                        help='Length ratio threshold for highlighting (default: 1.5)')
    return parser.parse_args()

def load_data(file_path):
    columns = ['qseqid', 'sseqid', 'qlen', 'slen', 'length', 'pident', 'evalue', 'bitscore']
    df = pd.read_csv(file_path, sep=r'\s+', header=None, names=columns)
    return df[['qlen', 'slen']].astype(int)

def calculate_density(x, y, bins=50):
    """Calculate 2D density using histogram"""
    hist, xedges, yedges = np.histogram2d(x, y, bins=bins)
    x_idx = np.clip(np.digitize(x, xedges), 0, hist.shape[0]-1)
    y_idx = np.clip(np.digitize(y, yedges), 0, hist.shape[1]-1)
    return hist[x_idx, y_idx]

def create_plot(df, output_base, ratio_threshold):
    plt.figure(figsize=(10, 10))
    sns.set_theme(style="whitegrid")
    
    df['density'] = calculate_density(df['qlen'], df['slen'])
    df['ratio'] = df['qlen'] / df['slen']
    df['outlier'] = (df['ratio'] > ratio_threshold) | (df['ratio'] < 1/ratio_threshold)
    
    corr_matrix = np.corrcoef(df['qlen'], df['slen'])
    r_squared = corr_matrix[0,1] ** 2
    
    ax = sns.scatterplot(
        data=df,
        x='qlen',
        y='slen',
        hue='density',
        style='outlier',
        markers={False: 'o', True: '^'},
        palette="viridis",
        size=0.7,
        edgecolor='none',
        legend='auto'
    )
    
    min_len = df[['qlen', 'slen']].min().min()
    max_len = df[['qlen', 'slen']].max().max()
    ax.plot([min_len, max_len], [min_len, max_len], '--', color='gray', lw=1, label='Equal length')
    
    ax.text(0.95, 0.95, f'R² = {r_squared:.2f}', transform=ax.transAxes,
            fontsize=12, ha='right', va='top', 
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
    legend_elements = [
        plt.Line2D([], [], marker='o', color='w', markerfacecolor='gray', 
                  markersize=8, label='Normal proteins'),
        plt.Line2D([], [], marker='^', color='w', markerfacecolor='gray', 
                  markersize=8, label=f'Length outliers (ratio > {ratio_threshold}x)'),
        plt.Line2D([], [], color='gray', linestyle='--', lw=1, label='Equal length')
    ]
    
    ax.legend(handles=legend_elements, loc='upper left')
    
    norm = plt.Normalize(df['density'].min(), df['density'].max())
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, label='Point Density')
    
    ax.set_xlabel('Query Protein Length', fontsize=12)
    ax.set_ylabel('Subject Protein Length', fontsize=12)
    ax.set_xlim(min_len*0.9, max_len*1.1)
    ax.set_ylim(min_len*0.9, max_len*1.1)
    ax.set_aspect('equal')
    
    plt.title('Protein Lengths Scatterplot', fontweight='bold')
    plt.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_base}.pdf', bbox_inches='tight')
    plt.close()

def main():
    args = parse_arguments()
    df = load_data(args.input)
    create_plot(df, args.output, args.ratio)

if __name__ == '__main__':
    main()