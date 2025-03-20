import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Generate boxplots from psauron output CSV')
    parser.add_argument('--input', help='Path to the input CSV file')
    parser.add_argument('--output', '-o', default='PSAURON_plot', 
                       help='Base name for output files')
    args = parser.parse_args()

    # Read the CSV file, skipping the first two rows
    df = pd.read_csv(args.input, skiprows=2)

    # Separate scores by True/False status
    true_scores = df[df['psauron_is_protein'] == True]['in-frame_score']
    false_scores = df[df['psauron_is_protein'] == False]['in-frame_score']

    # Calculate statistics
    true_count = len(true_scores)
    false_count = len(false_scores)
    true_median = true_scores.median()
    false_median = false_scores.median()

    # Prepare data for Seaborn
    df_melted = df.melt(id_vars=['psauron_is_protein'], value_vars=['in-frame_score'], 
                        var_name='metric', value_name='score')
    
    # Set Seaborn theme for better aesthetics
    sns.set_theme(style="whitegrid", palette="pastel")

    # Create the plot
    # Create the plot
    plt.figure(figsize=(10, 6))
    ax = sns.boxplot(
        x='psauron_is_protein', 
        y='score', 
        hue='psauron_is_protein',
        data=df_melted, 
        width=0.5, 
        linewidth=1.5, 
        fliersize=5, 
        palette={True: '#1f77b4', False: '#ff7f0e'},  
        legend=False
    )

    # Add custom legend outside the plot
    handles = [plt.Rectangle((0,0),1,1, color='#1f77b4'), plt.Rectangle((0,0),1,1, color='#ff7f0e')]
    labels = [f'True (n={true_count}, median={true_median:.3f})', 
              f'False (n={false_count}, median={false_median:.3f})']
    plt.legend(handles, labels, title='Protein Prediction', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Formatting
    plt.title('Distribution of In-Frame Scores by PSAURON', pad=20, fontsize=14, fontweight='bold')
    plt.xlabel('Protein Prediction')
    plt.ylabel('In-Frame Score')
    plt.ylim(0, 1.05)
    plt.tight_layout()

    # Save plots
    output_base = args.output
    plt.savefig(f'{output_base}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_base}.pdf', bbox_inches='tight')
    plt.close()

if __name__ == '__main__':
    main()