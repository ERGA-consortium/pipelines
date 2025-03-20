import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import re
import argparse

def parse_assessment_data(filename):
    data = {
        'completeness': {},
        'consistency': {},
        'hogs': 0,
        'proteins': 0
    }
    
    current_section = None
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            if line.startswith('COMPLETENESS ASSESSMENT'):
                current_section = 'completeness'
            elif line.startswith('CONSISTENCY ASSESSMENT'):
                current_section = 'consistency'
            elif line.startswith('Number of conserved HOGs'):
                data['hogs'] = int(re.search(r'\d+', line).group())
            elif line.startswith('Number of proteins in the whole proteome'):
                data['proteins'] = int(re.search(r'\d+', line).group())
            elif ':' in line and current_section:
                parts = re.split(r':|\(|\)', line)
                category = parts[0].strip()
                values = [x.strip() for x in parts[1:] if x]
                
                if len(values) >= 2:
                    count = int(re.search(r'\d+', values[0]).group())
                    percent = float(re.search(r'\d+\.\d+', values[1]).group())
                    data[current_section][category] = (count, percent)
    
    return data

def create_combined_plot(data, basename):
    sns.set_theme(style="whitegrid", font_scale=1.1)
    palette = sns.color_palette("husl", 8)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    fig.suptitle('OMARK Assessment Result', y=1.02, fontsize=18, fontweight='bold')
    
    # Completeness Plot
    comp_categories = ['Single', 'Duplicated', 'Duplicated, Unexpected', 
                      'Duplicated, Expected', 'Missing']
    comp_df = pd.DataFrame([{
        'Category': cat.replace(', ', ',\n'),
        'Percentage': data['completeness'][cat][1],
        'Count': data['completeness'][cat][0]
    } for cat in comp_categories])
    
    sns.barplot(x='Category', y='Percentage', hue='Category', data=comp_df, ax=ax1, 
                palette=palette[:5], edgecolor='black', linewidth=1, legend=False, dodge=False)
    ax1.set_title(f'Completeness Assessment\n({data["hogs"]:,} conserved HOGs)', 
                 fontsize=14, fontweight='bold')
    ax1.set_ylim(0, 100)
    ax1.set_ylabel('Percentage (%)', fontweight='bold')
    ax1.set_xlabel('')
    
    # Add annotations
    for p in ax1.patches:
        ax1.annotate(f"{p.get_height():.1f}%\n({comp_df['Count'].iloc[int(p.get_x() + p.get_width()/2)]:,})",
                    (p.get_x() + p.get_width() / 2., p.get_height() + 2),
                    ha='center', va='center', fontsize=10, color='black')

    # Consistency Plot
    cons_categories = [
        'Total Consistent', 'Consistent, partial hits', 'Consistent, fragmented',
        'Total Inconsistent', 'Inconsistent, partial hits', 'Inconsistent, fragmented',
        'Total Contaminants', 'Total Unknown'
    ]
    cons_labels = [
        'Consistent', 'Consistent\npartial', 'Consistent\nfragmented',
        'Inconsistent', 'Inconsistent\npartial', 'Inconsistent\nfragmented',
        'Contaminants', 'Unknown'
    ]
    
    cons_df = pd.DataFrame([{
        'Category': label,
        'Percentage': data['consistency'][cat][1],
        'Count': data['consistency'][cat][0]
    } for cat, label in zip(cons_categories, cons_labels)])
    
    sns.barplot(x='Category', y='Percentage', hue='Category', data=cons_df, ax=ax2,
                palette=palette, edgecolor='black', linewidth=1, legend=False, dodge=False)
    ax2.set_title(f'Consistency Assessment\n({data["proteins"]:,} total proteins)', 
                fontsize=14, fontweight='bold')
    ax2.set_ylim(0, 100)
    ax2.set_ylabel('')
    ax2.set_xlabel('')
    
    # Add annotations
    for p in ax2.patches:
        if p.get_height() > 0:
            ax2.annotate(f"{p.get_height():.1f}%\n({cons_df['Count'].iloc[int(p.get_x() + p.get_width()/2)]:,})",
                        (p.get_x() + p.get_width() / 2., p.get_height() + 2),
                        ha='center', va='center', fontsize=10, color='black')

    # Rotate x-labels and adjust layout
    for ax in [ax1, ax2]:
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3)
        for spine in ax.spines.values():
            spine.set_color('black')
            spine.set_linewidth(0.5)

    plt.tight_layout(pad=3)
    
    # Save outputs
    plt.savefig(f"{basename}.png", dpi=300, bbox_inches='tight')
    plt.savefig(f"{basename}.pdf", bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Generate combined assessment plots from analysis results')
    parser.add_argument('-i', '--input', required=True, help='Input text file with assessment results')
    parser.add_argument('-o', '--output', default='Assessment_plot', 
                       help='Output basename for the plot (default: Assessment_plot)')
    args = parser.parse_args()
    
    data = parse_assessment_data(args.input)
    create_combined_plot(data, args.output)

if __name__ == "__main__":
    main()