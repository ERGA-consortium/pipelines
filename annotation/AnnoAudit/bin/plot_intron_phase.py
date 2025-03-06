import argparse
import pickle
import matplotlib.pyplot as plt
import seaborn as sns

def plot_distribution(seq_lengths, image_prefix):
    # Prepare the plot
    plt.figure(figsize=(12, 6))
    
    # Plot the histograms for each phase
    for phase, lengths in seq_lengths.items():
        sns.histplot(lengths, kde=True, bins=100, label=phase, fill=False)
    
    plt.xlabel('Intron Length')
    plt.ylabel('Number of Sequences')
    plt.title('Intron Length Distribution by Phase')
    plt.legend(title='Phase')
    plt.xlim(40, 121)
    
    # Save the plot as PDF and PNG
    plt.savefig(image_prefix + '.pdf')
    plt.savefig(image_prefix + '.png')

def main(short_with_stop, short_without_stop):
    with open(short_with_stop, 'rb') as handle:
        short_with_stop = pickle.load(handle)
    with open(short_without_stop, 'rb') as handle:
        short_without_stop = pickle.load(handle)

    plot_distribution(short_with_stop, "Short_intron_with_stop_distribution")
    plot_distribution(short_without_stop, "Short_intron_without_stop_distribution")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the sequence length distribution from two pickle files categorized by phase.')
    parser.add_argument("-a", '--short_with_stop', type=str, help='Path to the input pickle file.')
    parser.add_argument("-b", "--short_without_stop",  type=str, help='Path to the input pickle file.')

    args = parser.parse_args()
    main(args.short_with_stop, args.short_without_stop)