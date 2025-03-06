import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def main(input_file):
    # Read the data into a DataFrame
    qlen = []
    slen = []

    # Read the file and parse the data
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            qlen.append(int(parts[2]))  # qlen is the third column
            slen.append(int(parts[3]))  # slen is the fourth column

    # Plot the distributions
    plt.figure(figsize=(12, 6))
    sns.histplot(qlen, color='blue', label='Reference_proteome', kde=True, bins=100)
    sns.histplot(slen, color='red', label='Predicted_proteome', kde=True, bins=100, alpha=0.5)
    plt.xlabel('Protein Length')
    plt.ylabel('Frequency')
    plt.title('Distribution of Protein Lengths')
    plt.legend()

    plt.savefig('Protein_distribution.pdf')
    plt.savefig('Protein_distribution.png')

    plt.xlim(0, 2000)
    plt.savefig('Protein_distribution_2000.pdf')
    plt.savefig('Protein_distribution_2000.png')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot the distribution of protein lengths from Diamond blastp output.')
    parser.add_argument("-i", '--input_file', type=str, help='Path to the input file containing Diamond blastp output.')
    
    args = parser.parse_args()
    main(args.input_file)