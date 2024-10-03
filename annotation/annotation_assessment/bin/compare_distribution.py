import math
import argparse

def read_columns(filename):
    col3 = []
    col4 = []
    with open(filename, 'r') as file:
        for line in file:
            parts = line.split()
            col3.append(int(parts[2]))
            col4.append(int(parts[3]))
    return col3, col4

def kl_divergence(p, q):
    return sum(p[i] * math.log(p[i] / q[i]) for i in range(len(p)) if p[i] != 0)

def js_divergence(p, q):
    m = [(p[i] + q[i]) / 2 for i in range(len(p))]
    return 0.5 * (kl_divergence(p, m) + kl_divergence(q, m))

def cdf(data):
    """Calculate the cumulative distribution function for a dataset."""
    n = len(data)
    sorted_data = sorted(data)
    cdf_values = [(i + 1) / n for i in range(n)]
    return sorted_data, cdf_values

def wasserstein_distance_custom(distribution1, distribution2):
    """Calculate the 1D Wasserstein distance between two distributions."""
    # Get the CDFs of both distributions
    sorted_d1, cdf_d1 = cdf(distribution1)
    sorted_d2, cdf_d2 = cdf(distribution2)

    # Combine the sorted data
    combined_sorted_data = sorted(set(sorted_d1 + sorted_d2))

    # Calculate the cumulative sums for both CDFs
    cdf_d1_dict = dict(zip(sorted_d1, cdf_d1))
    cdf_d2_dict = dict(zip(sorted_d2, cdf_d2))

    def get_cdf_value(cdf_dict, x):
        """Helper function to get the CDF value for a given x."""
        keys = sorted(cdf_dict.keys())
        for key in keys:
            if x <= key:
                return cdf_dict[key]
        return 1.0

    cdf_d1_values = [get_cdf_value(cdf_d1_dict, x) for x in combined_sorted_data]
    cdf_d2_values = [get_cdf_value(cdf_d2_dict, x) for x in combined_sorted_data]

    # Calculate the Wasserstein distance
    wasserstein_dist = 0
    for i in range(1, len(combined_sorted_data)):
        delta_x = combined_sorted_data[i] - combined_sorted_data[i - 1]
        avg_cdf_diff = abs(cdf_d1_values[i] - cdf_d2_values[i] + cdf_d1_values[i - 1] - cdf_d2_values[i - 1]) / 2
        wasserstein_dist += delta_x * avg_cdf_diff

    return wasserstein_dist

def calculate_divergence(brh_out_file):
    lengths1, lengths2 = read_columns(brh_out_file)

    # Calculate probability distributions from length
    p = [count / sum(lengths1) for count in lengths1]
    q = [count / sum(lengths2) for count in lengths2]
    
    # Calculate divergences on normalized data
    kl_div_normed = round(kl_divergence(p, q), 4)
    js_div_normed = round(js_divergence(p, q), 4)
    wasserstein_dist = round(wasserstein_distance_custom(lengths1, lengths2), 6)

    with open("Distribution_output.txt", "w") as f_out:
        print(f"KL_divergence_normed\t{kl_div_normed}", file=f_out)
        print(f"JS_divergence_normed\t{js_div_normed}", file=f_out)
        print(f"Wasserstein_distance\t{wasserstein_dist}", file=f_out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate KL and JS divergence based on sequence lengths.")
    parser.add_argument("--brh_out", type=str, help="Path to the best reciprocal output file.")
    args = parser.parse_args()
    calculate_divergence(args.brh_out)
