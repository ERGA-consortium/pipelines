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

def calculate_divergence(brh_out_file):
    lengths1, lengths2 = read_columns(brh_out_file)

    # Calculate probability distributions from length
    p = [count / sum(lengths1) for count in lengths1]
    q = [count / sum(lengths2) for count in lengths2]
    
    # Calculate divergences on normalized data
    kl_div_normed = round(kl_divergence(p, q), 4)
    js_div_normed = round(js_divergence(p, q), 4)

    # Calculate divergences on raw data
    kl_div_raw = round(kl_divergence(p, q), 4)
    js_div_raw = round(js_divergence(p, q), 4)

    with open("Distribution_output.txt", "w") as f_out:
        print(f"KL_divergence_normed\t{kl_div_normed}", file=f_out)
        print(f"JS_divergence_normed\t{js_div_normed}", file=f_out)
        print(f"KL_divergence_raw\t{kl_div_raw}", file=f_out)
        print(f"JS_divergence_raw\t{js_div_raw}", file=f_out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate KL and JS divergence based on sequence lengths.")
    parser.add_argument("--brh_out", type=str, help="Path to the best reciprocal output file.")
    args = parser.parse_args()
    calculate_divergence(args.brh_out)
