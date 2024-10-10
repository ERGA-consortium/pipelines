import argparse
import json
import csv

def parse_statistics_output(stats_file):
    """Reads an input file and returns a dictionary of data."""

    basic_stats = {}
    with open(stats_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            key = row[0].strip()
            value = row[1].strip()
            basic_stats[key] = value
    return basic_stats

def parse_busco_output(busco_file):
    """Parses BUSCO output and returns a dictionary of statistics."""

    with open(busco_file, 'r') as f:
        data = json.load(f)

    results = data["results"]
    busco_stats = {
        "lineage_dataset": data["lineage_dataset"]["name"],
        "complete": str(results["Complete"]) + "%",
        "single_copy": str(results["Single copy"]) + "%",
        "multi_copy": str(results["Multi copy"]) + "%",
        "fragmented": str(results["Fragmented"]) + "%",
        "missing": str(results["Missing"]) + "%",
        "num_markers": results["n_markers"],
        "domain": results["domain"],
    }

    return busco_stats 

def parse_omark_output(omark_output):
    """Parses omark output file and returns a dictionary of statistics."""

    omark_stats = {}
    with open(omark_output, 'r') as f:
        for line in f:
            # Completeness
            if line.startswith('The clade used was:'):
                omark_stats['OMA_clade'] = line.split(':')[1].strip()
            elif line.startswith('Number of conserved HOGs:'):
                omark_stats['num_conserved_hogs'] = int(line.split(':')[1].strip())
            elif line.startswith('Single:'):
                omark_stats['single'] = line.split(':')[1].strip()
            elif line.startswith('Duplicated:'):
                omark_stats['duplicated'] = line.split(':')[1].strip()
            elif line.startswith('Duplicated, Unexpected:'):
                omark_stats['duplicated_unexpected'] = line.split(':')[1].strip()
            elif line.startswith('Duplicated, Expected:'):
                omark_stats['duplicated_expected'] = line.split(':')[1].strip()
            elif line.startswith('Missing:'):
                omark_stats['missing'] = line.split(':')[1].strip()

            # Consistency
            elif line.startswith('Number of proteins in the whole proteome:'):
                omark_stats['num_proteins_in_proteome'] = line.split(':')[1].strip()
            elif line.startswith('Total Consistent:'):
                omark_stats['total_consistent'] = line.split(':')[1].strip()
            elif line.startswith('Consistent, partial hits:'):
                omark_stats['consistent_partial_hits'] = line.split(':')[1].strip()
            elif line.startswith('Consistent, fragmented:'):
                omark_stats['consistent_fragmented'] = line.split(':')[1].strip()
            elif line.startswith('Total Inconsistent:'):
                omark_stats['total_inconsistent'] = line.split(':')[1].strip()
            elif line.startswith('Inconsistent, partial hits:'):
                omark_stats['inconsistent_partial_hits'] = line.split(':')[1].strip()
            elif line.startswith('Inconsistent, fragmented:'):
                omark_stats['inconsistent_fragmented'] = line.split(':')[1].strip()
            elif line.startswith('Total Contaminants:'):
                omark_stats['total_contaminants'] = line.split(':')[1].strip()
            elif line.startswith('Contaminants, partial hits:'):
                omark_stats['contaminants_partial_hits'] = line.split(':')[1].strip()
            elif line.startswith('Contaminants, fragmented:'):
                omark_stats['contaminants_fragmented'] = line.split(':')[1].strip()
            elif line.startswith('Total Unknown:'):
                omark_stats['total_unknown'] = line.split(':')[1].strip()
    return omark_stats

def parse_brh_output(brh_output):
    """Parses BRH output and returns statistics."""

    brh_stats = {}
    num_lines = 0
    num_splitting_genes_08 = 0
    num_splitting_genes_05 = 0
    num_fusion_genes_15 = 0
    num_fusion_genes_12 = 0

    with open(brh_output, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            num_lines += 1
            qlen = int(row[2])
            slen = int(row[3])
            if slen < qlen * 0.8:
                num_splitting_genes_08 += 1
            elif slen < qlen * 0.5:
                num_splitting_genes_05 += 1
            elif slen > qlen * 1.5:
                num_fusion_genes_15 += 1
            elif slen > qlen * 1.2:
                num_fusion_genes_12 += 1

    brh_stats['num_best_reciprocal_hits'] = str(num_lines)
    brh_stats['num_splitting_genes_08'] = str(num_splitting_genes_08) + " (" + str(round(num_splitting_genes_08/num_lines*100, 2)) + "%)"
    brh_stats['num_splitting_genes_05'] = str(num_splitting_genes_05) + " (" + str(round(num_splitting_genes_05/num_lines*100, 2)) + "%)"
    brh_stats['num_fusion_genes_12'] = str(num_fusion_genes_12) + " (" + str(round(num_fusion_genes_12/num_lines*100, 2)) + "%)"
    brh_stats['num_fusion_genes_15'] = str(num_fusion_genes_15) + " (" + str(round(num_fusion_genes_15/num_lines*100, 2)) + "%)"
    return brh_stats

def parse_flagstat(file_path):

    with open(file_path, 'r') as f:
        data = json.load(f)

    passed = data["QC-passed reads"]
    rna_stats = {
        "mapping_rate": str(passed["mapped %"]) + "%",
        "primary_mapping_rate": str(passed["primary mapped %"]) + "%",
        "properly_paired": str(passed["properly paired %"]) + "%"
    }

    return rna_stats

def parse_feature(file_path):
    feature_stats = {}
    with open(file_path, "r") as f_in:
        for line in f_in:
            data = line.rstrip().split(sep="\t")
            feature_stats[data[0]] = str(data[1])
    return feature_stats

def combine_results(stats1, stats2, stats3, stats4, stats5, stats6, stats7):
    combined_stats = {**stats1, **stats2, **stats3, **stats4, **stats5, **stats6, **stats7}
    protein_combined = {**stats4, **stats5}
    rnaseq_combined = {**stats6, **stats7}

    with open("Evaluation_output.txt", "w") as f_out:
        
        max_key_length = max(len(str(key)) for key in combined_stats)
        max_value_length = max(len(str(value)) for value in combined_stats.values())
        format_string = f"|{{:<{max_key_length}}} | {{:<{max_value_length}}} |"

        print(format_string.format("General Statistics", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in stats1.items():
            print(format_string.format(key, value), file=f_out)
        print("\n", file=f_out)

        print(format_string.format("BUSCO", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in stats2.items():
            print(format_string.format(key, value), file=f_out)
        print("\n", file=f_out)

        print(format_string.format("OMARK", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in stats3.items():
            print(format_string.format(key, value), file=f_out)
        print("\n", file=f_out)

        print(format_string.format("Best Reciprocal Hits", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in protein_combined.items():
            print(format_string.format(key, value), file=f_out)
        print("\n", file=f_out)

        print(format_string.format("RNASeq", "Value"), file=f_out)
        print("-" * (max_key_length + max_value_length + 6), file=f_out)
        for key, value in rnaseq_combined.items():
            print(format_string.format(key, value), file=f_out)

    pretty_json_obj = json.dumps(combined_stats, indent=4)
    with open("Evaluation_output.json", "w") as json_out:
        json_out.write(pretty_json_obj)

def main():
    parser = argparse.ArgumentParser(description="Combine and analyze input files.")
    parser.add_argument("-s", "--statistics_output", help="Path to the general statistics file")
    parser.add_argument("-b", "--busco_output", help="Path to the busco output file")
    parser.add_argument("-o", "--omark_output", help="Path to the omark output file")
    parser.add_argument("-r", "--brh_output", help="Path to the brh output file")
    parser.add_argument("-i", "--idxstats_output", help="Path to the idxstats output file")
    parser.add_argument("-c", "--compare_distribution_output", help="Path to the compare distribution output")
    parser.add_argument("-f", "--feature_output", help="Path to the featureCounts output")
    args = parser.parse_args()

    statistics_out = parse_statistics_output(args.statistics_output)
    busco_out = parse_busco_output(args.busco_output)
    omark_out = parse_omark_output(args.omark_output)
    brh_out = parse_brh_output(args.brh_output)
    transcriptome_out = parse_flagstat(args.idxstats_output)
    prot_distribution_out = parse_statistics_output(args.compare_distribution_output)
    feature_out = parse_feature(args.feature_output)

    combine_results(statistics_out, busco_out, omark_out, brh_out, prot_distribution_out, transcriptome_out, feature_out)

if __name__ == '__main__':
    main()