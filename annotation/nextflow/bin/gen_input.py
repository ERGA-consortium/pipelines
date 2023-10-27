import os
import csv
import argparse

# Function to determine if a FASTQ file is single-end or paired-end
def determine_read_type(file_name):
    if "_1.fastq.gz" in file_name:
        return "paired"
    elif "_2.fastq.gz" in file_name:
        return "paired"
    else:
        return "single"

def main(directory, output_csv):
    # Initialize a dictionary to store sample information
    sample_info = {}

    # Iterate over files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".fastq.gz") or filename.endswith(".fq.gz"):
            sample_id = filename.split(".")[0].split("_")[0]
            read_type = determine_read_type(filename)
            if sample_id not in sample_info:
                sample_info[sample_id] = {"sample_id": sample_id, "R1_path": "", "R2_path": "", "read_type": read_type}
            if read_type == "paired":
                if "_1.fastq.gz" in filename:
                    sample_info[sample_id]["R1_path"] = os.path.join(directory, filename)
                elif "_2.fastq.gz" in filename:
                    sample_info[sample_id]["R2_path"] = os.path.join(directory, filename)
            else:
                sample_info[sample_id]["R1_path"] = os.path.join(directory, filename)

    # Write the sample information to the CSV file
    with open(output_csv, mode='w', newline='') as csvfile:
        csv_writer = csv.DictWriter(csvfile, fieldnames=["sample_id", "R1_path", "R2_path", "read_type"])
        csv_writer.writeheader()
        csv_writer.writerows(sample_info.values())

    print(f"CSV file '{output_csv}' has been generated.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Parse FASTQ files in a directory and generate a CSV file.")
    parser.add_argument("directory", help="Directory containing FASTQ files")
    parser.add_argument("output_csv", help="Output CSV file")
    args = parser.parse_args()

    main(os.path.abspath(args.directory), args.output_csv)

