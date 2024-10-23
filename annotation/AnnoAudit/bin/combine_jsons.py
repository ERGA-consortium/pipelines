import os
import json
import pandas as pd

def load_json_files(json_files):
    data = {}
    for json_file in json_files:
        with open(json_file, 'r') as f:
            data[os.path.basename(json_file)] = json.load(f)
    return data

def extract_data(data):
    extracted_data = {}

    for filename, content in data.items():
        for section, stats in content.items():
            for key, value in stats.items():
                if key not in extracted_data:
                    extracted_data[key] = {}
                extracted_data[key][filename] = value

    return extracted_data

def create_dataframe(extracted_data):
    df = pd.DataFrame.from_dict(extracted_data, orient='index')
    return df

def save_to_files(df, tsv_file, excel_file):
    df.to_csv(tsv_file, sep='\t')
    df.to_excel(excel_file, index=True)

def main(json_folder, tsv_file, excel_file):
    # Get list of all JSON files in the specified folder
    json_files = [os.path.join(json_folder, file) for file in os.listdir(json_folder) if file.endswith('.json')]
    
    # Load all JSON files
    data = load_json_files(json_files)
    
    # Extract data into a structured format
    extracted_data = extract_data(data)
    
    # Create a DataFrame from the structured data
    df = create_dataframe(extracted_data)
    
    # Save the DataFrame to files
    save_to_files(df, tsv_file, excel_file)

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Combine data from multiple JSON files into a table.')
    parser.add_argument('json_folder', help='Folder containing JSON files')
    parser.add_argument('tsv_file', help='Output TSV file for the combined table')
    parser.add_argument('excel_file', help='Output Excel file for the combined table')
    
    args = parser.parse_args()
    
    main(args.json_folder, args.tsv_file, args.excel_file)


