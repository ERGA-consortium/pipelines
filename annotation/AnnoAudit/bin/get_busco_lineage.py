import csv
import argparse
from Bio import Entrez

def load_busco_lineages(busco_database):
    with open(busco_database, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        return {row[0].strip() for row in reader}

def get_taxonomy(taxon_id, query_email):
    Entrez.email = query_email 
    handle = Entrez.efetch(db="taxonomy", id=taxon_id, retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    lineage = records[0]["Lineage"].split("; ")
    return lineage

def get_busco_lineage(taxon_id, query_email, busco_lineages, odb_version):
    lineage = get_taxonomy(taxon_id, query_email)
    for rank in reversed(lineage):
        database = f"{rank.lower()}_{odb_version}"
        if database in busco_lineages:
            return database
    return None

def main():
    parser = argparse.ArgumentParser(description='Fetch protein sequences from NCBI for a given taxon ID.')
    parser.add_argument("-e", '--email', type=str, required=True, help='Email address to use for NCBI Entrez.')
    parser.add_argument("-t", '--taxon_id', type=int, required=True, help='The taxon ID to start the search.')
    parser.add_argument("-b", "--busco_lineage_database", type= str, required=True, help="Path to the database of current BUSCO lineages")
    parser.add_argument("-v", "--odb_version", type=str, required=True, help="odb database version 10 or 12")
    
    args = parser.parse_args()

    busco_lineages = load_busco_lineages(args.busco_lineage_database)
    taxon_id = args.taxon_id
    query_email = args.email
    odb_version = args.odb_version
    busco_lineage_output = get_busco_lineage(taxon_id, query_email, busco_lineages, odb_version)
    print(busco_lineage_output)

if __name__ == "__main__":
    main()
