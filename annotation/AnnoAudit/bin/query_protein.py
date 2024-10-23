import argparse
from Bio import Entrez
from Bio import SeqIO
import time

def fetch_protein_sequences(email, taxon_id, target_count=100000, batch_size=1000):
    Entrez.email = email
    sequences = []
    current_taxon_id = str(taxon_id)
    
    while len(sequences) < target_count and current_taxon_id is not None:
        # Initialize total fetched for current taxon ID
        total_fetched = 0
        print(f"Fetching protein of taxon ID {current_taxon_id}")
        
        while True:
            # Search for protein sequences for the current taxon ID
            try:
                search_handle = Entrez.esearch(db="protein", term=f"txid{current_taxon_id}[Organism:exp]", retstart=total_fetched, retmax=batch_size)
                search_results = Entrez.read(search_handle)
                search_handle.close()
            except Exception as e:
                print(f"Error duiring Entrez search and read: {e}")
                break
            
            # Fetch sequences in batches
            if 'IdList' in search_results:
                protein_ids = search_results['IdList']
                if not protein_ids:
                    break

                try:
                    fetch_handle = Entrez.efetch(db="protein", id=protein_ids, rettype="fasta", retmode="text")
                    sequences_batch = list(SeqIO.parse(fetch_handle, "fasta"))
                    fetch_handle.close()
                except Exception as e:
                    print(f"Error during fetch: {e}")
                    break

                sequences.extend(sequences_batch)
                total_fetched += len(protein_ids)
                
                # Sleep to avoid NCBI rate limit
                time.sleep(1)
                
                if len(sequences) >= target_count:
                    break
            else:
                break

        # If not enough sequences, move up the taxonomic tree
        if len(sequences) < target_count:
            print(f"Found {len(sequences)} protein in total. Climbing up the taxonomic tree from taxon ID {current_taxon_id}")
            try:
                taxonomy_handle = Entrez.efetch(db="taxonomy", id=current_taxon_id)
                taxonomy_record = Entrez.read(taxonomy_handle)
                taxonomy_handle.close()
            except Exception as e:
                print(f"Error fetching taxonomy information: {e}")
                break
            
            if taxonomy_record and 'LineageEx' in taxonomy_record[0]:
                lineage = taxonomy_record[0]['LineageEx']
                if lineage:
                    current_taxon_id = lineage[-1]['TaxId']
                else:
                    current_taxon_id = None
            else:
                current_taxon_id = None
    
    return sequences[:target_count]

def main():
    parser = argparse.ArgumentParser(description='Fetch protein sequences from NCBI for a given taxon ID.')
    parser.add_argument("-e", '--email', type=str, required=True, help='Email address to use for NCBI Entrez.')
    parser.add_argument("-t", '--taxon_id', type=int, required=True, help='The taxon ID to start the search.')
    parser.add_argument("-c", '--target_count', type=int, default=100000, help='The target number of protein sequences to retrieve (default: 100000).')
    parser.add_argument("-b", '--batch_size', type=int, default=10000, help='Number of sequences to fetch in each batch (default: 10000).')
    
    args = parser.parse_args()
    
    sequences = fetch_protein_sequences(args.email, args.taxon_id, args.target_count, args.batch_size)
    
    print(f"Retrieved {len(sequences)} protein sequences.")

    with open('protein_sequences.fasta', 'w') as output_handle:
        SeqIO.write(sequences, output_handle, 'fasta')

if __name__ == "__main__":
    main()
