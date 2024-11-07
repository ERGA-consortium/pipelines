import argparse
import asyncio
import aiohttp
from Bio import Entrez, SeqIO
from io import StringIO, BytesIO
import json
import hashlib
import math

def calculate_bloom_filter_size(n: int, p: float, k: int) -> int:
    if p <= 0 or p >= 1:
        raise ValueError("Probability p must be between 0 and 1 (exclusive).")
    # Calculate size of the Bloom filter
    m = - (n * k) / math.log(1 - (p ** (1 / k)))
    return math.ceil(m)

class BloomFilter:
    def __init__(self, size: int, hash_count: int):
        self.size = size
        self.bit_array = [0] * size
        self.hash_count = hash_count

    def _hashes(self, item: str):
        """Generate hash values for the item."""
        for i in range(self.hash_count):
            hash_result = hashlib.md5(f"{item}{i}".encode()).hexdigest()
            yield int(hash_result, 16) % self.size

    def add(self, item: str):
        """Add an item to the Bloom filter."""
        for hash_value in self._hashes(item):
            self.bit_array[hash_value] = 1

    def __contains__(self, item: str) -> bool:
        """Check if an item is possibly in the Bloom filter."""
        return all(self.bit_array[hash_value] for hash_value in self._hashes(item))

async def fetch_with_retries(session, url, retries=3, delay=1):
    for i in range(retries):
        try:
            async with session.get(url) as response:
                response.raise_for_status()
                return await response.text()
        except aiohttp.ClientError as e:
            if i < retries - 1:
                await asyncio.sleep(delay)
            else:
                print(f"Failed to fetch {url}: {e}")
                return None

async def fetch_protein_sequences(email, taxon_id, target_count=100000, batch_size=6000, chunk_size=300, max_concurrent_requests=10, bloom_filter_error_rate=0.005, num_hash=5):

    Entrez.email = email
    sequences = []
    current_taxon_id = str(taxon_id)

    n = target_count
    p = bloom_filter_error_rate        
    k = num_hash            
    bloom_filter_size = calculate_bloom_filter_size(n, p, k)
    fetched_ids = BloomFilter(size=bloom_filter_size, hash_count=k)  # Initialize Bloom filter
    
    async with aiohttp.ClientSession() as session:
        while len(sequences) < target_count and current_taxon_id is not None:
            total_fetched = 0
            print(f"Fetching proteins of taxon ID {current_taxon_id}")

            while True:
                search_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                              f"db=protein&term=txid{current_taxon_id}[Organism:exp]"
                              f"&retstart={total_fetched}&retmax={batch_size}&retmode=json")
                
                print(f"Searching with URL: {search_url}")
                search_results = await fetch_with_retries(session, search_url)
                if search_results is None:
                    break

                try:
                    search_data = json.loads(search_results)
                except Exception as e:
                    print(f"Error parsing search results: {e}")
                    break
                
                if 'esearchresult' in search_data and 'idlist' in search_data['esearchresult']:
                    protein_ids = search_data['esearchresult']['idlist']
                    if not protein_ids:
                        break

                    new_ids = [pid for pid in protein_ids if pid not in fetched_ids]
                    if not new_ids:
                        break

                    # Split the new_ids into chunks to avoid "Request-URI Too Long" error
                    id_chunks = [new_ids[i:i + chunk_size] for i in range(0, len(new_ids), chunk_size)]

                    # Fetch each chunk in parallel, up to max_concurrent_requests at a time
                    chunk_tasks = []
                    for chunk in id_chunks:
                        if len(chunk_tasks) >= max_concurrent_requests:
                            results = await asyncio.gather(*chunk_tasks)
                            for fetch_results in results:
                                if fetch_results:
                                    sequences_batch = list(SeqIO.parse(StringIO(fetch_results), "fasta"))
                                    for seq_record in sequences_batch:
                                        if seq_record.id not in fetched_ids:
                                            sequences.append(seq_record)
                                            fetched_ids.add(seq_record.id)  # Using Bloom filter to track fetched IDs
                            chunk_tasks = []

                        fetch_url = (f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
                                     f"db=protein&id={','.join(chunk)}&rettype=fasta&retmode=text")
                        print(f"Fetching with URL: {fetch_url}")
                        chunk_tasks.append(fetch_with_retries(session, fetch_url))

                    # Fetch any remaining chunks
                    if chunk_tasks:
                        results = await asyncio.gather(*chunk_tasks)
                        for fetch_results in results:
                            if fetch_results:
                                sequences_batch = list(SeqIO.parse(StringIO(fetch_results), "fasta"))
                                for seq_record in sequences_batch:
                                    if seq_record.id not in fetched_ids:
                                        sequences.append(seq_record)
                                        fetched_ids.add(seq_record.id)

                    total_fetched += len(protein_ids)

                    if len(sequences) >= target_count:
                        break
                else:
                    print(f"No IdList found in search results: {search_data}")
                    break

            if len(sequences) < target_count:
                print(f"Found {len(sequences)} proteins in total. Climbing up the taxonomic tree from taxon ID {current_taxon_id}")
                taxonomy_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={current_taxon_id}&retmode=xml"
                
                taxonomy_results = await fetch_with_retries(session, taxonomy_url)
                if taxonomy_results is None:
                    break

                try:
                    taxonomy_data = Entrez.read(BytesIO(taxonomy_results.encode('utf-8')))
                except Exception as e:
                    print(f"Error parsing taxonomy results: {e}")
                    break
                
                if taxonomy_data and 'LineageEx' in taxonomy_data[0]:
                    lineage = taxonomy_data[0]['LineageEx']
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
    parser.add_argument("-b", '--batch_size', type=int, default=6000, help='Number of sequences to fetch in each batch (default: 1000).')
    parser.add_argument("-k", '--chunk_size', type=int, default=300, help='Number of sequence IDs to fetch in each request (default: 300).')
    parser.add_argument("-m", '--max_concurrent_requests', type=int, default=10, help='Maximum number of concurrent requests (default: 10).')
    parser.add_argument("-n", '--num_hash', type=int, default=5, help="Number of hash function for the Bloom filter")
    parser.add_argument("-r", '--error_rate', type=float, default=0.005, help="False positive rate for the Bloom filter")
    
    args = parser.parse_args()
    
    loop = asyncio.get_event_loop()
    sequences = loop.run_until_complete(fetch_protein_sequences(args.email, args.taxon_id, args.target_count, args.batch_size, args.chunk_size, args.max_concurrent_requests,
                                                                args.error_rate, args.num_hash))
    
    print(f"Retrieved {len(sequences)} protein sequences.")

    with open('protein_sequences.fasta', 'w') as output_handle:
        SeqIO.write(sequences, output_handle, 'fasta')

if __name__ == "__main__":
    main()
