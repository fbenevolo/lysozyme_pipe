#!/usr/bin/env python3
"""
Script to download genome sequences from NCBI based on BVBRC metadata CSV.
Downloads FASTA files and creates a metadata file for downstream analysis.
"""
import urllib.request
import urllib.error
import pandas as pd
import sys
import subprocess
from pathlib import Path
import logging
import time
from typing import Dict, List, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def parse_bvbrc_csv(csv_path: str) -> pd.DataFrame:
    """Parse BVBRC genome metadata CSV file."""
    logger.info(f"Reading CSV file: {csv_path}")
    df = pd.read_csv(csv_path)
    logger.info(f"Found {len(df)} genomes in CSV")
    return df


def extract_accession_numbers(df: pd.DataFrame) -> pd.DataFrame:
    """Extract GenBank accession numbers from the DataFrame."""
    logger.info("Extracting accession numbers...")
    
    genome_info = []
    
    for idx, row in df.iterrows():
        genome_id = row.get('Genome ID', '')
        genome_name = row.get('Genome Name', '')
        strain = row.get('Strain', '')
        pathovar = row.get('Pathovar', '')
        genbank_accessions = row.get('GenBank Accessions', '')
        assembly_accession = row.get('Assembly Accession', '')
        
        # Parse GenBank accessions (can be comma-separated)
        if pd.notna(genbank_accessions) and genbank_accessions:
            accessions = [acc.strip() for acc in str(genbank_accessions).split(',')]
        elif pd.notna(assembly_accession) and assembly_accession:
            accessions = [assembly_accession.strip()]
        else:
            accessions = []
        
        if accessions:
            genome_info.append({
                'genome_id': genome_id,
                'genome_name': genome_name,
                'strain': strain,
                'pathovar': pathovar if pd.notna(pathovar) else 'Unknown',
                'accessions': accessions,
                'primary_accession': accessions[0]
            })
    
    result_df = pd.DataFrame(genome_info)
    logger.info(f"Extracted {len(result_df)} genomes with valid accessions")
    
    return result_df


def download_genome_efetch(accession: str, output_path: Path) -> Optional[Path]:
    """Download genome using NCBI E-utilities REST API directly via Python."""
    
    # The NCBI E-utilities endpoint for fetching sequences
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    
    try:
        # NCBI requires a User-Agent, and identifying your script is good practice
        headers = {'User-Agent': 'Python-urllib/3.x (BVBRC_Downloader)'}
        req = urllib.request.Request(url, headers=headers)
        
        with urllib.request.urlopen(req, timeout=300) as response:
            fasta_data = response.read().decode('utf-8')
            
            # Check if NCBI returned an empty response or an error string
            if not fasta_data.strip():
                logger.error(f"  Empty response from NCBI for {accession}")
                return None
            if "Error" in fasta_data[:100]: # Sometimes NCBI returns 200 OK but includes an error in text
                logger.error(f"  NCBI Error for {accession}: {fasta_data[:100].strip()}...")
                return None
                
            # Write out the FASTA file
            with open(output_path, 'w') as f:
                f.write(fasta_data)
                
        logger.info(f"  ✓ Downloaded: {output_path.name}")
        return output_path
        
    except urllib.error.HTTPError as e:
        logger.error(f"  HTTP Error for {accession}: {e.code} - {e.reason}")
        return None
    except urllib.error.URLError as e:
        logger.error(f"  Network/URL Error for {accession}: {e.reason}")
        return None
    except Exception as e:
        logger.error(f"  Download error for {accession}: {e}")
        return None

 
def download_genome_ncbi(accession: str, output_dir: Path, genome_id: str) -> Optional[Path]:
    """Download genome FASTA from NCBI."""
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Sanitize genome_id for filename - normalize dots to underscores
    safe_id = str(genome_id).replace(' ', '_').replace('/', '_').replace(':', '_').replace('.', '_')
    output_fasta = output_dir / f"{safe_id}.fasta"
    
    # Check if already downloaded
    if output_fasta.exists():
        logger.info(f"  Already exists: {output_fasta.name}")
        return output_fasta
    
    logger.info(f"  Downloading {accession}...")
    
    return download_genome_efetch(accession, output_fasta)


def create_metadata_file(genome_info: pd.DataFrame, output_path: Path):
    """Create a metadata TSV file for downstream analysis."""
    logger.info(f"Creating metadata file: {output_path}")
    
    metadata = genome_info[['genome_id', 'genome_name', 'strain', 'pathovar', 'primary_accession']].copy()
    # Ensure genome_id is string and normalize format (dots to underscores)
    metadata['genome_id'] = metadata['genome_id'].astype(str).str.replace('.', '_')
    
    # Ensure pathovar defaults to 'Unknown' if missing
    metadata['pathovar'] = metadata['pathovar'].fillna('Unknown')
    
    metadata.to_csv(output_path, sep='\t', index=False)
    
    logger.info(f"Metadata saved with {len(metadata)} genomes")
    
    # Log pathovar distribution
    pathovar_counts = metadata['pathovar'].value_counts()
    logger.info(f"\nPathovar distribution in metadata:")
    for pathovar, count in pathovar_counts.items():
        logger.info(f"  {pathovar}: {count}")


def main():
    """Main function."""
    if len(sys.argv) < 2:
        print("Usage: python download_genomes.py <BVBRC_CSV_file> [output_directory]")
        print("\nExample:")
        print("  python download_genomes.py 'BVBRC_genome (1).csv' genomes/")
        sys.exit(1)
    
    csv_path = sys.argv[1]
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("downloaded_genomes")
    
    logger.info("=" * 80)
    logger.info("BVBRC Genome Downloader")
    logger.info("=" * 80)
    
    # Parse CSV
    df = parse_bvbrc_csv(csv_path)
    
    # Extract accession info
    genome_info = extract_accession_numbers(df)
    
    if len(genome_info) == 0:
        logger.error("No valid accessions found in CSV!")
        sys.exit(1)
    
    # Create output directory
    fasta_dir = output_dir / "fasta"
    fasta_dir.mkdir(parents=True, exist_ok=True)
    
    # Download genomes
    logger.info(f"\nDownloading {len(genome_info)} genomes...")
    logger.info("-" * 80)
    
    successful = 0
    failed = 0
    
    for idx, row in genome_info.iterrows():
        genome_id = row['genome_id']
        accession = row['primary_accession']
        
        logger.info(f"[{idx+1}/{len(genome_info)}] {genome_id} ({accession})")
        
        result = download_genome_ncbi(accession, fasta_dir, genome_id)
        
        if result:
            successful += 1
        else:
            failed += 1
        
        # Rate limiting to avoid NCBI blocking
        time.sleep(1)
    
    logger.info("-" * 80)
    logger.info(f"Download complete: {successful} successful, {failed} failed")
    
    # Create metadata file
    metadata_path = output_dir / "metadata.tsv"
    create_metadata_file(genome_info, metadata_path)
    
    # Summary
    logger.info("\n" + "=" * 80)
    logger.info("SUMMARY")
    logger.info("=" * 80)
    logger.info(f"FASTA files: {fasta_dir}")
    logger.info(f"Metadata: {metadata_path}")
    logger.info(f"Total genomes: {successful}")
    
    # Pathovar statistics
    pathovar_counts = genome_info['pathovar'].value_counts()
    logger.info(f"\nPathovar distribution:")
    for pathovar, count in pathovar_counts.items():
        logger.info(f"  {pathovar}: {count}")
    
    logger.info("\nReady to run analysis:")
    logger.info(f"  python pipeline.py --input-dir {fasta_dir} -l data/lysozymes.fasta -o output/analysis")
    logger.info("=" * 80)


if __name__ == "__main__":
    main()
