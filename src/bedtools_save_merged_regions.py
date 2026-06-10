import logging
import pandas as pd
from pathlib import Path
from typing import List
from bedtools_merge import GenomicRegion, parse_merged_bed

logger = logging.getLogger(__name__)

def save_merged_regions(regions: List[GenomicRegion], output_path: Path) -> None:
    """
    Save merged regions to TSV file.
    
    Args:
        regions: List of genomic regions
        output_path: Path to output file
    """
    logger.debug(f"Saving {len(regions)} merged regions to: {output_path}")
    
    data = [region.to_dict() for region in regions]
    df = pd.DataFrame(data)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    
    logger.debug("Merged regions saved successfully")




import sys
if __name__ == "__main__":
    genomic_regions = parse_merged_bed(Path(sys.argv[1]), genome_id=sys.argv[3])
    save_merged_regions(genomic_regions, Path(sys.argv[2]))