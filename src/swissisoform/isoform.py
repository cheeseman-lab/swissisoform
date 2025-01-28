import pandas as pd

class AlternativeIsoform:
    def __init__(self, bed_path):
        """Initialize isoform handler from BED file
        
        Args:
            bed_path (str): Path to BED file with alternative start sites
        """
        self.isoforms = pd.read_csv(bed_path, sep='\t', 
                                  names=['chrom', 'start', 'end', 'info', 'score', 'strand'])
        self._parse_info_field()
    
    def _parse_info_field(self):
        """Parse the info field to extract gene, start codon, and isoform type"""
        def split_info(info):
            parts = info.split('_')
            return {
                'gene_id': parts[0],
                'gene_name': parts[1],
                'start_codon': parts[2],
                'isoform_type': parts[3]
            }
        
        info_df = pd.DataFrame(self.isoforms['info'].apply(split_info).tolist())
        self.isoforms = pd.concat([self.isoforms, info_df], axis=1)
    
    def get_isoforms(self, gene_name=None, isoform_type=None):
        """Get isoforms filtered by gene name and/or type
        
        Args:
            gene_name (str): Gene name to filter by
            isoform_type (str): Isoform type to filter by
            
        Returns:
            pd.DataFrame: Filtered isoforms
        """
        result = self.isoforms
        if gene_name:
            result = result[result['gene_name'] == gene_name]
        if isoform_type:
            result = result[result['isoform_type'].str.contains(isoform_type)]
        return result