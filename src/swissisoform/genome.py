from Bio import SeqIO
import pandas as pd

class GenomeHandler:
    def __init__(self, genome_path, gtf_path=None):
        """Initialize with genome FASTA and optional GTF path."""
        self.genome = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
        self.gtf_path = gtf_path
        if gtf_path:
            self.load_annotations(gtf_path)
    
    def load_annotations(self, gtf_path):
        """Load and parse GTF annotations into a pandas DataFrame."""
        features_list = []
        
        with open(gtf_path) as handle:
            for line in handle:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                
                if len(fields) != 9:
                    continue
                    
                # Parse attributes
                attrs = {}
                attribute_string = fields[8]
                attr_pairs = [x.strip() for x in attribute_string.rstrip(';').split(';')]
                
                for attr in attr_pairs:
                    if attr:
                        parts = attr.strip().split(' "')
                        if len(parts) == 2:
                            key = parts[0].strip()
                            value = parts[1].rstrip('"').strip()
                            attrs[key] = value
                
                feature_dict = {
                    'chromosome': fields[0],
                    'source': fields[1],
                    'feature_type': fields[2],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'score': fields[5],
                    'strand': fields[6],
                    'frame': fields[7],
                    **attrs
                }
                
                features_list.append(feature_dict)
        
        # Convert to DataFrame
        df = pd.DataFrame(features_list)
        
        # Reorder columns
        base_cols = ['chromosome', 'source', 'feature_type', 'start', 'end', 
                    'score', 'strand', 'frame']
        other_cols = [col for col in df.columns if col not in base_cols]
        self.annotations = df[base_cols + other_cols]
    
    def find_gene_features(self, gene_name):
        """Find all features associated with a gene name."""
        if not hasattr(self, 'annotations'):
            raise ValueError("No GTF file loaded")
        
        return self.annotations[self.annotations['gene_name'] == gene_name]
    
    def inspect_genome(self):
        """Print genome file structure."""
        for id, record in self.genome.items():
            print(f"Chromosome: {id}")
            print(f"Length: {len(record.seq)} bp")
    
    def inspect_annotations(self):
        """Print GTF structure using GFFExaminer."""
        if not hasattr(self, 'gtf_path'):
            raise ValueError("No GTF file path specified")
        
        from BCBio.GFF import GFFExaminer
        with open(self.gtf_path) as handle:
            examiner = GFFExaminer()
            print("\nFeature relationships:")
            print(examiner.parent_child_map(handle))
            handle.seek(0)
            print("\nFeature counts:")
            print(examiner.available_limits(handle))
    
    def get_sequence(self, chrom, start, end, strand='+'):
        """Get sequence for a genomic region."""
        if chrom not in self.genome:
            raise ValueError(f"Chromosome {chrom} not found in genome")
            
        seq = self.genome[chrom].seq[start-1:end]
        return seq if strand == '+' else seq.reverse_complement()
    
    def get_transcript_features(self, transcript_id):
        """Get all features associated with a transcript ID."""
        if not hasattr(self, 'annotations'):
            raise ValueError("No GTF file loaded")
        
        return self.annotations[self.annotations['transcript_id'] == transcript_id]
    
    def get_gene_stats(self, gene_name):
        """Get basic statistics for a gene."""
        features = self.find_gene_features(gene_name)
        if features.empty:
            return None
            
        stats = {
            'transcripts': len(features[features['feature_type'] == 'transcript']),
            'exons': len(features[features['feature_type'] == 'exon']),
            'cds': len(features[features['feature_type'] == 'CDS']),
            'utrs': len(features[features['feature_type'].isin(['5UTR', '3UTR'])]),
            'chromosomes': features['chromosome'].unique().tolist(),
            'total_span': features['end'].max() - features['start'].min() + 1
        }
        
        return stats
    
    def get_transcript_ids(self, gene_name):
        """Get all transcript IDs for a given gene name."""
        if not hasattr(self, 'annotations'):
            raise ValueError("No GTF file loaded")
        
        transcripts = self.annotations[
            (self.annotations['gene_name'] == gene_name) & 
            (self.annotations['feature_type'] == 'transcript')
        ]
        return transcripts[['transcript_id', 'chromosome', 'start', 'end', 'strand']].drop_duplicates()

    def get_transcript_sequence(self, transcript_id):
        """Get the full genomic sequence for a transcript."""
        if not hasattr(self, 'annotations'):
            raise ValueError("No GTF file loaded")
        
        # Get transcript features
        transcript = self.annotations[
            (self.annotations['transcript_id'] == transcript_id) & 
            (self.annotations['feature_type'] == 'transcript')
        ].iloc[0]
        
        # Get the sequence
        sequence = self.get_sequence(
            transcript['chromosome'],
            transcript['start'],
            transcript['end'],
            transcript['strand']
        )
        
        return {
            'transcript_id': transcript_id,
            'chromosome': transcript['chromosome'],
            'start': transcript['start'],
            'end': transcript['end'],
            'strand': transcript['strand'],
            'sequence': str(sequence)
        }

    def get_transcript_features_with_sequence(self, transcript_id):
        """Get all features and sequence for a transcript."""
        # Get transcript features
        features = self.get_transcript_features(transcript_id)
        sequence_info = self.get_transcript_sequence(transcript_id)
        
        return {
            'features': features,
            'sequence': sequence_info
        }