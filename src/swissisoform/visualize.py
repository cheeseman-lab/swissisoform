import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle

class GenomeVisualizer:
    def __init__(self, genome_handler):
        """
        Initialize the visualizer with a genome handler.
        
        Args:
            genome_handler: GenomeHandler instance providing access to genomic data
        """
        self.genome = genome_handler
        self.feature_colors = {
            'exon': '#4CAF50',    # green
            'CDS': '#2196F3',     # blue
            '5UTR': '#FFC107',    # yellow
            '3UTR': '#FFC107',    # yellow
            'transcript': '#666666', # gray
            'start_codon': '#FF0000', # red
            'stop_codon': '#FF0000'   # red
        }
        self.track_height = 0.5
        self.spacing = 0.2

    def visualize_transcript(self, gene_name, transcript_id, output_file=None):
        """
        Visualize a specific transcript.
        
        Args:
            gene_name (str): Name of the gene
            transcript_id (str): ID of the transcript to visualize
            output_file (str): Optional path to save the plot
        """
        # Get transcript features using genome handler
        transcript_data = self.genome.get_transcript_features_with_sequence(transcript_id)
        features = transcript_data['features']
        
        if features.empty:
            raise ValueError(f"No features found for transcript {transcript_id}")

        # Get transcript bounds
        transcript_info = features[features['feature_type'] == 'transcript'].iloc[0]
        transcript_start = transcript_info['start']
        transcript_end = transcript_info['end']
        span = transcript_end - transcript_start

        # Create figure
        fig, ax = plt.subplots(figsize=(12, 3))
        
        # Draw baseline (transcript)
        ax.hlines(y=1, xmin=transcript_start, xmax=transcript_end, 
                 color=self.feature_colors['transcript'], linewidth=1)

        # Plot features
        for _, feature in features.iterrows():
            if feature['feature_type'] in self.feature_colors:
                # Create rectangle for feature
                width = feature['end'] - feature['start']
                height = self.track_height
                y_pos = 0.75
                
                # Adjust height and position for different feature types
                if feature['feature_type'] in ['5UTR', '3UTR']:
                    height = self.track_height * 0.6
                    y_pos = 0.85
                
                rect = Rectangle((feature['start'], y_pos), width, height,
                               facecolor=self.feature_colors[feature['feature_type']])
                ax.add_patch(rect)

        # Add strand arrow
        strand = transcript_info['strand']
        if strand == '+':
            ax.arrow(transcript_end, 1, -span*0.02, 0, 
                    head_width=0.1, head_length=span*0.02, fc='black', ec='black')
        else:
            ax.arrow(transcript_start, 1, span*0.02, 0, 
                    head_width=0.1, head_length=span*0.02, fc='black', ec='black')

        # Add labels
        ax.text(transcript_start - (span * 0.1), 1, 
               f"{gene_name}\n{transcript_id}", 
               verticalalignment='center')

        # Remove y-axis
        ax.set_yticks([])
        
        # Format x-axis
        ax.set_xlim(transcript_start - (span * 0.05), transcript_end + (span * 0.05))
        ax.ticklabel_format(style='plain', axis='x')
        plt.xticks(rotation=45)

        # Add legend
        legend_elements = [plt.Rectangle((0,0), 1, 1, facecolor=color, label=feat)
                         for feat, color in self.feature_colors.items()
                         if feat in features['feature_type'].unique()]
        ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc='upper left')

        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, bbox_inches='tight', dpi=300)
            plt.close()
        else:
            plt.show()

    def visualize_gene_all_transcripts(self, gene_name, output_prefix):
        """
        Create visualizations for all transcripts of a gene.
        
        Args:
            gene_name (str): Name of the gene to visualize
            output_prefix (str): Prefix for output files
        """
        transcript_info = self.genome.get_transcript_ids(gene_name)
        
        for idx, row in transcript_info.iterrows():
            transcript_id = row['transcript_id']
            output_file = f"{output_prefix}_{transcript_id}.png"
            self.visualize_transcript(gene_name, transcript_id, output_file)