import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

class GenomeVisualizer:
    def __init__(self, genome_handler):
        self.genome = genome_handler
        self.feature_colors = {
            'exon': '#4CAF50',      # green
            'CDS': '#2196F3',       # blue
            'UTR': '#FFA500',       # orange
            'start_codon': '#FF0000', # red
            'stop_codon': '#800080',  # purple
            'truncation': '#FF1493',  # deep pink
            'alternative_start': '#FFD700'  # yellow
        }
        self.track_height = 0.15
        self.codon_height = 0.2

    def visualize_transcript(self, gene_name, transcript_id, alt_features=None, output_file=None):
        """
        Visualize transcript with optional alternative start sites.
        
        Args:
            gene_name (str): Name of the gene
            transcript_id (str): Transcript ID to visualize
            alt_features (pd.DataFrame, optional): Alternative start features
            output_file (str, optional): Path to save the visualization
        """
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
        fig, ax = plt.subplots(figsize=(15, 3))
        
        # Draw base transcript line (thicker)
        ax.hlines(y=0.4, xmin=transcript_start, xmax=transcript_end,
                 color='black', linewidth=1.5)

        # Plot regular transcript features
        for _, feature in features.iterrows():
            width = feature['end'] - feature['start']
            
            if feature['feature_type'] == 'exon':
                rect = Rectangle((feature['start'], 0.325),
                               width,
                               self.track_height,
                               facecolor=self.feature_colors['exon'],
                               alpha=0.3)
                ax.add_patch(rect)
            
            elif feature['feature_type'] == 'CDS':
                rect = Rectangle((feature['start'], 0.325),
                               width,
                               self.track_height,
                               facecolor=self.feature_colors['CDS'])
                ax.add_patch(rect)
            
            elif feature['feature_type'] in ['5UTR', '3UTR']:
                rect = Rectangle((feature['start'], 0.6),
                               width,
                               self.track_height,
                               facecolor=self.feature_colors['UTR'])
                ax.add_patch(rect)
            
            elif feature['feature_type'] in ['start_codon', 'stop_codon']:
                rect = Rectangle((feature['start'], 0.3),
                               width,
                               self.codon_height,
                               facecolor=self.feature_colors[feature['feature_type']])
                ax.add_patch(rect)

        # Plot alternative start sites if provided
        if alt_features is not None and not alt_features.empty:
            for _, alt_feature in alt_features.iterrows():
                alt_start = alt_feature['start']
                
                # Draw truncation bracket
                alt_end = alt_feature['end']
                bracket_height = 0.05  # Height of the vertical parts of bracket
                bracket_y = 0.2  # Base y-position of the bracket
                
                # Draw horizontal line of bracket
                ax.hlines(y=bracket_y, 
                         xmin=alt_start-0.5,
                         xmax=alt_end+0.5,
                         color=self.feature_colors['truncation'],
                         linewidth=1.5,
                         label='Truncation')
                
                # Draw vertical parts of bracket
                ax.vlines(x=[alt_start, alt_end-0.5],
                         ymin=bracket_y,
                         ymax=bracket_y + bracket_height,
                         color=self.feature_colors['truncation'],
                         linewidth=0.5)
                
                # Draw alternative start marker (vertical yellow line)
                # Calculate vertical position for alternative start marker
                black_bar_y = 0.4  # y-position of the black bar
                alt_start_height = self.codon_height  # match codon height
                alt_start_ymin = black_bar_y - (alt_start_height / 2)  # center around black bar
                alt_start_ymax = black_bar_y + (alt_start_height / 2)
                
                ax.vlines(x=alt_end,
                         ymin=alt_start_ymin,
                         ymax=alt_start_ymax,
                         color=self.feature_colors['alternative_start'],
                         linewidth=2,
                         label='Alternative start')
                
                # Add start codon label if available (now above the line)
                if 'start_codon' in alt_feature:
                    plt.text(alt_end - span*0.01, alt_start_ymax + 0.05, 
                            alt_feature['start_codon'],
                            fontsize=8,
                            rotation=45)

        # Center title
        plt.title(f"{gene_name} - {transcript_id}", pad=20, y=1.05)

        # Customize plot
        ax.set_ylim(0, 1)
        ax.set_xlim(transcript_start - 10, transcript_end + 10)
        
        # Calculate dynamic tick interval based on span
        tick_span = transcript_end - transcript_start
        if tick_span > 100000:
            tick_interval = 10000
        elif tick_span > 50000:
            tick_interval = 5000
        elif tick_span > 10000:
            tick_interval = 1000
        else:
            tick_interval = 500

        # Calculate tick positions
        base_position = transcript_start - (transcript_start % tick_interval)
        tick_positions = range(base_position, transcript_end + tick_interval, tick_interval)
        
        plt.xticks(tick_positions, [f"{pos:,}" for pos in tick_positions])
        plt.xticks(rotation=45)
        
        # Remove y-axis
        ax.set_yticks([])

        # Add legend without frame
        legend_elements = [
            plt.Rectangle((0,0), 1, 1, facecolor=self.feature_colors['exon'], 
                         alpha=0.3, label='Exon'),
            plt.Rectangle((0,0), 1, 1, facecolor=self.feature_colors['CDS'], 
                         label='CDS'),
            plt.Rectangle((0,0), 1, 1, facecolor=self.feature_colors['UTR'], 
                         label='UTR'),
            plt.Rectangle((0,0), 1, 1, facecolor=self.feature_colors['start_codon'], 
                         label='Start codon'),
            plt.Rectangle((0,0), 1, 1, facecolor=self.feature_colors['stop_codon'], 
                         label='Stop codon')
        ]
        
        # Add alternative start elements to legend if relevant
        if alt_features is not None and not alt_features.empty:
            legend_elements.extend([
                Line2D([0], [0], color=self.feature_colors['truncation'], 
                      label='Truncation', linewidth=1.5),
                Line2D([0], [0], color=self.feature_colors['alternative_start'], 
                      label='Alternative start', linewidth=2)
            ])

        ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), 
                 loc='upper left', borderaxespad=0., frameon=False)

        plt.tight_layout()
        
        # Remove spines        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)

        if output_file:
            plt.savefig(output_file, bbox_inches='tight', dpi=300,
                       facecolor='white', edgecolor='none')
            plt.close()
        else:
            plt.show()