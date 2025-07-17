#!/usr/bin/env python3
"""Summary analysis of SwissIsoform pipeline results.

This script analyzes mutation analysis results and localization predictions
to provide comprehensive summaries and identify genes with interesting
localization changes for all available datasets.
"""

import warnings
from pathlib import Path

# Suppress pandas warnings
warnings.filterwarnings("ignore", category=FutureWarning)

from swissisoform.summary import SummaryAnalyzer


def main():
    """Main analysis function that processes all available datasets."""
    print("SwissIsoform Pipeline Results Summary")
    print("=" * 70)

    # Initialize the analyzer
    analyzer = SummaryAnalyzer()

    # Process both datasets
    datasets = ["reduced", "full"]

    for dataset in datasets:
        print(f"\n{'=' * 20} ANALYZING {dataset.upper()} DATASET {'=' * 20}")

        # Check if this dataset has any data
        if not analyzer.dataset_has_data(dataset):
            print(f"Skipping {dataset} dataset - no data available")
            continue

        # Run analysis for this dataset
        try:
            analyzer.analyze_dataset(dataset)
            print(f"✅ {dataset} dataset analysis completed successfully!")
        except Exception as e:
            print(f"❌ Error analyzing {dataset} dataset: {e}")
            continue

    print(f"\n{'=' * 70}")
    print("🎉 Summary analysis completed for all datasets!")
    print("Results saved to: ../results/[dataset]/summary/")


if __name__ == "__main__":
    main()
