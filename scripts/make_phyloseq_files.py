#!/usr/bin/env python3
import pandas as pd
import glob
import os
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Combine Lemur classification outputs into OTU and Taxonomy tables for phyloseq using Target_ID as the OTU identifier."
    )
    parser.add_argument(
        "-i", "--input", required=True,
        help="Input directory containing the Lemur output files (tab-delimited)."
    )
    parser.add_argument(
        "--otu", default="otu_table.tsv",
        help="Output filename for the OTU table (default: otu_table.tsv)"
    )
    parser.add_argument(
        "--tax", default="tax_table.tsv",
        help="Output filename for the Taxonomy table (default: tax_table.tsv)"
    )
    parser.add_argument(
        "--pattern", default="*.tsv",
        help="File pattern to match input files (default: *.tsv)"
    )
    args = parser.parse_args()

    # Dictionaries to hold OTU abundance data and taxonomy info.
    otu_dict = {}  # Mapping Target_ID -> {sample: summed estimated counts}
    tax_dict = {}  # Mapping Target_ID -> taxonomy information (dictionary)
    sample_names = []

    file_list = glob.glob(os.path.join(args.input, args.pattern))
    if not file_list:
        print("No files found matching the pattern in the input directory.")
        return

    # Define the required taxonomy columns that are expected in each input file.
    # (These names should match the headers in your Lemur outputs.)
    required_tax_cols = ["species", "genus", "family", "order", "class", "phylum", "superkingdom"]
    # Define the desired output taxonomy column names and which input column provides each.
    desired_columns = {
        "Kingdom": "superkingdom",
        "Phylum": "phylum",
        "Class": "class",
        "Order": "order",
        "Family": "family",
        "Genus": "genus",
        "Species": "species"
    }

    for file in file_list:
        try:
            df = pd.read_csv(file, sep="\t", header=0)
        except Exception as e:
            print(f"Error reading {file}: {e}")
            continue

        # Verify required columns.
        if "Target_ID" not in df.columns:
            print(f"Column 'Target_ID' not found in {file}. Skipping this file.")
            continue
        if "estimated counts" not in df.columns:
            print(f"Column 'estimated counts' not found in {file}. Skipping this file.")
            continue

        # Use the file's basename (without extension) as the sample name.
        sample_name = os.path.splitext(os.path.basename(file))[0]
        sample_names.append(sample_name)

        # Check if all required taxonomy columns exist in the file.
        missing = [col for col in required_tax_cols if col not in df.columns]
        if missing:
            print(f"Missing taxonomy columns {missing} in file {file}. Skipping taxonomy info for this file.")
            continue

        # Group by Target_ID to sum the estimated counts for this sample.
        grouped = df.groupby("Target_ID", as_index=False)["estimated counts"].sum()

        # Get taxonomy info for each OTU by taking the first occurrence for the required taxonomy columns.
        tax_info = df.groupby("Target_ID", as_index=False)[required_tax_cols].first().set_index("Target_ID")

        # Reformat taxonomy info to the desired columns (and order).
        tax_info_new = pd.DataFrame(index=tax_info.index)
        for new_col, orig_col in desired_columns.items():
            tax_info_new[new_col] = tax_info[orig_col]

        # Update dictionaries with the data from this file.
        for _, row in grouped.iterrows():
            otu = row["Target_ID"]
            count_val = row["estimated counts"]
            if otu not in otu_dict:
                otu_dict[otu] = {}
            otu_dict[otu][sample_name] = count_val

            # Save taxonomy info only the first time this OTU is encountered.
            if otu not in tax_dict and otu in tax_info_new.index:
                tax_dict[otu] = tax_info_new.loc[otu].to_dict()

    # Create OTU table DataFrame.
    otu_df = pd.DataFrame.from_dict(otu_dict, orient="index")
    # Reindex columns to ensure a consistent sample order.
    otu_df = otu_df.reindex(columns=sample_names).fillna(0)
    otu_df.index.name = "Target_ID"

    # Create taxonomy table DataFrame.
    tax_df = pd.DataFrame.from_dict(tax_dict, orient="index")
    tax_df.index.name = "Target_ID"

    # Write the tables to tab-delimited files.
    otu_df.to_csv(args.otu, sep="\t")
    tax_df.to_csv(args.tax, sep="\t")
    print(f"OTU table written to {args.otu}")
    print(f"Taxonomy table written to {args.tax}")

if __name__ == "__main__":
    main()