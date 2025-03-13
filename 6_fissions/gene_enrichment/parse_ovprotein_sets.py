import argparse
import os

# Function to check if two ranges overlap
def is_overlapping(region_start, region_end, feature_start, feature_end):
    return max(region_start, feature_start) <= min(region_end, feature_end)

# Function to parse a BED file
def parse_bed_file(bed_file):
    regions = []
    with open(bed_file, 'r') as rfile:
        for line in rfile:
            fields = line.strip().split()
            regions.append({
                'chrom': fields[0],
                'start': int(fields[1]),
                'end': int(fields[2]),
                'new': fields[3],
            })
            print(fields[0], int(fields[1]))
    return regions

# Function to parse a features file
def parse_features_file(features_file):
    features = []
    with open(features_file, 'r') as ffile:
        for line in ffile:
            fields = line.strip().split()
            features.append({
                'feat': fields[0],
                'chrom': fields[1],
                'start': int(fields[2]),
                'end': int(fields[3]),
            })
    return features

# Main function to process the input files and create the output
def process_files(regions_file, features_file, output_dir):
    # Parse input files
    regions = parse_bed_file(regions_file)
    features = parse_features_file(features_file)
    # Ensure the directory exists, create it if it doesn't
    os.makedirs(output_dir, exist_ok=True)
    # Process each region and find overlapping features
    for region in regions:
        new_file = f"{region['new']}.txt"
        output_file_path = os.path.join(output_dir, new_file)
        overlapping_features = []
        
        for feature in features:
            if region['chrom'] == feature['chrom'] and is_overlapping(region['start'], region['end'], feature['start'], feature['end']):
                overlapping_features.append(f"{feature['feat']} {feature['chrom']} {feature['start']} {feature['end']}")
        # Write the overlapping features to the new file
        if overlapping_features:
            with open(output_file_path, 'w') as outfile:
                for feature in overlapping_features:
                    outfile.write(f"{feature}\n")
        else:
            print('Chr:', region['chrom'], region['start'], region['end'], 'does not overlap with any proteins.')
# Setup argparse for command-line arguments
def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description="Process regions and features to find overlaps")
    
    # Add arguments for the input files
    parser.add_argument('--region', required=True, help="Path to the regions file (e.g., regions.bed)")
    parser.add_argument('--features', required=True, help="Path to the features file (e.g., featured.fa)")
    parser.add_argument('--outdir', required=True, help="Path to the output directory")

    # Parse the command-line arguments
    args = parser.parse_args()

    # Process the files with the provided arguments
    process_files(args.region, args.features, args.outdir)

# Run the main function if this script is executed
if __name__ == "__main__":
    main()
