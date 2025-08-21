#!/lustre/scratch123/tol/teams/blaxter/users/cw22/software/conda_install/miniconda3/bin/python3

import argparse

def filter_lines(input_file, output_file):
    # Define the set of valid 5-mers
    valid_kmers = {"TTAGG", "TAGGT", "AGGTT", "GGTTA", "GTTAG",
                   "CCTAA", "CTAAC", "TAACC", "AACCT", "ACCTA"}

    chr = None
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Check for 'Sequence:' line
            if line.startswith('Sequence:'):
                # Extract the value next to 'Sequence:'
                chr = line.split()[1]
                continue  # Skip the 'Sequence:' line
            
            # Skip empty lines and parameter lines
            if line.strip() == '' or line.startswith('Parameters:'):
                continue

            # Process the remaining lines
            columns = line.strip().split()
            if len(columns) >= 14 and columns[13] in valid_kmers:
                outfile.write(f"{chr} {line}") # add the chromosome identifier

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter TRF output to include the TTAGG telomeric repeat')
    parser.add_argument('input_file', type=str, help='The input TRF file')
    parser.add_argument('output_file', type=str, help='The output filtered TFT file')
    
    args = parser.parse_args()
    
    filter_lines(args.input_file, args.output_file)
