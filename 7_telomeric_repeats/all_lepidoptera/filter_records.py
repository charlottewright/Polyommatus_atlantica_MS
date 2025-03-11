
import pandas as pd
import sys

input_file=sys.argv[1]
output_file=sys.argv[2]

# Read the TSV file into a DataFrame
df = pd.read_csv(input_file, sep='\t')

# Filter the table to keep only rows where the "Assembly Level" is "Chromosome" or "Complete" 
filtered_df = df[(df['Assembly_Level'] == 'Chromosome') | (df['Assembly_Level'] == 'Complete')]
# Filter to keep the row with the highest "Scaffold N50" for each "Organism Name"
filtered_df = filtered_df.loc[filtered_df.groupby('Organism_Name')['Assembly_Stats_Scaffold_N50'].idxmax()]

# Write the resulting DataFrame to a new TSV file
filtered_df.to_csv(output_file, sep='\t', index=False)

print(f"Filtered data has been written to {output_file}")
