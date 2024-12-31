
#%%
import os

def parse_gff3(gff3_file):
    """Parse GFF3 file and return a dictionary with chromosome as keys and list of tuples (start, end, gene_name) as values."""
    genes,gene2mRNA, mRNA2CDS = {}, {}, {}
    with open(gff3_file, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue

            chrom = parts[0]
            feature_type = parts[2]
            start = int(parts[3])
            end = int(parts[4])
            attributes = parts[8]
            
            if feature_type == 'gene':
                gene_name = attributes.split(';')[0].replace("ID=gene:", "")
                if chrom not in genes:
                    genes[chrom] = []
                genes[chrom].append((start, end, gene_name))
            if feature_type == 'mRNA':
                mRNA_name = attributes.split(';')[0].replace("ID=transcript:", "")
                gene_name = attributes.split(';')[1].replace("Parent=gene:", "")
            #    print('mRNA:', mRNA_name, 'gene_name:', gene_name)
                if mRNA_name not in gene2mRNA:
                    gene2mRNA[gene_name] =mRNA_name
            if feature_type == 'CDS':
                transcript_name = attributes.split(';')[1].replace("Parent=transcript:", "")
                CDS_name = attributes.split(';')[0].replace("ID=CDS:", "")
                if transcript_name not in mRNA2CDS:
                    mRNA2CDS[transcript_name] = CDS_name
    CDS_dict = {}
    for chr,gene_set in genes.items():
        CDS_dict[chr] = []
        for start, end, gene_name in gene_set:
            mRNA = gene2mRNA[gene_name]
            CDS_name= mRNA2CDS[mRNA]
            new_tuple = (start, end, CDS_name)
            print(new_tuple)
            CDS_dict[chr].append((start, end, CDS_name))
        #cds_tuples.append(new_tuple)
    return CDS_dict

def find_genes_in_region(genes_dict, chrom, start, end):
    """Find all genes in the specified region."""
    if chrom not in genes_dict:
        return []
    
    gene_list = []
    for gene_start, gene_end, gene_name in genes_dict[chrom]:
        if gene_start <= end and gene_end >= start:
            gene_list.append(gene_name)
    return gene_list

def process_input_file(input_file, gff3_file, parent_output_dir):
    """Process the input file and save genes for each region into separate files within a parent directory."""
    genes_dict = parse_gff3(gff3_file)
    
    # Create the parent output directory if it doesn't exist
    os.makedirs(parent_output_dir, exist_ok=True)
    
    with open(input_file, 'r') as file:
        next(file)  # Skip header
        prev_chr_set_id = None
        file_count = 1
        
        for line in file:
            parts = line.strip().split()
            chr_set_id = parts[0]
            chrom = parts[1]
            start = int(parts[2])
            end = int(parts[3])
            if chr_set_id != prev_chr_set_id:
                # Reset numbering and create new directory for the chr_set_ID
                file_count = 1
                prev_chr_set_id = chr_set_id
                chr_set_dir = os.path.join(parent_output_dir, f'chr_set_{chr_set_id}')
                os.makedirs(chr_set_dir, exist_ok=True)
            
            gene_list = find_genes_in_region(genes_dict, chrom, start, end)
            output_file = os.path.join(chr_set_dir, f'sim{chr_set_id}_chr{file_count}.txt')
            with open(output_file, 'w') as out_file:
                for gene in gene_list:
                    out_file.write(f'{gene}\n')
            
            file_count += 1

#%%
# Example usage:
#input_file = '../Analysis/enrichment_analysis/test.tsv'
#parent_output_dir = '../Analysis/enrichment_analysis/test_output/'  # Replace with your desired parent output directory

input_file = '../Analysis/enrichment_analysis/P_icarus_similated_sets_227_autosomes.120924.tsv' # file containing locations of simulated chr
gff3_file = '../Analysis/enrichment_analysis/Polyommatus_icarus-GCA_937595015.1-2022_06-genes.filtered.modified.gff3'  # Replace with your GFF3 file name
parent_output_dir = '../Analysis/enrichment_analysis/simulated_chr_sets_redone/'  # Replace with your desired parent output directory

process_input_file(input_file, gff3_file, parent_output_dir)


