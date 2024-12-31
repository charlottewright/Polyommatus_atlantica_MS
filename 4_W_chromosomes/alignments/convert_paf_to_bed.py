import csv
import sys
def paf_to_bed(paf_file, bed_file):
    with open(paf_file, 'r') as paf, open(bed_file, 'w') as bed:
        reader = csv.reader(paf, delimiter='\t')
        writer = csv.writer(bed, delimiter='\t', lineterminator='\n')
        
        for row in reader:
            target_name = row[5]
            target_start = row[7]
            target_end = row[8]
            writer.writerow([target_name, target_start, target_end])

def filter_paf(paf_file, overlapping_file, output_file):
    overlapping_intervals = set()

    with open(overlapping_file, 'r') as overlaps:
        reader = csv.reader(overlaps, delimiter='\t')
        for row in reader:
            overlapping_intervals.add((row[0], row[1], row[2]))

    with open(paf_file, 'r') as paf, open(output_file, 'w') as output:
        reader = csv.reader(paf, delimiter='\t')
        writer = csv.writer(output, delimiter='\t', lineterminator='\n')

        for row in reader:
            target_name = row[5]
            target_start = row[7]
            target_end = row[8]
            if (target_name, target_start, target_end) not in overlapping_intervals:
                writer.writerow(row)

# Run script
paf_file = sys.argv[1]
bed_file = sys.argv[2]
overlapping_file = "overlapping.bed"
output_file = "filtered_output.paf"

paf_to_bed(paf_file, bed_file)
