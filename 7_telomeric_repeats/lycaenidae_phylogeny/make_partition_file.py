#!/usr/bin/env python3

import argparse
import os
import sys

def convert_file(input_file, output_file):
    with open(input_file, 'r') as infile:
        with open(output_file, 'w') as outfile:
            part_counter = 1
            for line in infile:
                parts = line.split('=')
                range_value = parts[1].strip()
                output_line = f"DNA, part{part_counter} = {range_value}\n"
                outfile.write(output_line)
                part_counter += 1

input_file=sys.argv[1]
output_file=sys.argv[2]

convert_file(input_file, output_file)
