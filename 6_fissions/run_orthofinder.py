#!/bin/bash

# a = the_number_of_orthofinder_threads b- make 16 or 1/8 of threads. Use d by orthofinder which is parallised but needs a lot of RAM
# t = number_of_threads - make as much as possible

PEP_FOLDER=$1

~/Scripts/OrthoFinder/orthofinder.py -f $PEP_FOLDER -t 10 -a 4
