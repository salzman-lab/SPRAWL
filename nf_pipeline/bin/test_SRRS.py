#!/usr/local/bin/python
import SRRS
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--input_bam")
args = parser.parse_args()

print('made it here')
print(args.input_bam)

f_size = os.path.getsize(args.input_bam)

with open('num_alignments.txt','w') as f_out:
    f_out.write('{} has size {:,}\n'.format(args.input_bam,f_size))

