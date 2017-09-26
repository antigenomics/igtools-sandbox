import argparse
import random
import numpy as np
import multiprocessing as mp


def read(str, er):
    np.random.seed()
    index = str.split(':')[0]
    seq = str.split(':')[1]
    errors =  np.random.random(size=len(seq)) < er
    new_seq =  ''.join([seq[i] if not errors[i] else np.random.choice(['A', 'C', 'G', 'T']) for i in range(len(errors))])
    quality_char = chr(int(-10 * np.log10(er) + 33))
    return '@' + index + '\n' + new_seq + '\n+\n' + quality_char*len(seq) + '\n'


parser = argparse.ArgumentParser(description="Simulate sequensing process on generated IGH sequences.")
parser.add_argument("-e", dest="error_rate", type=float)
parser.add_argument("-i", dest="input", type=str)
parser.add_argument("-d", dest="depth", type=int)

args = parser.parse_args()

# read mutated IGH repertoire
r = open(args.input)
seqs = r.read().splitlines()
r.close()

all_reads = []

for i in range(args.depth):
    pool = mp.Pool()
    reads = [pool.apply(read, args=(s, args.error_rate)) for s in seqs]
    all_reads += reads

# write results
o = open('reads.fastq', 'w')
o.write(''.join(all_reads))
o.close()