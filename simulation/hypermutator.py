import argparse
import random
import numpy as np
import multiprocessing as mp

class Clonotype:

    def __init__(self, seq, index, parent, muts, root, germline):
        self.seq = seq
        self.index = index
        self.parent = parent
        self.muts = muts
        self.root = root
        self.germline = germline


def parse_model():
    f = open(args.model)
    substitution_matrix = np.zeros((4, 4))
    for i, line in enumerate(f):
        substitution_matrix[i, :] = np.array([float(x) for x in line.split()])
    return substitution_matrix


def seq_to_array(seq):
    seq = np.array(list(seq))
    array = np.zeros(len(seq), dtype = int)
    array[seq == 'A'] = 0
    array[seq == 'C'] = 1
    array[seq == 'G'] = 2
    array[seq == 'T'] = 3
    return array


def array_to_seq(array):
    seq = np.chararray(shape=array.shape)
    seq[array == 0] = 'A'
    seq[array == 1] = 'C'
    seq[array == 2] = 'G'
    seq[array == 3] = 'T'
    return ''.join(seq)


def matrix_repr(seq):
    matrix = np.zeros((seq.shape[0], 4))
    matrix[np.arange(seq.shape[0]), seq] = 1
    return matrix


def random_choice_prob_index(a, axis=1):
    r = np.expand_dims(np.random.rand(a.shape[1-axis]), axis=axis)
    return (a.cumsum(axis=axis) > r).argmax(axis=axis)


def mutate(cl, model):
    np.random.seed()

    substitution_probas = np.dot(matrix_repr(cl.seq), model)

    new_seq = random_choice_prob_index(substitution_probas)
    mut_idxs = np.where(new_seq != cl.germline.seq)[0]

    if mut_idxs.shape[0] > 0:
        return Clonotype(seq=new_seq, index=0, parent=cl.index, muts=mut_idxs, root=cl.root, germline=cl.germline)


def write_mutations(cl):
    old = array_to_seq(cl.germline.seq[cl.muts])
    new = array_to_seq(cl.seq[cl.muts])
    muts_str = [str(idx) + ':' + x + '>' + y for idx, x, y in zip(cl.muts, old, new)]
    return muts_str


parser = argparse.ArgumentParser(description="Introduce hypermutations into generated IGH sequences.")
parser.add_argument("-m", dest="model", type=str)
parser.add_argument("-s", dest="seqs", type=str)
parser.add_argument("-t", dest="tree_height", type=int)

args = parser.parse_args()

# read SHM pattern
m = parse_model()

# read generated IGH repertoire and introduce "old" mutations to the future tree roots
rep = []
germlines = []

r = open(args.seqs)
r.readline()

for line in r:
    index = line.split(';')[0]
    seq = line.split(';')[1]
    seq = seq_to_array(seq)

    g = Clonotype(seq, index, index, np.array([]), index, germline=None)
    g.germline = g
    germlines.append(g)

    cl = mutate(g, m)
    cl.index = index
    cl.germline = g
    rep.append(cl)

# go through repertoire several times and introduce mutations
for s in range(args.tree_height):
    pool = mp.Pool()
    new_rep = [pool.apply(mutate, args=(cl, m)) for cl in rep]
    new_rep = [cl for cl in new_rep if cl is not None]
    for i, cl in enumerate(new_rep):
        cl.index = len(rep) + 1 + i
    rep = rep + new_rep

# write results
f1 = open('mutated_seqs_test.csv', 'w')

f1.write( '\n'.join( ['{}:{}'.format(x.index, array_to_seq(x.seq)) for x in rep] ) )
f1.close()

f2 = open('mutated_seqs_test_info.csv', 'w')
f2.write('index\tparent\troot\tsequence\tmutations\n')
f2.write( ''.join( ['{}\t{}\t{}\t{}\t{}\n'.format(x.index, x.parent, x.root, array_to_seq(x.seq), write_mutations(x)) for x in rep] ) )
f2.close()