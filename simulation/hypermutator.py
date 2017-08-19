import argparse
import random
import numpy as np


class Clonotype:

    def __init__(self, sequence, index, parent, mutations):
        self.seq = sequence
        self.i = index
        self.p = parent
        self.muts = mutations

    def mutate(self, model, rep, new_rep):

        mutations = []
        new_seq = self.seq
        indels = []

        for i in range(len(self.seq)):
            if random.random() > model['substitution.ratio']:
                nucl = self.seq[i]
                new_nucl = np.random.choice(model['substitution.pattern'][nucl][0],
                                            p = model['substitution.pattern'][nucl][1])
                new_seq = new_seq[:i] + new_nucl + new_seq[i + 1:]
                mutations.append('S' + str(i) + ':' + self.seq[i] + '>' + new_nucl)

            indel = np.random.choice(model['indel.pattern'][0], p = model['indel.pattern'][1])
            if indel != 0:
                indels.append([i, indel])
                mutations.append(('D' if indel < 0 else 'I') + str(i) + ':' + str(indel).strip('-'))

        shift = 0
        for i in indels:
            if i[1] < 0:
                new_seq = new_seq[:i[0] + shift] + new_seq[i[0] + i[1] + shift:]
            else:
                insertion = ''.join(map(str, np.random.choice(['A', 'C', 'G', 'T'], i[1])))
                new_seq = new_seq[:i[0] + shift] + insertion + new_seq[i[0] + shift:]
            shift += i[1]

        if mutations:
            new_rep.append(Clonotype(new_seq, len(rep + new_rep) + 1, self.i, mutations))


def parse_model():

    f = open(args.model)
    mr = f.read().strip('#').split('#')
    f.close()
    model = {}
    for i in mr:
        params = i.strip('\n').split('\n')
        model[params[0]] = params[1:]

    model['substitution.ratio'] = float(model['substitution.ratio'][0])

    model['substitution.pattern'] = {['A', 'C', 'G', 'T'][i] :
                                         [['A', 'C', 'G', 'T'][:i] + ['A', 'C', 'G', 'T'][i+1:],
                                          model['substitution.pattern'][i].split()] for i in range(4)}
    model['substitution.pattern']['N'] = [['A', 'C', 'G', 'T'], [0.25] * 4]

    for i in ['A', 'C', 'G', 'T']:
        model['substitution.pattern'][i][1] = [float(x) for x in model['substitution.pattern'][i][1]]

    model['indel.pattern'] = [x.split() for x in model['indel.pattern']]
    model['indel.pattern'][0] = [int(x) for x in model['indel.pattern'][0]]
    model['indel.pattern'][1] = [float(x) for x in model['indel.pattern'][1]]

    return model


parser = argparse.ArgumentParser(description="Introduce hypermutations into generated IGH sequences.")
parser.add_argument("-m", dest="model", type=str)
parser.add_argument("-s", dest="seqs", type=str)
parser.add_argument("-t", dest="tree_height", type=int)

args = parser.parse_args()

# read SHM pattern
m = parse_model()

# read generated IGH repertoire
r = open(args.seqs)
rep = []
for line in r.read().splitlines()[1:]:
    index = int(line.split(';')[0]) + 1
    seq = line.split(';')[1]
    rep.append(Clonotype(seq, index, index, []))

# go through repertoire several times and introduce mutations
new_rep = []
for i in range(args.tree_height):
    for c in rep:
        c.mutate(m, rep, new_rep)
    rep = rep + new_rep

# write results
f1 = open('mutated_seqs.csv', 'w')
f1.write( '\n'.join( ['{}:{}'.format(x.i, x.seq) for x in rep] ) )
f1.close()

f2 = open('mutated_seqs_info.csv', 'w')
f2.write('index\tparent\tsequence\tmutations\n')
f2.write( ''.join( ['{}\t{}\t{}\t{}\n'.format(x.i, x.p, x.seq, ','.join(x.muts)) for x in rep] ) )
f2.close()