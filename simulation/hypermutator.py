import random
import numpy as np
import pandas as pd
import multiprocessing as mp


class Clonotype:

    def __init__(self, seq, index, parent, root, proliferation, mutations):
        self.seq = seq
        self.index = index
        self.parent = parent
        self.root = root
        self.proliferation = proliferation
        self.mutations = mutations


def merge_repertoires(df1, df2, column_list):
    # Reset indices in the listed columns in the second df
    df1_max = df1[column_list].max().max()
    df2_min = df2[column_list].min().min()
    for c in column_list:
        df2[c] = df2[c].apply(lambda x: x - df2_min + df1_max + 1)
    df = df1.append(df2)
    df = df.reset_index(drop=True)

    return df


def correct_mut_rate(mut_rate):
    # This coefficient was calculated empirically
    return mut_rate * 0.543


def parse_model_subs(mut_rate):
    # The assessed frequency of substitutions among all mutations is 0.992
    rate = 0.992 * mut_rate
    #rate = mut_rate
    f = open('models/hypermutator/model_subs_steele.txt')
    substitution_matrix = np.zeros((4, 4))
    for i, line in enumerate(f):
        row = [float(x) * rate * 4 for x in line.split()]
        row = row[:i] + [1 - sum(row)] + row[i:]
        substitution_matrix[i, :] = np.array(row)
    return substitution_matrix


def parse_model_indels(mut_rate):
    # The assessed frequency of indels among all mutations is 0.008
    rate = 0.008 * mut_rate
    indel_matrix = np.loadtxt('models/hypermutator/model_indels.txt')
    indel_matrix[1, :] *= rate
    indel_matrix = np.c_[indel_matrix, np.array([0, 1 - rate])]
    return indel_matrix


def generate_insertion(n):
    np.random.seed()
    return np.random.choice(4, n)


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
    return ''.join(seq.decode("utf-8"))


def matrix_repr(seq):
    matrix = np.zeros((seq.shape[0], 4))
    matrix[np.arange(seq.shape[0]), seq] = 1
    return matrix


def random_choice_prob_index(a, axis=1):
    r = np.expand_dims(np.random.rand(a.shape[1-axis]), axis=axis)
    return (a.cumsum(axis=axis) > r).argmax(axis=axis)


def mutate(cl, model_subs, model_indels):
    np.random.seed()

    if np.random.random() < cl.proliferation:

        # introduce substitutions firstly
        substitution_probas = np.dot(matrix_repr(cl.seq), model_subs)

        new_seq = random_choice_prob_index(substitution_probas)
        sub_idxs = np.where(new_seq != cl.seq)[0]

        # then introduce indels
        indel_lengths = np.random.choice(model_indels[0,:], new_seq.shape[0], p=model_indels[1,:])
        indel_idxs = np.where(indel_lengths != 0)[0]
        indel_lengths = indel_lengths[indel_idxs].astype(int)

        for i in range(indel_idxs.shape[0]):
            idx = indel_idxs[i] + indel_lengths[:i].sum()
            if indel_lengths[i] > 0:
                insertion = generate_insertion(indel_lengths[i])
                new_seq = np.hstack((new_seq[:idx], insertion, new_seq[idx:]))
            else:
                new_seq = np.hstack((new_seq[:idx-1], new_seq[idx+indel_lengths[i]:]))

        mut_idxs = np.hstack((sub_idxs, indel_idxs))
        mut_num = mut_idxs.shape[0]

        if mut_num > 0:
            return Clonotype(seq=new_seq, index=0, parent=cl.index, root=cl.root,
                             proliferation = cl.proliferation, mutations=cl.mutations + mut_num)


def df_to_clonotypes(df):
    # Create a list of Clonotype objects from a data frame
    rep = []

    for i in range(len(df.index)):
        index = df['index'][i]
        seq = df['sequence'][i]
        seq = seq_to_array(seq)

        # If no proliferation probabilities provided, assign initial proliferation probabilities equal to 1,
        # because all clonotypes must mutate firstly to simulate "old" mutations
        if 'proliferation' in df:
            proliferation = df['proliferation'][i]
        else:
            proliferation = 1

        g = Clonotype(seq=seq, index=index, parent=index, root=index,
                      proliferation = proliferation, mutations=0)
        rep.append(g)

    return rep


def clonotypes_to_df(rep):
    # Create a data frame from a list of Clonotype objects
    return pd.concat([pd.DataFrame([[cl.index, cl.parent, cl.root, cl.mutations, array_to_seq(cl.seq)]],
                                   columns=['index', 'parent', 'root', 'mutations', 'sequence']) for cl in rep],
                     ignore_index=True)


def mutate_roots(rep, ms, mi):
    # Introducing "old" mutations into sequences
    # Applied to tree roots and singletons
    # Does not create new Clonotypes
    new_rep = []

    for g in rep:
        cl = mutate(g, ms, mi)
        if cl is None:
            cl = g
        cl.index = g.index
        new_rep.append(cl)
        del cl

    return new_rep


def mutate_repertoire(final_mutation_rate, rep_df):
    # Get correct mutation rate
    mutation_rate = correct_mut_rate(final_mutation_rate)
    #print('mut rate corrected')
    #print(mutation_rate)

    # Split the repertoire to two random samples (1/6 and 5/6) for tree roots and singletons
    # Note: we want final repertoire to contain half singletons and half tree nodes;
    # N tree roots give 6N tree nodes in current algorithm

    r = 5/6
    single_idxs = np.random.choice(rep_df.index, int(len(rep_df.index)*r), replace=False)
    truefalses = np.array([False] * len(rep_df.index))
    truefalses[single_idxs] = True

    df_single = rep_df[truefalses]
    df_single = df_single.reset_index(drop=True)
    #print(len(df_single.index))

    df_tree = rep_df[~truefalses]
    df_tree = df_tree.reset_index(drop=True)
    #print(len(df_tree.index))

    rep_singlet = df_to_clonotypes(df_single)
    rep_tree = df_to_clonotypes(df_tree)

    # Read SHM and indels patterns with mut.rate 0.01 and mutate tree roots
    coeff = 0.9
    ms_root = parse_model_subs(coeff * mutation_rate)
    mi_root = parse_model_indels(coeff * mutation_rate)
    rep_tree = mutate_roots(rep_tree, ms_root, mi_root)

    # Reassign proliferation probabilities that were assigned to 1 by default,
    # because all clonotypes must mutate firstly to simulate "old" mutations
    probs = np.random.triangular(left=0, mode=0, right=1, size=len(rep_tree))
    for i, cl in enumerate(rep_tree):
        cl.proliferation = probs[i]

    # Read SHM and indels patterns that will be used for singletons and trees
    ms = parse_model_subs(mutation_rate)
    mi = parse_model_indels(mutation_rate)

    # Mutate singletons
    rep_singlet = mutate_roots(rep_singlet, ms, mi)

    # Iteratively mutate rep_tree_roots creating novel tree nodes
    for s in range(5):
        pool = mp.Pool()
        new_rep = [pool.apply(mutate, args=(cl, ms, mi)) for cl in rep_tree]
        new_rep = [cl for cl in new_rep if cl is not None]
        for i, cl in enumerate(new_rep):
            cl.index = len(rep_tree) + 1 + i
        rep_tree = rep_tree + new_rep
    #print(len(rep_tree))

    # Rep to df
    rep_singlet_df = clonotypes_to_df(rep_singlet)
    rep_singlet_df = rep_singlet_df.merge(rep_df, left_on='root', right_on='index')
    rep_singlet_df['index'] = rep_singlet_df['index_x']
    rep_singlet_df['sequence'] = rep_singlet_df['sequence_x']
    rep_singlet_df['single'] = ['single']*len(rep_singlet_df.index)

    rep_tree_df = clonotypes_to_df(rep_tree)
    rep_tree_df = rep_tree_df.merge(rep_df, left_on='root', right_on='index')
    rep_tree_df['index'] = rep_tree_df['index_x']
    rep_tree_df['sequence'] = rep_tree_df['sequence_x']
    rep_tree_df['single'] = ['tree'] * len(rep_tree_df.index)

    print('rep_singlet')
    print(rep_singlet_df.head())

    print('rep_tree')
    print(rep_tree_df.head())

    rep = merge_repertoires(rep_tree_df, rep_singlet_df, ['index', 'parent', 'root'])
    rep = rep.drop(['sequence_x', 'sequence_y', 'index_x', 'index_y'], axis=1)

    return rep
