import random
import numpy as np
import pandas as pd
import multiprocessing as mp


class Clonotype:

    def __init__(self, seq, index, parent, root, proliferation, sub_total, sub_new, indel_total, indel_new, sub_list, region_ends):
        self.seq = seq
        self.index = index
        self.parent = parent
        self.root = root
        self.proliferation = proliferation
        self.sub_total = sub_total
        self.sub_new = sub_new
        self.indel_total = indel_total
        self.indel_new = indel_new
        self.sub_list = sub_list
        self.region_ends = region_ends

    def update_regions(self, coords_list, len_list):
        for i, e in enumerate(coords_list):
            indel_in_reg = 0
            while e > self.region_ends[indel_in_reg]:
                indel_in_reg += 1
                if indel_in_reg == 6:
                    break
            for r in range(indel_in_reg, 6):
                self.region_ends[r] += int(len_list[i])

def merge_repertoires(df1, df2, column_list):
    # Reset indices in the listed columns in the second df
    df1_max = df1[column_list].max().max()
    df2_min = df2[column_list].min().min()
    for c in column_list:
        df2[c] = df2[c].apply(lambda x: x - df2_min + df1_max + 1)
    df = df1.append(df2)
    df = df.reset_index(drop=True)

    return df


def correct_mut_rate(mut_rate, max_tree_height):
    # These coefficients were calculated empirically
    if max_tree_height == 15:
        coeff = 0.1572 - 0.0004
    elif max_tree_height == 5:
        coeff = 0.543
    return mut_rate * coeff


def fwr_mut_rate(mut_rate, ratio=3):
    # CDR mutation rate = ratio * FWR mutation rate
     cdr_len_prop = 0.15
     fwr_len_prop = 1 - cdr_len_prop

     fwr_mut_rate = mut_rate / (ratio * cdr_len_prop + fwr_len_prop)
     return fwr_mut_rate


def subs_prob_matrix(mut_rate):
    # The assessed frequency of substitutions among all mutations is 0.992
    rate = 0.992 * mut_rate
    # rate = mut_rate
    f = open('models/hypermutator/model_subs_steele.txt')
    substitution_matrix = np.zeros((4, 4))
    for i, line in enumerate(f):
        row = [float(x) * rate * 4 for x in line.split()]
        row = row[:i] + [1 - sum(row)] + row[i:]
        substitution_matrix[i, :] = np.array(row)
    return substitution_matrix


def indels_prob_matrix(mut_rate):
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


def combine_mat_by_region(matrix_cdr, matrix_fwr, ends):
    fwr1 = matrix_fwr[range(0, ends[0] + 1)]
    cdr1 = matrix_cdr[range(ends[0] + 1, ends[1] + 1)]
    fwr2 = matrix_fwr[range(ends[1] + 1, ends[2] + 1)]
    cdr2 = matrix_cdr[range(ends[2] + 1, ends[3] + 1)]
    fwr3 = matrix_fwr[range(ends[3] + 1, ends[4] + 1)]
    cdr3 = matrix_cdr[range(ends[4] + 1, ends[5] + 1)]
    fwr4 = matrix_fwr[range(ends[5] + 1, matrix_fwr.shape[0])]

    return np.vstack((fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4))


def mutate(cl, mutation_model):
    np.random.seed()

    if np.random.random() < cl.proliferation:

        # introduce substitutions firstly
        subs_probas_cdr = np.dot(matrix_repr(cl.seq), mutation_model['subs_cdr'])
        subs_probas_fwr = np.dot(matrix_repr(cl.seq), mutation_model['subs_fwr'])
        subs_probas = combine_mat_by_region(subs_probas_cdr, subs_probas_fwr, cl.region_ends)

        new_seq = random_choice_prob_index(subs_probas)
        sub_idxs = np.where(new_seq != cl.seq)[0]

        # then introduce indels
        indel_len_cdr = np.random.choice(mutation_model['indel_cdr'][0, :], new_seq.shape[0],
                                         p=mutation_model['indel_cdr'][1, :])
        indel_len_fwr = np.random.choice(mutation_model['indel_fwr'][0, :], new_seq.shape[0],
                                         p=mutation_model['indel_fwr'][1, :])
        indel_len_cdr = np.transpose(indel_len_cdr[np.newaxis, :])
        indel_len_fwr = np.transpose(indel_len_fwr[np.newaxis, :])

        indel_len = combine_mat_by_region(indel_len_cdr, indel_len_fwr, cl.region_ends)
        indel_idxs = np.where(indel_len != 0)[0]
        indel_len = np.transpose(indel_len[indel_idxs])[0,:].astype(int)

        for i in range(indel_idxs.shape[0]):
            idx = indel_idxs[i] + indel_len[:i].sum()
            if indel_len[i] > 0:
                insertion = generate_insertion(indel_len[i])
                new_seq = np.hstack((new_seq[:idx], insertion, new_seq[idx:]))
            else:
                new_seq = np.hstack((new_seq[:idx - 1], new_seq[idx - indel_len[i] - 1:]))

        mut_idxs = np.hstack((sub_idxs, indel_idxs))
        mut_num = mut_idxs.shape[0]
        sub_unique = np.unique(np.hstack((cl.sub_list, sub_idxs)))

        if mut_num > 0:
            cl_new = Clonotype(seq=new_seq, index=0, parent=cl.index, root=cl.root,
                             proliferation = cl.proliferation,
                             sub_total = sub_unique.shape[0], 
                             sub_new = sub_idxs.shape[0],
                             indel_total = indel_idxs.shape[0] + cl.indel_total, 
                             indel_new = indel_idxs.shape[0],
                             sub_list = sub_unique,
                             region_ends = cl.region_ends)
            cl_new.update_regions(indel_idxs, indel_len)
            return cl_new


def df_to_clonotypes(df):
    # Create a list of Clonotype objects from a data frame
    rep = []

    for i in range(len(df.index)):
        # If no proliferation probabilities provided, assign initial proliferation probabilities equal to 1,
        # because all clonotypes must mutate firstly to simulate "old" mutations
        if 'proliferation' in df:
            proliferation = df['proliferation'][i]
        else:
            proliferation = 1

        g = Clonotype(seq=seq_to_array( df['sequence'][i] ),
                      index=df['index'][i],
                      parent=df['index'][i],
                      root=df['index'][i],
                      proliferation=proliferation,
                      sub_total=0,
                      sub_new=0,
                      indel_total=0,
                      indel_new=0,
                      sub_list=[],
                      region_ends = [df['cdr1.start'][i] - 1,
                                     df['cdr1.end'][i],
                                     df['cdr2.start'][i] - 1,
                                     df['cdr2.end'][i],
                                     df['cdr3.start'][i] - 1,
                                     df['cdr3.end'][i]])
        rep.append(g)

    return rep


def clonotypes_to_df(rep):
    # Create a data frame from a list of Clonotype objects
    return pd.concat([pd.DataFrame([[cl.index,
                                     cl.parent,
                                     cl.root,
                                     ','.join([str(x) for x in cl.sub_list]),
                                     array_to_seq(cl.seq),
                                     cl.sub_total, cl.sub_new, cl.indel_total, cl.indel_new] + cl.region_ends],
                                   columns=['index', 'parent', 'root', 'sub_list', 'sequence',
                                            'sub_total', 'sub_new', 'indel_total', 'indel_new',
                                            'fwr1_ends', 'cdr1_ends', 'fwr2_ends', 'cdr2_ends',
                                            'fwr3_ends', 'cdr3_ends']) for cl in rep], ignore_index=True)


def mutate_roots(rep, mutation_model):
    # Introducing "old" mutations into sequences of tree roots and singletons
    # Does not create new Clonotypes
    new_rep = []

    for g in rep:
        cl = mutate(g, mutation_model)
        if cl is None:
            cl = g
        cl.index = g.index
        new_rep.append(cl)
        del cl

    return new_rep


def mutate_repertoire(final_mutation_rate, rep_df):
    rep_cl = df_to_clonotypes(rep_df)

    # Reassign proliferation probabilities that were assigned to 1 by default
    probs = np.random.triangular(left=0, mode=0, right=1, size=len(rep_cl))
    for i, cl in enumerate(rep_cl):
        cl.proliferation = probs[i]

    # Load model of SHM and indels in CDRs and FWRs
    cdr_to_fwr_mutations = 3
    # max_tree_height = 5
    max_tree_height = 15
    mutation_rate = correct_mut_rate(final_mutation_rate, max_tree_height)
    fmr = fwr_mut_rate(mutation_rate, cdr_to_fwr_mutations)
    cmr = fmr * cdr_to_fwr_mutations

    mut_model = {'subs_cdr': subs_prob_matrix(cmr),
                 'subs_fwr': subs_prob_matrix(fmr),
                 'indel_cdr': indels_prob_matrix(cmr),
                 'indel_fwr': indels_prob_matrix(fmr)}

    # Iteratively mutate rep_cl_roots creating novel tree nodes
    for s in range(max_tree_height):
        print(s)
        pool = mp.Pool(2)
        new_rep = [pool.apply(mutate, args=(cl, mut_model)) for cl in rep_cl]
        new_rep = [cl for cl in new_rep if cl is not None]
        for i, cl in enumerate(new_rep):
            cl.index = len(rep_cl) + 1 + i
        rep_cl = rep_cl + new_rep
    print(len(rep_cl))
    print(rep_cl[0])

    rep_df_mutated = clonotypes_to_df(rep_cl)
    rep_df = rep_df.rename(index=str, columns={"index": "root"})
    rep_df = rep_df[['root','V','J','v.name','j.name']]
    rep_df_mutated = rep_df_mutated.merge(rep_df, left_on='root', right_on='root')

    print('rep_df_mutated')
    print(rep_df_mutated.head())

    return rep_df_mutated
