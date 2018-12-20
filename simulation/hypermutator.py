import numpy as np
import pandas as pd
import multiprocessing as mp
from copy import deepcopy

class Clonotype:
    """ Clonotype object simulates single B cell clonotype

    Arguments:
    seq = nucleotide sequence (str or np.array)
    seq_aa_fwr = concatenated amino acid sequences of FWRs
    seq_aa_cdr = concatenated amino acid sequences of CDRs
    index = unique identifier
    parent = index of parent Clonotype
    root = index of tree root Clonotype
    proliferation = proliferation probability
    sub_list = list of nucleotide coordinates of all substitutions
    indel_list = list of nucleotide coordinates of all insertions and deletions
    indel_lens = list of lengths of insertions and deletions (corresponds to coordinates in indel_list)
    sub_total = total number of subsitutions
    sub_new = number of substitutions that are not present in parent Clonotype
    indel_total = total number of insertions and deletions
    indel_new = number of indels that are not present in parent Clonotype
    region_ends = np.array of coordinates of FWR1, CDR1, FWR2, CDR2, FWR3, CDR3, FWR4 ends; it is initially obtained from segm_prop.txt and updated after each proliferation event
    segment_ends = list of coordinates of V, D and J segments
    fwr_r = number of replacement substitutions in FWR
    cdr_r = number of replacement substitutions in CDR
    fwr_s = number of silent substitutions in FWR
    cdr_s = number of silent substitutions in CDR

    """

    def __init__(self, seq, seq_aa_fwr, seq_aa_cdr, index, parent, root, proliferation, sub_list, indel_list, indel_lens, sub_total, sub_new, indel_total, indel_new,
                 region_ends, segment_ends, fwr_r, cdr_r, fwr_s, cdr_s):
        self.seq = seq
        self.seq_aa_fwr = seq_aa_fwr
        self.seq_aa_cdr = seq_aa_cdr
        self.index = index
        self.parent = parent
        self.root = root
        self.proliferation = proliferation
        self.sub_total = sub_total
        self.sub_new = sub_new
        self.indel_total = indel_total
        self.indel_new = indel_new
        self.sub_list = sub_list
        self.indel_list = indel_list
        self.indel_lens = indel_lens
        self.region_ends = region_ends
        self.segment_ends = segment_ends
        self.fwr_r = fwr_r
        self.cdr_r = cdr_r
        self.fwr_s = fwr_s
        self.cdr_s = cdr_s

    def update_regions(self, indel_coords_list, indel_len_list):
        # Shift region coordinates if any indels occured

        for i, e in enumerate(indel_coords_list):
            indel_in_reg = 0
            while e > self.region_ends[indel_in_reg]:
                indel_in_reg += 1
                if indel_in_reg == 6:
                    break
            for r in range(indel_in_reg, 6):
                self.region_ends[r] += int(indel_len_list[i])


    def update_segments(self, indel_coords_list, indel_len_list):
        # Shift segment coordinates if any indels occured

        for i, e in enumerate(indel_coords_list):
            indel_in_segm = 0
            while e > self.segment_ends[indel_in_segm]:
                indel_in_segm += 1
                if indel_in_segm == 4:
                    break
            for r in range(indel_in_segm, 4):
                self.segment_ends[r] += int(indel_len_list[i])


def correct_mut_rate(mut_rate, max_tree_height):
    # These coefficients were calculated empirically

    if max_tree_height == 15:
        coeff = 0.1572 - 0.0004
    elif max_tree_height == 5:
        coeff = 0.543
    return mut_rate * coeff


def fwr_mut_rate(mut_rate, ratio=3.2):
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
    rate = 0.007 * mut_rate
    indel_matrix = np.loadtxt('models/hypermutator/model_indels_ant.txt')
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
    # Create one-hot encoding matrix
    matrix = np.zeros((seq.shape[0], 4))
    matrix[np.arange(seq.shape[0]), seq] = 1
    return matrix


def random_choice_prob_index(a, axis=1):
    r = np.expand_dims(np.random.rand(a.shape[1-axis]), axis=axis)
    return (a.cumsum(axis=axis) > r).argmax(axis=axis)


def combine_mat_by_region(matrix_cdr, matrix_fwr, ends):
    # Cut CDR and FWR matrices by region coordinates and concatenate them
    fwr1 = matrix_fwr[range(0, ends[0] + 1)]
    cdr1 = matrix_cdr[range(ends[0] + 1, ends[1] + 1)]
    fwr2 = matrix_fwr[range(ends[1] + 1, ends[2] + 1)]
    cdr2 = matrix_cdr[range(ends[2] + 1, ends[3] + 1)]
    fwr3 = matrix_fwr[range(ends[3] + 1, ends[4] + 1)]
    cdr3 = matrix_cdr[range(ends[4] + 1, ends[5] + 1)]
    fwr4 = matrix_fwr[range(ends[5] + 1, matrix_fwr.shape[0])]

    return np.vstack((fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4))


def translate(array, array_out=True):
    # Translates nucleotide numeric array into amino acid string

    code = {(0, 0, 0): 'K', (0, 0, 1): 'N', (0, 0, 2): 'K', (0, 0, 3): 'N',
            (0, 1, 0): 'T', (0, 1, 1): 'T', (0, 1, 2): 'T', (0, 1, 3): 'T',
            (0, 2, 0): 'R', (0, 2, 1): 'S', (0, 2, 2): 'R', (0, 2, 3): 'S',
            (0, 3, 0): 'I', (0, 3, 1): 'I', (0, 3, 2): 'M', (0, 3, 3): 'I',
            (1, 0, 0): 'Q', (1, 0, 1): 'H', (1, 0, 2): 'Q', (1, 0, 3): 'H',
            (1, 1, 0): 'P', (1, 1, 1): 'P', (1, 1, 2): 'P', (1, 1, 3): 'P',
            (1, 2, 0): 'R', (1, 2, 1): 'R', (1, 2, 2): 'R', (1, 2, 3): 'R',
            (1, 3, 0): 'L', (1, 3, 1): 'L', (1, 3, 2): 'L', (1, 3, 3): 'L',
            (2, 0, 0): 'E', (2, 0, 1): 'D', (2, 0, 2): 'E', (2, 0, 3): 'D',
            (2, 1, 0): 'A', (2, 1, 1): 'A', (2, 1, 2): 'A', (2, 1, 3): 'A',
            (2, 2, 0): 'G', (2, 2, 1): 'G', (2, 2, 2): 'G', (2, 2, 3): 'G',
            (2, 3, 0): 'V', (2, 3, 1): 'V', (2, 3, 2): 'V', (2, 3, 3): 'V',
            (3, 0, 0): '*', (3, 0, 1): 'Y', (3, 0, 2): '*', (3, 0, 3): 'Y',
            (3, 1, 0): 'S', (3, 1, 1): 'S', (3, 1, 2): 'S', (3, 1, 3): 'S',
            (3, 2, 0): '*', (3, 2, 1): 'C', (3, 2, 2): 'W', (3, 2, 3): 'C',
            (3, 3, 0): 'L', (3, 3, 1): 'F', (3, 3, 2): 'L', (3, 3, 3): 'F'}
    aa = sorted(set(code.values()))
    aa_code = {y:x for x, y in enumerate(aa)}
    code2 = {x:aa_code[code[x]] for x in code.keys()}

    codons = np.split(array, range(3, array.shape[0], 3))
    if codons[-1].shape[0] < 3:
        codons.pop(-1)
    if array_out:
        seq_aa = np.array([code2[tuple(x)] for x in codons])
    else:
        seq_aa = ''.join([code[tuple(x)] for x in codons])
    return seq_aa


def split_mat_by_region(array, region_ends):
    # Cut matrix by region coordinates and concatenate CDRs and FWRs separately
    regions = np.split(array, region_ends)
    return (np.hstack( np.array(regions)[[0,2,4,6]] ), np.hstack( np.array(regions)[[1,3,5]] ))


def mutate(cl, mutation_model):
    # Stochastic events of Clonotype proliferation and mutation
    np.random.seed()

    if np.random.random() < cl.proliferation:

        # first introduce substitutions
        subs_probas_cdr = np.dot(matrix_repr(cl.seq), mutation_model['subs_cdr'])
        subs_probas_fwr = np.dot(matrix_repr(cl.seq), mutation_model['subs_fwr'])
        subs_probas = combine_mat_by_region(subs_probas_cdr, subs_probas_fwr, cl.region_ends)

        seq_new = random_choice_prob_index(subs_probas)
        sub_idxs = np.where(seq_new != cl.seq)[0]

        # count R and S
        if sub_idxs.shape[0] > 0:
            try:
                sub_map_fwr, sub_map_cdr = split_mat_by_region(seq_new != cl.seq, cl.region_ends)
            except TypeError:
                print(cl.seq.shape[0])
                print(seq_new.shape[0])
                raise

            sub_fwr, sub_cdr = sub_map_fwr.sum(), sub_map_cdr.sum()
            seq_fwr, seq_cdr = split_mat_by_region(seq_new, cl.region_ends)
            if sub_fwr:
                seq_aa_fwr_tmp = translate(seq_fwr)
                fwr_r = np.where(seq_aa_fwr_tmp != cl.seq_aa_fwr)[0].shape[0]
            else:
                fwr_r = 0

            if sub_cdr:
                seq_aa_cdr_tmp = translate(seq_cdr)
                cdr_r = np.where(seq_aa_cdr_tmp != cl.seq_aa_cdr)[0].shape[0]
            else:
                cdr_r = 0
        else:
            cdr_r, fwr_r, sub_cdr, sub_fwr = 0, 0, 0, 0

        cdr_s, fwr_s = sub_cdr - cdr_r, sub_fwr - fwr_r

        # introduce indels
        indel_len_cdr = np.random.choice(mutation_model['indel_cdr'][0, :], seq_new.shape[0],
                                         p=mutation_model['indel_cdr'][1, :])
        indel_len_cdr = np.transpose(indel_len_cdr[np.newaxis, :])
        indel_len_fwr = np.zeros(indel_len_cdr.shape) # no indels in FWR

        indel_len = combine_mat_by_region(indel_len_cdr, indel_len_fwr, cl.region_ends)
        indel_idxs = np.where(indel_len != 0)[0]
        indel_len = np.transpose(indel_len[indel_idxs])[0,:].astype(int)

        for i in range(indel_idxs.shape[0]):
            idx = indel_idxs[i] + indel_len[:i].sum()
            if indel_len[i] > 0:
                insertion = generate_insertion(indel_len[i])
                seq_new = np.hstack((seq_new[:idx], insertion, seq_new[idx:]))
            else:
                seq_new = np.hstack((seq_new[:idx - 1], seq_new[idx - indel_len[i] - 1:]))

        mut_num = np.hstack((sub_idxs, indel_idxs)).shape[0]

        if mut_num > 0:
            sub_unique = np.unique(np.hstack((cl.sub_list, sub_idxs)))

            cl_new = Clonotype(seq=seq_new, index=0, parent=cl.index, root=cl.root,
                               proliferation = cl.proliferation,
                               sub_total = sub_unique.shape[0],
                               sub_new = sub_idxs.shape[0],
                               indel_total = indel_idxs.shape[0] + cl.indel_total,
                               indel_new = indel_idxs.shape[0],
                               sub_list = sub_unique,
                               indel_list = indel_idxs,
                               indel_lens = indel_len,
                               region_ends = deepcopy(cl.region_ends),
                               segment_ends = deepcopy(cl.segment_ends),
                               seq_aa_fwr = [],
                               seq_aa_cdr = [],
                               fwr_r = cl.fwr_r + fwr_r,
                               cdr_r = cl.cdr_r + cdr_r,
                               fwr_s = cl.fwr_s + fwr_s,
                               cdr_s = cl.cdr_s + cdr_s)

            # update region coords and aa sequences if indels occured
            if indel_idxs.shape[0] > 0:
                cl_new.update_regions(indel_idxs, indel_len)
                cl_new.update_segments(indel_idxs, indel_len)
                cl_new.seq_aa_fwr, cl_new.seq_aa_cdr = (translate(reg) for reg in split_mat_by_region(cl_new.seq, cl_new.region_ends))
            else:
                cl_new.seq_aa_fwr = seq_aa_fwr_tmp if sub_fwr else cl.seq_aa_fwr
                cl_new.seq_aa_cdr = seq_aa_cdr_tmp if sub_cdr else cl.seq_aa_cdr

            # check that each region length > 12 nt
            region_len = np.hstack(([cl_new.region_ends[0]+1], \
                cl_new.region_ends[1:] - cl_new.region_ends[:5], \
                [cl_new.seq.shape[0] - cl_new.region_ends[5]]))

            if sum(region_len >= 12) == 7:
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

        region_ends = np.array([df['cdr1.start'][i] - 1, df['cdr1.end'][i], df['cdr2.start'][i] - 1,
                       df['cdr2.end'][i], df['cdr3.start'][i] - 1, df['cdr3.end'][i]])
        segment_ends = [df['v.end'][i], df['n1.end'][i], df['d.end'][i], df['n2.end'][i]]
        seq = seq_to_array(df['sequence'][i])
        seq_aa_fwr, seq_aa_cdr = (translate(x) for x in split_mat_by_region(seq, region_ends))

        g = Clonotype(seq=seq,
                      seq_aa_fwr = seq_aa_fwr,
                      seq_aa_cdr = seq_aa_cdr,
                      index=df['index'][i],
                      parent=df['index'][i],
                      root=df['index'][i],
                      proliferation=proliferation,
                      sub_total=0,
                      sub_new=0,
                      indel_total=0,
                      indel_new=0,
                      sub_list=[],
                      indel_list=[],
                      indel_lens=[],
                      region_ends = region_ends,
                      segment_ends = segment_ends,
                      fwr_r = 0,
                      cdr_r = 0,
                      fwr_s = 0,
                      cdr_s = 0)
        rep.append(g)

    return rep


def clonotypes_to_df(rep):
    # Create a data frame from a list of Clonotype objects
    return pd.concat([pd.DataFrame([[cl.index, cl.parent, cl.root,
                                     array_to_seq(cl.seq), 
                                     ','.join([str(int(x)) for x in cl.sub_list]),
                                     ','.join([str(int(x)) for x in cl.indel_list]),
                                     ','.join([str(int(x)) for x in cl.indel_lens]),
                                     cl.sub_total, cl.sub_new, cl.indel_total, cl.indel_new, cl.fwr_r, cl.cdr_r, cl.fwr_s, cl.cdr_s] + list(cl.region_ends) + cl.segment_ends],
                                   columns=['index', 'parent', 'root', 'sequence', 'sub_list', 'indel_list', 'indel_len_list',
                                            'sub_total', 'sub_new', 'indel_total', 'indel_new', 'fwr_r', 'cdr_r', 'fwr_s', 'cdr_s',
                                            'fwr1_end', 'cdr1_end', 'fwr2_end', 'cdr2_end',
                                            'fwr3_end', 'cdr3_end', 'v_end','n_end','d_end','n2_end']) for cl in rep], ignore_index=True)


def mutate_repertoire(final_mutation_rate, rep_df):
    # Perform several (max_tree_height) iterations of Clonotype mutation

    rep_df['proliferation'] = np.random.triangular(left=0, mode=0, right=1, size=rep_df.shape[0])
    rep_cl = df_to_clonotypes(rep_df)

    # Load model of SHM and indels in CDRs and FWRs
    cdr_to_fwr_mutations = 3
    max_tree_height = 15
    #mutation_rate = correct_mut_rate(final_mutation_rate, max_tree_height)
    mutation_rate = final_mutation_rate
    fmr = fwr_mut_rate(mutation_rate, cdr_to_fwr_mutations)
    cmr = fmr * cdr_to_fwr_mutations

    mut_model = {'subs_cdr': subs_prob_matrix(cmr), 
                 'subs_fwr': subs_prob_matrix(fmr),
                 'indel_cdr': indels_prob_matrix(cmr),
                 'indel_fwr': indels_prob_matrix(fmr)}

    # Iteratively mutate rep_cl_roots creating novel tree nodes
    for s in range(max_tree_height):
        print(s)
        pool = mp.Pool()
        new_rep = [pool.apply(mutate, args=(cl, mut_model)) for cl in rep_cl]
        new_rep = [cl for cl in new_rep if cl is not None]
        new_rep = [cl for cl in new_rep if cl.fwr_r < 10 & 2 < cl.cdr_r < 15]

        for i, cl in enumerate(new_rep):
            cl.index = len(rep_cl) + 1 + i
        rep_cl = rep_cl + new_rep

    # for s in range(max_tree_height):
    #     print(s)
    #     new_rep = []
    #     for c in rep_cl:
    #         new_rep.append(mutate(c, mut_model))
    #     new_rep = [cl for cl in new_rep if cl is not None]
    #     for i, cl in enumerate(new_rep):
    #         cl.index = len(rep_cl) + 1 + i
    #     rep_cl = rep_cl + new_rep

    rep_df_mutated = clonotypes_to_df(rep_cl)
    rep_df = rep_df.rename(index=str, columns={"index": "root", "V":"v.igor", "J":"j.igor"})
    rep_df = rep_df[['root','v.igor','j.igor','v.name','j.name']]
    rep_df_mutated = rep_df_mutated.merge(rep_df, left_on='root', right_on='root')

    return rep_df_mutated
