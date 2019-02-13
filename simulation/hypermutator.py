import numpy as np
import pandas as pd
import multiprocessing as mp
from copy import deepcopy
from itertools import product


def fivemer_dict():
    nt = [[0,1,2,3]]*5
    fivemer_list = list(product(*nt))
    return {x:i for i,x in enumerate(fivemer_list)}


__fivemer_dict__ = fivemer_dict()


class Sequence:

    def _make_random_end(self):
        nts = [[0,1,2,3]]*2
        dimers = list(product(*nts))
        ind = int(np.random.choice(16, 1))
        return np.array(dimers[ind])

    def _init_kmer_array(self):
        # the first and the last two positions are assigned to random 5mers
        global __fivemer_dict__
        array = np.hstack((self._make_random_end(), self.nt_array, self._make_random_end()))
        fivemer_list = [array[x-2:x+3] for x in range(2, array.shape[0]-2)]
        fivemers_code = [__fivemer_dict__[tuple(x)] for x in fivemer_list]
        self.kmer_array = np.array(fivemers_code)

    def __init__(self, nt_array=None, kmer_array=None, add_kmers=False):
        self.nt_array = nt_array
        if add_kmers:
            self._init_kmer_array()
        else:
            self.kmer_array = kmer_array

    def translate(self, array_out=True):
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

        codons = np.split(self.nt_array, range(3, self.nt_array.shape[0], 3))
        if codons[-1].shape[0] < 3:
            codons.pop(-1)

        if array_out:
            aa = sorted(set(code.values()))
            aa_code = {y:x for x, y in enumerate(aa)}
            code2 = {x:aa_code[code[x]] for x in code.keys()}
            seq_aa = np.array([code2[tuple(x)] for x in codons])
        else:
            seq_aa = ''.join([code[tuple(x)] for x in codons])
        return seq_aa

    def _one_hot_matrix(self, arr_type='nt'):
        if arr_type == 'nt':
            n, array = 4, self.nt_array
        elif arr_type == 'kmer':
            n, array = 1024, self.kmer_array
        matrix = np.zeros((array.shape[0], n))
        matrix[np.arange(array.shape[0]), array] = 1
        return matrix

    def split_by_region(self, region_ends):
        # Cut matrix by region coordinates and concatenate CDRs and FWRs separately
        nt_regions = np.split(self.nt_array, region_ends+1)
        kmer_regions = np.split(self.kmer_array, region_ends+1)
        return (Sequence(nt_array=np.hstack(np.array(nt_regions)[[0,2,4,6]]),
                         kmer_array=np.hstack(np.array(kmer_regions)[[0,2,4,6]])),
                Sequence(nt_array=np.hstack(np.array(nt_regions)[[1,3,5]]),
                         kmer_array=np.hstack(np.array(kmer_regions)[[1,3,5]])))

    def _ranks_to_probs_mat(self, mutability_ranks, subs_matrix, mutation_rate):
        # transform ranks to probabilities proportional to ranks keeping mean(p) = mutation_rate
        coeff = mutability_ranks.shape[0]*mutation_rate/mutability_ranks.sum()
        mutation_probs = mutability_ranks * coeff
        mutation_probs = mutation_probs[..., np.newaxis]
        # calculate probabilities of substitutions
        mutation_probs_mat = mutation_probs * subs_matrix
        # add probabilities of position staying unmutated instead of zeros
        unmut_probs_mat = self._one_hot_matrix('nt') * (1 - mutation_probs)
        mutation_probs_mat += unmut_probs_mat
        return mutation_probs_mat

    @staticmethod
    def _random_choice_prob_index(a, axis=1):
        r = np.expand_dims(np.random.rand(a.shape[1-axis]), axis=axis)
        return (a.cumsum(axis=axis) > r).argmax(axis=axis)

    def _aa_substitution(self, Sequence2):
        aa1, aa2 = self.translate(), Sequence2.translate()
        n_aa_subs = sum(aa1 != aa2)
        return (n_aa_subs, aa2)

    def introduce_subs(self, model, region):
        if model['name'] == 'steele':
            mutability_ranks = model['mutability'][self.nt_array]
            subs_matrix = model['substitution'][self.nt_array,:]
        elif model['name'] == 's5f':
            mutability_ranks = model['mutability'][self.kmer_array]
            subs_matrix = model['substitution'][self.kmer_array,:]
        probs = self._ranks_to_probs_mat(mutability_ranks, subs_matrix, model['subs_rates'][region])
        seq_with_subs = self._random_choice_prob_index(probs)
        n_nt_subs = sum(seq_with_subs != self.nt_array)
        aa_subs = self._aa_substitution(Sequence(nt_array=seq_with_subs))
        return {'nt_array':seq_with_subs, 'n_nt_subs':n_nt_subs, 'n_aa_subs':aa_subs[0], 'aa_array':aa_subs[1]}


class Clonotype:
    """
    seq = Sequence
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
    region_ends = np.array; it is initially obtained from segm_prop.txt and updated after each proliferation event
    segment_ends = list of coordinates of V, D and J segments
    fwr_r = number of replacement substitutions in FWR
    cdr_r = number of replacement substitutions in CDR
    fwr_s = number of silent substitutions in FWR
    cdr_s = number of silent substitutions in CDR

    """

    def __init__(self, seq=None, seq_aa_fwr=None, seq_aa_cdr=None, index=None, parent=None, root=None, generation=None, 
                 proliferation=None, sub_list=None, indel_list=None, indel_lens=None, sub_total=None, sub_new=None, 
                 indel_total=None, indel_new=None, region_ends=None, segment_ends=None, fwr_r=None, cdr_r=None, fwr_s=None, cdr_s=None):
        self.seq = seq
        self.seq_aa_fwr = seq_aa_fwr
        self.seq_aa_cdr = seq_aa_cdr
        self.index = index
        self.parent = parent
        self.root = root
        self.generation = generation
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

    def update_regions(self, indel_list):
        # Shift region coordinates if any indel occured
        indel_coords_list, indel_len_list = indel_list
        for i, e in enumerate(indel_coords_list):
            indel_in_reg = 0
            while e > self.region_ends[indel_in_reg]:
                indel_in_reg += 1
                if indel_in_reg == 6:
                    break
            for r in range(indel_in_reg, 6):
                self.region_ends[r] += int(indel_len_list[i])

    def update_segments(self, indel_list):
        # Shift segment coordinates if any indel occured
        indel_coords_list, indel_len_list = indel_list
        for i, e in enumerate(indel_coords_list):
            indel_in_segm = 0
            while e > self.segment_ends[indel_in_segm]:
                indel_in_segm += 1
                if indel_in_segm == 4:
                    break
            for r in range(indel_in_segm, 4):
                self.segment_ends[r] += int(indel_len_list[i])

    def combine_full_mat_by_region(self, matrix_cdr, matrix_fwr):
        # Combine full length matrices
        fwr = np.split(matrix_fwr, self.region_ends+1)
        cdr = np.split(matrix_cdr, self.region_ends+1)
        return np.hstack((fwr[0], cdr[1], fwr[2], cdr[3], fwr[4], cdr[5], fwr[6]))

    def combine_short_mat_by_region(self, matrix_fwr, matrix_cdr):
        # Combine matrices containing only CDR/FWR parts
        regions = ["cdr1", "fwr2", "cdr2", "fwr3", "cdr3"]
        lens = [self.region_ends[x]-self.region_ends[x-1] for x in range(1,6)]
        reg_len = dict(zip(regions, lens))
        reg_len['fwr1'] = self.region_ends[0]

        fwr1 = matrix_fwr[range(0, reg_len['fwr1'] + 1)]
        cdr1 = matrix_cdr[range(0, reg_len['cdr1'])]
        fwr2 = matrix_fwr[range(reg_len['fwr1'] + 1, reg_len['fwr1'] + reg_len['fwr2'] + 1)]
        cdr2 = matrix_cdr[range(reg_len['cdr1'], reg_len['cdr1'] + reg_len['cdr2'])]
        fwr3 = matrix_fwr[range(reg_len['fwr1'] + reg_len['fwr2'] + 1, reg_len['fwr1'] + reg_len['fwr2'] + reg_len['fwr3'] + 1)]
        cdr3 = matrix_cdr[range(reg_len['cdr1'] + reg_len['cdr2'], matrix_cdr.shape[0])]
        fwr4 = matrix_fwr[range(reg_len['fwr1'] + reg_len['fwr2'] + reg_len['fwr3'] + 1, matrix_fwr.shape[0])]
        return np.hstack((fwr1, cdr1, fwr2, cdr2, fwr3, cdr3, fwr4))

    def split_mat_by_region(self, array):
        # Cut matrix by region coordinates and concatenate CDRs and FWRs separately
        regions = np.split(array, self.region_ends+1)
        return (np.hstack( np.array(regions)[[0,2,4,6]] ), np.hstack( np.array(regions)[[1,3,5]] ))

    def substitution(self, fwr, cdr, model):
        fwr_subs_res = fwr.introduce_subs(model, 'fwr')
        cdr_subs_res = cdr.introduce_subs(model, 'cdr')
        try:
            seq_with_subs = self.combine_short_mat_by_region(fwr_subs_res['nt_array'], cdr_subs_res['nt_array'])
        except TypeError:
            print(probs_by_region[0].shape[0])
            print(probs_by_region[1].shape[0])
            raise
        return {'full_seq':seq_with_subs, 'fwr_result':fwr_subs_res, 'cdr_result':cdr_subs_res}

    def make_indel_list(self, fwr, cdr, model):
        indel_len_cdr = np.random.choice(model['indel']['cdr'][0, :], cdr.nt_array.shape[0],
                                             p=model['indel']['cdr'][1, :])
        indel_len_fwr = np.random.choice(model['indel']['fwr'][0, :], fwr.nt_array.shape[0],
                                             p=model['indel']['fwr'][1, :])
        indel_len = self.combine_short_mat_by_region(indel_len_fwr, indel_len_cdr)
        indel_idxs = np.where(indel_len != 0)[0]
        indel_len = indel_len[indel_idxs].astype(int)
        return (indel_idxs, indel_len)

    def generate_insertion(self, n):
        np.random.seed()
        return np.random.choice(4, n)

    def introduce_indels(self, seq, indel_list):
        indel_idxs, indel_len = indel_list
        for i in range(indel_idxs.shape[0]):
                idx = indel_idxs[i] + indel_len[:i].sum()
                if indel_len[i] > 0:
                    insertion = self.generate_insertion(indel_len[i])
                    seq = np.hstack((seq[:idx], insertion, seq[idx:]))
                else:
                    seq = np.hstack((seq[:idx - 1], seq[idx - indel_len[i] - 1:]))
        return seq

    def make_mutant(self, model):
        np.random.seed()
        fwr, cdr = self.seq.split_by_region(self.region_ends)

        substitution = self.substitution(fwr, cdr, model)
        n_subs = substitution['cdr_result']['n_nt_subs'] + substitution['fwr_result']['n_nt_subs']
        indel_list = self.make_indel_list(fwr, cdr, model)

        if n_subs + indel_list[0].shape[0] > 0:
            if n_subs > 0:
                sub_idxs = np.where(substitution['full_seq'] != self.seq.nt_array)[0]
                sub_unique = np.unique(np.hstack((self.sub_list, sub_idxs)))
            else:
                sub_idxs, sub_unique = np.array([]), np.array([])

            cl_new = Clonotype(index=0, parent=self.index, root=self.root,
                               generation = self.generation + 1,
                               proliferation = self.proliferation,
                               sub_total = sub_unique.shape[0],
                               sub_new = sub_idxs.shape[0],
                               sub_list = sub_unique,
                               region_ends = deepcopy(self.region_ends),
                               segment_ends = deepcopy(self.segment_ends),
                               fwr_r = self.fwr_r + substitution['fwr_result']['n_aa_subs'],
                               cdr_r = self.cdr_r + substitution['cdr_result']['n_aa_subs'],
                               fwr_s = self.fwr_s + substitution['fwr_result']['n_nt_subs'] - substitution['fwr_result']['n_aa_subs'],
                               cdr_s = self.cdr_s + substitution['cdr_result']['n_nt_subs'] - substitution['cdr_result']['n_aa_subs'])

            if indel_list[0].shape[0] > 0:
                new_seq = self.introduce_indels(substitution['full_seq'], indel_list)
                cl_new.seq = Sequence(nt_array=new_seq, add_kmers=True)
                cl_new.indel_total = indel_list[0].shape[0] + self.indel_total
                cl_new.indel_new = indel_list[0].shape[0]
                cl_new.indel_list, cl_new.indel_lens = indel_list
                cl_new.update_regions(indel_list)
                cl_new.update_segments(indel_list)
                cl_new.seq_aa_fwr, cl_new.seq_aa_cdr = (reg.translate() for reg in cl_new.seq.split_by_region(cl_new.region_ends))
            else:
                cl_new.seq = Sequence(nt_array=substitution['full_seq'], add_kmers=True)
                cl_new.indel_total = self.indel_total
                cl_new.indel_new = 0
                cl_new.indel_list, cl_new.indel_lens = (np.array([]),np.array([]))
                cl_new.seq_aa_fwr = substitution['fwr_result']['aa_array']
                cl_new.seq_aa_cdr = substitution['cdr_result']['aa_array']

            # check that each region length > 12 nt
            region_len = np.hstack(([cl_new.region_ends[0]+1], \
                cl_new.region_ends[1:] - cl_new.region_ends[:5], \
                [cl_new.seq.nt_array.shape[0] - cl_new.region_ends[5]]))
            if sum(region_len >= 12) == 7:
                return cl_new


def proliferate(Clonotype, model):
    if np.random.random() < Clonotype.proliferation:
        return Clonotype.make_mutant(model)


def correct_mut_rate(mut_rate, max_tree_height):
    # OBSOLETE FUNCTION
    # These coefficients were calculated empirically
    if max_tree_height == 15:
        coeff = 0.1572 - 0.0004
    elif max_tree_height == 5:
        coeff = 0.543
    return mut_rate * coeff


def fwr_cdr_mut_rates(mut_rate, ratio=3.2):
    # Get CDR and FWR mutation rate from the average mutation rate
    cdr_len_prop = 0.15
    fwr_len_prop = 1 - cdr_len_prop
    fwr_mut_rate = mut_rate / (ratio * cdr_len_prop + fwr_len_prop)
    return (fwr_mut_rate, fwr_mut_rate * ratio)


def load_subs_model(model_name):
    model = {'mutability': np.loadtxt('models/model_mutability_' + model_name + '.txt'),
             'substitution': np.loadtxt('models/model_substitution_' + model_name + '.csv', delimiter='\t')}
    return model


def load_indel_model(indel_rate):
    indel_matrix = np.loadtxt('models/model_indels_ant.txt')
    indel_matrix[1, :] *= indel_rate
    indel_matrix = np.c_[indel_matrix, np.array([0, 1 - indel_rate])]  # probability of making no indel
    return indel_matrix


def make_model(model_name, mutation_rate):
    model = {}
    mut_rates = fwr_cdr_mut_rates(mutation_rate)
    subs_rates = [x*0.992 for x in mut_rates]
    indel_rates = [x*0.008 for x in mut_rates]
    model['subs_rates'] = dict(zip(['fwr', 'cdr'], subs_rates))
    model.update(load_subs_model(model_name))
    model['indel'] = {'fwr':load_indel_model(indel_rates[0]),'cdr':load_indel_model(indel_rates[1])}
    model['name'] = model_name
    return model


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


def df_to_clonotypes(df):
    # Create a list of Clonotype objects from a data frame
    rep = []

    for i in range(len(df.index)):
        region_ends = np.array([df['cdr1.start'][i] - 1, df['cdr1.end'][i], df['cdr2.start'][i] - 1,
                       df['cdr2.end'][i], df['cdr3.start'][i] - 1, df['cdr3.end'][i]])
        segment_ends = [df['v.end'][i], df['n1.end'][i], df['d.end'][i], df['n2.end'][i]]
        seq = Sequence(seq_to_array(df['sequence'][i]), add_kmers=True)
        seq_aa_fwr, seq_aa_cdr = (reg.translate() for reg in seq.split_by_region(region_ends))

        g = Clonotype(seq=seq,
                      seq_aa_fwr=seq_aa_fwr,
                      seq_aa_cdr=seq_aa_cdr,
                      index=df['index'][i],
                      parent=df['index'][i],
                      root=df['index'][i],
                      generation=0,
                      proliferation=df['proliferation'][i],
                      sub_total=0,
                      sub_new=0,
                      indel_total=0,
                      indel_new=0,
                      sub_list=[],
                      indel_list=[],
                      indel_lens=[],
                      region_ends=region_ends,
                      segment_ends=segment_ends,
                      fwr_r=0,
                      cdr_r=0,
                      fwr_s=0,
                      cdr_s=0)
        rep.append(g)

    return rep


def clonotypes_to_df(rep):
    # Create a data frame from a list of Clonotype objects
    return pd.concat([pd.DataFrame([[cl.index,
                                     cl.parent,
                                     cl.root,
                                     cl.generation,
                                     array_to_seq(cl.seq.nt_array),
                                     ','.join([str(int(x)) for x in cl.sub_list]),
                                     ','.join([str(int(x)) for x in cl.indel_list]),
                                     ','.join([str(int(x)) for x in cl.indel_lens]),
                                     cl.sub_total,
                                     cl.sub_new,
                                     cl.indel_total,
                                     cl.indel_new,
                                     cl.fwr_r,
                                     cl.cdr_r,
                                     cl.fwr_s,
                                     cl.cdr_s] +
                                    list(cl.region_ends) +
                                    cl.segment_ends],
                                   columns=['index', 'parent', 'root', 'generation', 'sequence',
                                            'sub_list', 'indel_list', 'indel_len_list',
                                            'sub_total', 'sub_new', 'indel_total', 'indel_new',
                                            'fwr_r', 'cdr_r', 'fwr_s', 'cdr_s',
                                            'fwr1_end', 'cdr1_end', 'fwr2_end', 'cdr2_end',
                                            'fwr3_end', 'cdr3_end', 'v_end','n_end','d_end','n2_end']) for cl in rep],
                     ignore_index=True)


def mutate_repertoire(final_mutation_rate, model_name, rep_df):
    # Perform several (max_tree_height) iterations of Clonotype mutation

    #p = np.random.triangular(left=0, mode=0, right=1, size=rep_df.shape[0])
    p = np.random.exponential(size=rep_df.shape[0])
    #p = np.array(range(len(rep_df.index)))
    #p = np.exp(np.exp(-p**5))*(p**5)
    p = p/max(p)
    rep_df['proliferation'] = p
    rep_cl = df_to_clonotypes(rep_df)

    max_tree_height = 10
    #mutation_rate = correct_mut_rate(final_mutation_rate, max_tree_height)
    mutation_rate = final_mutation_rate
    mut_model = make_model(model_name, mutation_rate)

    # Iteratively mutate rep_cl creating novel tree nodes
    for s in range(max_tree_height):
        print(s)
        pool = mp.Pool()
        new_rep = [pool.apply(proliferate, args=(cl, mut_model)) for cl in rep_cl]
        new_rep = [cl for cl in new_rep if cl is not None]
        new_rep = [cl for cl in new_rep if cl.fwr_r < 8]
        new_rep = [cl for cl in new_rep if ((cl.cdr_r < 2) or (cl.cdr_s == 0) or (cl.cdr_r/cl.cdr_s > 2))]

        for i, cl in enumerate(new_rep):
            cl.index = len(rep_cl) + 1 + i
        rep_cl = rep_cl + new_rep

    # for s in range(max_tree_height):
    #     print(s)
    #     new_rep = []
    #     for c in rep_cl:
    #         new_rep.append(c.make_mutant(mut_model))
    #     new_rep = [cl for cl in new_rep if cl is not None]
    #     for i, cl in enumerate(new_rep):
    #         cl.index = len(rep_cl) + 1 + i
    #     rep_cl = rep_cl + new_rep

    rep_df_mutated = clonotypes_to_df(rep_cl)
    rep_df = rep_df.rename(index=str, columns={"index": "root", "V":"v.igor", "J":"j.igor"})
    rep_df = rep_df[['root','v.igor','j.igor','v.name','j.name']]
    rep_df_mutated = rep_df_mutated.merge(rep_df, left_on='root', right_on='root')

    return rep_df_mutated