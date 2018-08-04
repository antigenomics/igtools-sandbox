import argparse
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE, call
from copy import deepcopy
from scipy.stats import gamma
from hypermutator import *
#from sequencer import *


class IgorError(Exception):
    pass


def run_igor(size, model, dir, batch):
    # Note: do not forget to add correct path to igor to PATH variable

    if model.startswith('bcr'):
        chain = 'heavy_naive'
    else:
        chain = model[4:]

    cmd = ["igor",
           "-generate", str(size),
           "--noerr",
           "-set_wd", dir,
           "-chain", chain,
           "-batch", batch,
           "-set_custom_model", "/home/anna/simulation/models/igor/" + model + "_" + batch + "/model_parms.txt",
           "/home/anna/simulation/models/igor/" + model + "_" + batch + "/model_marginals.txt",
           "-species", "human"]

    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    with open(dir + 'igor_log.txt', 'a') as f:
            f.write('COMMAND:\n' + ' '.join(cmd) + '\n\n')
            f.write('STDOUT:\n' + stdout.decode("utf-8") + '\n')
            f.write('STDERR:\n' + stderr.decode("utf-8") + '\n')
            f.write('\n')

    if 'cannot' in stderr.decode("utf-8"):
        raise IgorError('IGoR failed, look igor_log.txt for more information')


def read_igor_output(dir, model):
    s = pd.read_csv(dir + 'generated_seqs_noerr.csv', sep=';')
    r = pd.read_csv(dir + 'generated_realizations_noerr.csv', sep=';')
    if model.startswith('bcr'):
        v_col = 'GeneChoice_V_gene_Undefined_side_prio8_size97'
        j_col = 'GeneChoice_J_gene_Undefined_side_prio7_size7'
    else:
        v_col = 'GeneChoice_V_gene_Undefined_side_prio7_size103'
        j_col = 'GeneChoice_J_gene_Undefined_side_prio6_size68'

    r['V'] = r[v_col].apply(lambda x: x.strip('()'))
    r['J'] = r[j_col].apply(lambda x: x.strip('()'))
    s['sequence'] = s['nt_sequence']
    s['index'] = s['seq_index']

    return pd.concat([s[['index', 'sequence']], r[['V', 'J']]], axis=1)


def set_counts(df):
    # the proportion of counts > 140 is about 0.5 if scale == 200
    threshold = 140

    df['clone_count'] = np.random.exponential(scale=200, size=len(df.index))
    df['clone_count'] = df['clone_count'].apply(lambda x: round(x) - threshold)
    df = df.loc[df['clone_count'] > 0]

    return df


def make_rep(model, size, dir, mutation_rate=0.05, allelic_rate=0.5):
    if model.startswith('bcr'):
        # Here BCR repertoire is made up of 0.5 singletons, 0.5 * 1/5 tree roots and 0.5 * 4/5 tree nodes
        # Then the repertoire will be downsampled to just a half
        # Therefore, we need to generate n = 2 * size * 6/10 initial sequences
        n = int(round(size * 0.6 * 2))
    else:
        # TCR repertoire is just halved after downsampling
        n = size * 2

    # Generate non-allelic initial repertoire
    try:
        run_igor(size=int(round(n * (1 - allelic_rate))),
                 model=model,
                 dir=dir,
                 batch='default')
    except IgorError:
        raise

    rep_default = read_igor_output(dir + 'default_generated/', model)
    rep_default['descr'] = ['default'] * len(rep_default.index)
    print('default')
    print(len(rep_default.index))

    # Generate allelic initial repertoire
    try:
        run_igor(size = int(round(n * allelic_rate)),
                 model = model,
                 dir = dir,
                 batch = 'allelic')
    except IgorError:
        raise

    rep_allelic = read_igor_output(dir + 'allelic_generated/', model)
    rep_allelic['descr'] = ['allelic'] * len(rep_allelic.index)
    print('allelic')
    print(len(rep_allelic.index))

    # Copying the second df to ignore corrections made during merging 
    rep = merge_repertoires(rep_default, deepcopy(rep_allelic), ['index'])
    print('after merging')
    print(len(rep.index))

    # Introduce SHMs and simulate cell proliferation
    if args.model.startswith('bcr'):
        rep = mutate_repertoire(mutation_rate, rep)

    # Add random counts and downsample
    rep = set_counts(rep)
    rep = rep.reset_index(drop=True)

    return rep


def write_fasta(df):
    lines = []

    for i in range(len(df.index)):
        for j in range(df['clone_count'][i]):
            lines.append('>' + df['descr'][i] + '_clone_' + str(df['index'][i]) + '_' + str(j) + '\n' + df['sequence'][i] + '\n')

    with open(args.dir+'clones.fasta', 'w') as f:
        f.write(''.join(lines))


def run_sequencer(dir, sequencing_mode):
    call(['mkdir', dir + 'reads'])

    if sequencing_mode == 'RNASeq':
        l = 75
        ss = 'NS50'
    elif sequencing_mode == 'MiSeq':
        l = 250
        ss = 'MSv1'
    elif sequencing_mode == 'HiSeq':
        l = 150
        ss = 'HSXn'

    call(['art_illumina',
          '-p',
          '-i', dir+'clones.fasta',
          '-l', str(l),
          '-ss', ss,
          '-f','10',
          '-m','1000',
          '-s','0',
          '-qL', '32',
          '-o', dir+'reads/'])


def run_mixcr(dir):
    call(["mkdir", dir+"mixcr"])

    call(["mixcr", "align", "-a", "-g", "-t", "50", dir+"reads/1.fq", dir+"reads/2.fq", dir+"mixcr/alignments.vdjca"])
    call(["mixcr", "assemble",
          "-t", "50",
          '-OassemblingFeatures="[FR1+CDR1+FR2+CDR2+FR3+CDR3+FR4]"',
          "--index", dir+"mixcr/index",
          dir+"mixcr/alignments.vdjca",
          dir+"mixcr/clones.clns"])
    call(["mixcr", "exportAlignments", "-p", "full", "-cloneId",
          dir+"mixcr/index", dir + "mixcr/alignments.vdjca",
          dir + "mixcr/alignments.txt"])
    call(["mixcr", "exportClones",
          "--preset", "full",
          "-readIds", dir+"mixcr/index",
          dir+"mixcr/clones.clns", dir+"mixcr/clones.txt"])


parser = argparse.ArgumentParser(description="Simulate B or T cell receptor repertoire")
parser.add_argument("-m", dest="model", type=str, help='Name of the model (bcr_heavy/tcr_alpha)', default='bcr_heavy')
parser.add_argument("-a", dest="allelic_rate", type=float, help='Share of allelic variants (0.5 for heterozygosity)', default=0.5)
parser.add_argument("-n", dest="size", type=int, help='Approximate size of the repertoire (minimal)')
parser.add_argument("-d", dest="dir", type=str, help='Output directory')
parser.add_argument("-mr", dest="mutation_rate", type=float, help='Bulk mutation rate')
parser.add_argument("-s", dest="sequencing_mode", type=str, help='Type of sequencing data (MiSeq/HiSeq/RNASeq')
args = parser.parse_args()

rep = make_rep(args.model, args.size, args.dir, args.mutation_rate)
print('Repertoire is generated')
rep.to_csv(args.dir + 'clones_init.csv', sep='\t', index=False)
write_fasta(rep)

run_sequencer(args.dir, args.sequencing_mode)
print('Sequencing is simulated')

run_mixcr(args.dir)
print('MiXCR processing is done')
