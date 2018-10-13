import argparse
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE, call
from hypermutator import *


class IgorError(Exception):
    pass

def run_igor(size, model, dir):
    # Note: do not forget to add correct path to igor to the PATH variable

    if model.startswith('bcr'):
        chain = 'heavy_naive'
    else:
        chain = model[4:]

    cmd = ["igor",
           "-generate", str(size),
           "--noerr",
           "-set_wd", dir,
           "-chain", chain,
           "-species", "human",
           "-set_custom_model", dir+"model_parms.txt", dir+"model_marginals.txt"
           ]

    process = Popen(cmd, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()

    with open(dir + 'igor_log.txt', 'a') as f:
            f.write('COMMAND:\n' + ' '.join(cmd) + '\n\n')
            f.write('STDOUT:\n' + stdout.decode("utf-8") + '\n')
            f.write('STDERR:\n' + stderr.decode("utf-8") + '\n')
            f.write('\n')

    if 'cannot' in stderr.decode("utf-8"):
        raise IgorError('IGoR failed, look igor_log.txt for more information')


def find_root_number(size, mutation_rate):
    # empirical formula for tree_height = 15
    if mutation_rate < 0.08:
        return int(size / (15240 * mutation_rate - 151.8))
    else:
        return size / 1000


def read_igor_output(dir):
    s = pd.read_csv(dir + 'generated_seqs_noerr.csv', sep=';')
    r = pd.read_csv(dir + 'generated_realizations_noerr.csv', sep=';')

    r['V'] = r[r.columns[1]].apply(lambda x: int(x.strip('()')))
    r['J'] = r[r.columns[2]].apply(lambda x: int(x.strip('()')))
    s['sequence'] = s['nt_sequence']
    s['index'] = s['seq_index']

    df = pd.concat([s[['index', 'sequence']], r[['V', 'J']]], axis=1)
    df = add_region_coords(df)
    return df


def add_region_coords(df):
    props = pd.read_csv('segment_prop.csv', sep='\t')

    props_v = props[props['segment'] == 'Variable']
    props_v = props_v.drop(['sequence', 'gene', 'segment'], axis=1)
    props_v = props_v.rename(index=str, columns={"id": "v.name", "reference_point": "cdr3.start"})

    props_j = props[props['segment'] == 'Joining']
    props_j['cdr3.end.reverse'] = props_j['sequence'].str.len() - props_j['reference_point']
    props_j = props_j.drop(['sequence', 'gene', 'segment', 'reference_point', 'cdr1.start', 'cdr1.end', 'cdr2.start', 'cdr2.end'], axis=1)
    props_j = props_j.rename(index=str, columns={"id": "j.name"})

    df = df.merge(props_v, left_on = 'V', right_on='igor_id')
    df = df.drop(['igor_id'], axis=1)
    df = df.merge(props_j, left_on='J', right_on='igor_id')
    df = df.drop(['igor_id'], axis=1)
    df['cdr3.end'] = df['sequence'].str.len() - df['cdr3.end.reverse']
    return df


def set_counts(df, sampling_rate):
    n = len(df.index)
    counts = np.random.multinomial(n=int(sampling_rate * n), pvals=[1./n]*n, size=1)
    df['clone_count'] = pd.Series(counts[0, :])
    df = df.loc[df['clone_count'] > 0]
    return df


def make_rep(model, size, dir, mr):
    mutation_rate=0.05 if mr == None else mr
    smpl_rate = 0.01
    prop_after_sampling = 0.092 * (smpl_rate ** 2) + smpl_rate * 0.928 + 0.005

    if model.startswith('bcr'):
        n = find_root_number(size, mutation_rate) * (1 / prop_after_sampling)
    else:
        n = size * (1 / prop_after_sampling)

    # Generate initial repertoire
    try:
        run_igor(size=n,
                 model=model,
                 dir=dir)
    except IgorError:
        raise

    rep = read_igor_output(dir + 'generated/')

    # Introduce SHMs and simulate cell proliferation
    if args.model.startswith('bcr'):
        rep = mutate_repertoire(mutation_rate, rep)

    # Add random counts and downsample
    rep.to_csv(args.dir + 'clones_before_sampling.csv', sep='\t')
    rep = set_counts(rep, smpl_rate)
    rep = rep.reset_index(drop=True)

    return rep


def write_fasta(df):
    lines = []

    for i in range(len(df.index)):
        for j in range(df['clone_count'][i]):
            lines.append('>' + 'clone_' + str(df['index'][i]) + '_' + str(j) + '\n' + df['sequence'][i] + '\n')

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
          "-OcloneClusteringParameters=null",
          "-OmaxBadPointsPercent=0",
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
parser.add_argument("-n", dest="size", type=int, help='Approximate size of the repertoire (minimal)')
parser.add_argument("-d", dest="dir", type=str, help='Output directory')
parser.add_argument("-mr", dest="mutation_rate", type=float, help='Bulk mutation rate')
parser.add_argument("-s", dest="sequencing_mode", type=str, help='Type of sequencing data (MiSeq/HiSeq/RNASeq')
args = parser.parse_args()

rep = make_rep(args.model, args.size, args.dir, args.mutation_rate)
print('Repertoire is generated')
rep.to_csv(args.dir + 'clones_init.csv', sep='\t')
write_fasta(rep)

run_sequencer(args.dir, args.sequencing_mode)
print('Sequencing is simulated')

run_mixcr(args.dir)
print('MiXCR processing is done')
