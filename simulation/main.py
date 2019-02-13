import argparse
import os
import numpy as np
import pandas as pd
from subprocess import Popen, PIPE, call
from Bio.Seq import Seq
from hypermutator import *


class IgorError(Exception):
    pass


def run_igor(size, model, dir):
    # do not forget to add correct path to igor to the PATH variable
    if os.path.isdir(dir + "generated"):
        call(["rm", "-r", dir+"generated"])

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


def read_parms_file(dir):
    # Extract feature length from IGoR model to further get segment coordinates
    # Returns a dictionary {feature_name:{id:length}}

    with open(dir+'model_parms.txt') as f:
        e_list = f.read().split('@')[1].split('#')[1:] # split event list
        v = e_list[0].split('%')[1:] # split V choices
        v = {x.split(';')[2].rstrip('\n') : len(x.split(';')[1]) for x in v}

        d = e_list[1].split('%')[1:] # split D choices
        d = {x.split(';')[2].rstrip('\n') : len(x.split(';')[1]) for x in d}

        dct = {'Vlen':v, 'Dlen':d}

        names = ['delV3', 'delD3', 'delD5', 'delJ5', 'insVD', 'insDJ']
        for i in range(3,9):
            # split all the remaining featuresdct[names[i-3]] = {x.split(';')[1].rstrip('\n') : x.split(';')[0] for x in l}

            l = e_list[i].split('%')[1:]
            dct[names[i-3]] = {x.split(';')[1].rstrip('\n') : int(x.split(';')[0]) for x in l}

        return dct


def add_region_coords(df):
    # Extract CDR/WFR coordinates from segment_prop.csv and add to df with generated clonotypes

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


def read_igor_output(dir):
    # Read files with event choices and generated sequences
    # Add region and segment mapping

    s = pd.read_csv(dir + 'generated/generated_seqs_noerr.csv', sep=';')
    r = pd.read_csv(dir + 'generated/generated_realizations_noerr.csv', sep=';')

    r['V'] = r[r.columns[1]].apply(lambda x: int(x.strip('()')))
    r['J'] = r[r.columns[2]].apply(lambda x: int(x.strip('()')))
    r['Vlen'] = r[r.columns[1]].apply(lambda x: x.strip('()'))
    r['Dlen'] = r[r.columns[3]].apply(lambda x: x.strip('()'))
    r['delD5'] = r[r.columns[4]].apply(lambda x: x.strip('()'))
    r['delD3'] = r[r.columns[5]].apply(lambda x: x.strip('()'))
    r['delV3'] = r[r.columns[6]].apply(lambda x: x.strip('()'))
    r['delJ5'] = r[r.columns[7]].apply(lambda x: x.strip('()'))
    r['insVD'] = r[r.columns[8]].apply(lambda x: x.strip('()'))
    r['insDJ'] = r[r.columns[10]].apply(lambda x: x.strip('()'))

    s['sequence'] = s['nt_sequence']
    s['index'] = s['seq_index']

    df = pd.DataFrame()
    parms_dict = read_parms_file(dir)
    for i in parms_dict.keys():
        df[i] = r[i].map(parms_dict[i])

    df = pd.concat([s[['index', 'sequence']], r[['V', 'J']], df], axis=1)
    df = add_region_coords(df)

    df['v.end'] = df['Vlen'] - df['delV3']
    df['n1.end'] = df['v.end'] + df['insVD']
    df['d.end'] = df['n1.end'] + df['Dlen'] - df['delD3'] - df['delD5']
    df['n2.end'] = df['d.end'] + df['insDJ']

    # Filter out-of-frames
    in_frame = [len(x)%3==0 for x in df['sequence']]
    df = df[in_frame]
    df = df.reset_index(drop=True)

    return df


def set_counts(df, sampling_rate):
    n = len(df.index)
    counts = np.random.multinomial(n=int(sampling_rate * n), pvals=[1./n]*n, size=1)
    df['clone_count'] = pd.Series(counts[0, :])
    df = df.loc[df['clone_count'] > 0]
    return df


def make_rep(chain, model, size, dir, mr, smpl_rate):
    mutation_rate = 0.005 if mr is None else mr
    smpl_rate = 1 if smpl_rate is None else smpl_rate
    #prop_after_sampling = 0.092 * (smpl_rate ** 2) + smpl_rate * 0.928 + 0.005

    # if chain.startswith('bcr'):
    #     n = find_root_number(size, mutation_rate) / prop_after_sampling
    # else:
    #n = size / prop_after_sampling
    n = size

    # Generate initial repertoire
    in_frame_ratio = 3
    try:
        run_igor(size=n*in_frame_ratio, model=chain, dir=dir)
    except IgorError:
        raise

    rep = read_igor_output(dir)

    # Introduce SHMs and simulate cell proliferation
    if chain.startswith('bcr'):
        rep = mutate_repertoire(mutation_rate, model, rep)

    # Add random counts and downsample
    rep.to_csv(args.dir + 'clones_before_sampling.csv', sep='\t')
    rep = set_counts(rep, smpl_rate)
    rep = rep.reset_index(drop=True)

    return rep


def write_fasta(df, path):
# Creates fasta input for sequencing simulator

    lines = []

    if 'clone_count' in df:
        cc = df['clone_count']
    else:
        cc = [1] * len(df.index)

    for i in range(len(df.index)):
        for j in range(cc[i]):
            lines.append('>' + 'clone_' + str(df['index'][i]) + '_' + str(j) + '\n' + df['sequence'][i] + '\n')

    with open(path, 'w') as f:
        f.write(''.join(lines))


def write_fastq(df, dir):
# Creates two fastq files (MiXCR input)

    with open(dir+'reads/1.fq', 'w') as f:
        for i in range(len(df.index)):
            f.write('@' + 'clone_' + str(df['index'][i]) + '\n' + df['sequence'][i] + '\n+\n' + 'C'*len(df['sequence'][i]) + '\n')

    with open(dir+'reads/2.fq', 'w') as f:
        for i in range(len(df.index)):
            f.write('@' + 'clone_' + str(df['index'][i]) + '\n' + str(Seq(df['sequence'][i]).reverse_complement()) + '\n+\n' + 'C'*len(df['sequence'][i]) + '\n')


def run_sequencer(dir, sequencing_mode):
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


def make_reads_files(df, dir, sequencing_mode):
    if os.path.isdir(dir + "reads") == False:
        call(['mkdir', dir + 'reads'])

    if sequencing_mode in ['RNASeq', 'MiSeq', 'HiSeq']:
        write_fasta(df, dir+'clones.fasta')
        run_sequencer(args.dir, args.sequencing_mode)
    else:
        write_fastq(df, args.dir)


def run_mixcr(dir):
    if not os.path.isdir(dir + "mixcr"):
        call(["mkdir", dir+"mixcr"])

    receptor = 'igh' if args.model.startswith('bcr') else 'tcr'

    call(["mixcr", "analyze", "amplicon",
        "-s", "hsa",
        "--starting-material", "dna",
        "--5-end", "no-v-primers",
        "--3-end", "c-primers",
        "--adapters", "no-adapters",
        "--receptor-type", receptor,
        "--contig-assembly",
        "--align", '"--library default"',
        "--assemble", '"-OmaxBadPointsPercent=0"',
        "--assemble", '"-OcloneClusteringParameters=null"',
        "--assemble", '"-OassemblingFeatures="[FR1+CDR1+FR2+CDR2+FR3+CDR3+FR4]""',
        dir+"reads/1.fq", dir+"reads/2.fq", dir+"mixcr/analysis"])


parser = argparse.ArgumentParser(description="Simulate B or T cell receptor repertoire")
parser.add_argument("-c", dest="chain", type=str, help='Immunoglobulin chain (bcr_heavy/tcr_alpha/tcr_beta)', default='bcr_heavy')
parser.add_argument("-m", dest="model", type=str, help='Substitution model (steele/s5f)', default='steele')
parser.add_argument("-n", dest="size", type=int, help='Approximate size of the repertoire (minimal)')
parser.add_argument("-d", dest="dir", type=str, help='Output directory')
parser.add_argument("-mr", dest="mutation_rate", type=float, help='Bulk mutation rate')
parser.add_argument("-sr", dest="smpl_rate", type=float, help='Sampling rate')
parser.add_argument("-s", dest="sequencing_mode", type=str, help='Type of sequencing data (MiSeq/HiSeq/RNASeq/None')
args = parser.parse_args()


rep = make_rep(args.chain, args.model, args.size, args.dir, args.mutation_rate, args.smpl_rate)
print('Repertoire is generated')
#rep.to_csv(args.dir + 'clones_init.csv', sep='\t')


#make_reads_files(rep, args.dir, args.sequencing_mode)
print('Reads are generated')


#run_mixcr(args.dir)
print('MiXCR processing is done')
