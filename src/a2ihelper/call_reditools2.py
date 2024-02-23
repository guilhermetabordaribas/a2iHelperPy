import os
import subprocess
from multiprocessing import Pool
import warnings
import itertools
import gzip
# import pandas as pd

def run_per_gene_position(gene_position: str, in_bam_file: str, path_out_res: str, ref_genome_file: str, path_reditools: str, reditools_options: str) -> None:
    """
    Run reditools2.0 via Python

    Parameters
    ----------
    gene_position : str
        coordinate of chromossome:start-end positions to search editing sites.
        Example: 'chr2:122147686-122153083'. One can get the coordinates from gene symbol using get_genes_positions function.
    in_bam_file : str
        full input BAM file path
    path_out_res: str
        output directory
    ref_genome_file: str
        full reference FASTA file path. Must the same used to build the aligned bam file.
    path_reditools: str
        full directory where reditools.py is installed. Usually is in similiar path '/../reditools2.0/src/cineca'
    reditools_options: str
        optional arguments to run reditools2. All the options are expalined in https://github.com/BioinfoUNIBA/REDItools2

    Returns
    -------
    None
        it doesn't return nothing, just run reditools2

    Example
    -------
        Using the toyfile from https://github.com/guilhermetabordaribas/a2iHelperPy

        >>> gene_position = 'chr2:122147686-122153083'
        >>> in_bam_file = '/.../sample1.sortedByCoord.out.bam'
        >>> path_out_res = '/.../out/'
        >>> ref_genome_file = '/.../GRCh38.p14.genome.fa'
        >>> path_reditools = '/.../reditools2.0/src/cineca/'
        >>> reditools_options = '--strict'
        >>> run_per_gene_position(gene_position, in_bam_file, path_out_res, ref_genome_file, path_reditools, reditools_options='--strict')
    """

    out_file = os.path.join( path_out_res, os.path.basename(in_bam_file)+'_'+gene_position.replace(':','_').replace('-','_')+'_RES.tsv' )
    cmd_list = ['python', 'reditools.py', '-f', in_bam_file, '-r', ref_genome_file, '-o', out_file, '-g', gene_position]
    if reditools_options:
        cmd_list += reditools_options.split(' ')
    subprocess.call(cmd_list, cwd=path_reditools, stdout=subprocess.PIPE)

def run_per_gene_position_list(genes_positions: list, in_bam_file: str, path_out_res: str, ref_genome_file: str, path_reditools: str, reditools_options: str, n_jobs=4) -> None:
    """
    Run run_per_gene_position for a list of gense coordinates (genes_positions)

    Parameters
    ----------
    gene_position: list
        list of coordinates of chromossome:start-end positions to search editing sites.
        Example: ['chr2:122147686-122153083', 'chr18:60803848-60812646', 'chr6:65671590-65712326'].
        One can get the coordinates from gene symbol using get_genes_positions function.
    in_bam_file: str
        full input file BAM path
    path_out_res: str
        output directory
    ref_genome_file: str
        full reference FASTA file path. Must the same used to build the aligned bam file.
    path_reditools: str
        full directory where reditools.py is installed. Usually is in similiar path '/../reditools2.0/src/cineca'
    reditools_options: str
        optional arguments to run reditools2. All the options are expalined in https://github.com/BioinfoUNIBA/REDItools2
    n_jobs: int
        number of jobs in parallel

    Returns
    -------
    None
        it doesn't return nothing, just run reditools2 for a list o coordinates

    Example
    -------
        Using the toyfile from https://github.com/guilhermetabordaribas/a2iHelperPy

        >>> genes_positions = ['chr2:122147686-122153083', 'chr6:65671590-65712326', 'chr15:78191114-78206400']
        >>> in_bam_file = '/.../sample1.sortedByCoord.out.bam'
        >>> path_out_res = '/.../out/'
        >>> ref_genome_file = '/.../GRCh38.p14.genome.fa'
        >>> path_reditools = '/.../reditools2.0/src/cineca/'
        >>> reditools_options = '--strict'
        >>> run_per_gene_position_list(genes_positions, in_bam_file, path_out_res, ref_genome_file, path_reditools, reditools_options='--strict', n_jobs=4)
    """

    arguments_list = zip(genes_positions, itertools.repeat(in_bam_file), itertools.repeat(path_out_res), itertools.repeat(ref_genome_file), itertools.repeat(path_reditools), itertools.repeat(reditools_options))
    with Pool(processes=n_jobs) as p:
        p.starmap(run_per_gene_position, arguments_list)

def get_genes_positions(genes, path_ref_annotation, gzip_file=True) -> list:
    """
    Return the coordinates of a gene symbol from a GTF file. It can be used as input to run_per_gene_position_list.

    Parameters
    ----------
    genes: list
        list of genes to get coordinates
    path_ref_annotation: str
        full reference GTF file path.

    Returns
    -------
        dict
            a dict of coordinates of each gene symbol (chr:start-end).

    Example
    -------
        Using the GTF file from gencode

        >>> get_genes_positions(['B2m', 'Apol1'], '/.../GRCh38.p14.genome.fa')
        ['chr2:122147686-122153083', 'chr18:60803848-60812646']
    """

    genes_positions_list = []
    gens_aux = []
    if genes:
        if gzip_file:
            # for g in genes:
            with gzip.open(path_ref_annotation,'r') as f_gtf:
                for line in f_gtf:
                    if not line.startswith('#'.encode()):
                        l = line.decode().split('\t')
                        dict_g = { i.split(' ')[0]:i.split(' ')[1] for i in [j.strip() for j in l[-1].replace(';\n','').replace('"','').split(';')] }
                        if (l[2]=='gene') and (dict_g['gene_name'] in genes):
                            g_list.append(l[0]+':'+l[3]+'-'+l[4])
                            gens_aux.append( dict_g['gene_name'] )
        else:
            # for g in genes:
            with open(path_ref_annotation,'r') as f_gtf:
                for line in f_gtf:
                    if not line.startswith('#'):
                        l = line.split('\t')
                        dict_g = { i.split(' ')[0]:i.split(' ')[1] for i in [j.strip() for j in l[-1].replace(';\n','').replace('"','').split(';')] }
                        if (l[2]=='gene') and (dict_g['gene_name'] in genes):
                            genes_positions_list.append(l[0]+':'+l[3]+'-'+l[4])
                            gens_aux.append( dict_g['gene_name'] )

    if not genes_positions_list:
        warnings.warn('*Returning empty list* -> Positions of genes were not found in the '+ path_ref_annotation+'. Please verify genes names or gtf file.')

    return dict(zip(gens_aux,genes_positions_list))


def indexing_ref(path_ref_genome):
    pass
