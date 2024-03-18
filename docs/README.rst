a2ihelper
=========

**a2ihelper** is a Python library for downstream analysis of A-to-I editing. It is build under Python 3.10 to analyze statistics `REDItools2 <https://github.com/BioinfoUNIBA/REDItools2>`_ output files.
It is possible to analyze frequency or proportional data. The principal statistics tests used in 68 papers are included in a2iHelperPy.

Installing
----------
.. code-block:: bash

   pip install a2ihelper

Quickstart
----------
You need:
  - Mapped BAM files.
  - GTF annotation files (the same version used to mapping).
  - List of genes or coordinates of interest.

**Get coordinates from list of genes:**

.. code-block:: python

   genes = ['Notch1','B2m','Hdc','Il1b']

   path_ref_annotation='~/../reference_files/gencode.annotation.gtf'

   gene_coord = a2i.call_reditools2.get_genes_positions(genes, path_ref_annotation, gzip_file=False)
   gene_coord

.. code-block:: python

   {'Notch1': 'chr2:26457903-26516663', 'B2m': 'chr2:122147686-122153083',
   'Hdc': 'chr2:126593667-126619299', 'Il1b': 'chr2:129364570-129371139'}

**Running REDItools2 from a2ihelper**

.. code-block:: python

   genes_positions = gene_coord.values() # coordinates from gene_coord dictionary
   in_bam_file_list = ['/mnt/d/rna_editing/bams/'+f for f in os.listdir('/mnt/d/rna_editing/bams/') if f.endswith('bam')] #list of bam files
   path_out_res = '/mnt/d/rna_editing/a2i_helper_output/res/' # directory for output files
   ref_genome_file = '/mnt/d/reference_files/mus_musculus/GRCm39.genome.fa' # Reference Genome
   path_reditools = '/mnt/d/reditools2.0/src/cineca/' # directory where reditools.py is installed
   reditools_options = '--strict' # optional arguments separeted per single space

   for in_bam_file in in_bam_file_list:
      a2i.call_reditools2.run_per_gene_position_list(genes_positions, in_bam_file, path_out_res,
                                                     ref_genome_file, path_reditools,
                                                     reditools_options='', n_jobs=10)


It'll generate the RES files in the path_out_res directory.

Now you need to prepare your metadata file. The first column must be the path to the RES files, the second columns must be the sample name, third must be the region, fourth the condition, and the fifth the coordinates.

**Analyzing ONE region**

.. code-block:: python

   df, df_a, df_g = a2i.editing.merge_files_one_region(meta[meta.region=='B2m'])

The function *merge_files_one_region()* returns three DataFrames. First is the frequency of editing, second is the Adenine counts and the third the Guanine (Inosine) counts.

**Analyzing ALL region**

.. code-block:: python

   df, df_a, df_g, region_list = a2i.editing.merge_files_all_regions(meta)

The function *merge_files_all_regions()* returns three DataFrames and one list with the sequence of genes counts. First DataFrame is the frequency of editing, second is the Adenine counts and the third the Guanine (Inosine) counts.
