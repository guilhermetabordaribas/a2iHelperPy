a2ihelper
=========

**a2ihelper** is a Python library for downstream analysis of A-to-I editing. It is build under Python 3.10.13 to analyze statistics `REDItools2 <https://github.com/BioinfoUNIBA/REDItools2>`_ output files.
It is possible to analyze frequency or proportional data. The principal statistics tests used in 68 papers are included in a2iHelperPy.

Requirements
------------

  - numpy (version >= 1.26.4)
  - pandas (version >= 2.2.0)
  - scipy (version >= 1.12.0)
  - scikit-posthocs (version >= 0.9.0)


  - **Optional plot:**

    - matplotlib (version >= 3.8.3)
    - seaborn (version >= 0.13.2)


  - **Optionals Editing detection:**

    - reditools (version = 2.0)
    - pysam (versoin >= 0.22.0)
    - sortedcontainers (version >= 2.4.0)
    - psutil (version >= 5.9.8)
    - netifaces (version >= 0.11.0)

Installing
----------
We recommend use pip to install the package

.. code-block:: bash

   pip install a2ihelper

But you can also use the git clone and install via setup.py

Quickstart
----------
You need:
  - Mapped BAM files.
  - GTF annotation files (the same version used to mapping).
  - List of genes or coordinates of interest.

**Get coordinates from list of genes:**

You need to inform a list of genes that you want to get the coordinates and the GTF annotated path file (Make sure to set ``gzip_file=True`` if the file is "gzipped"). It is important to use the same file used to mapping.

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

**Statistics for FREQUENCY**

*Mann-Whitney U test*

.. code-block:: python

   df_pv = a2i.editing.mannwhitney_test(df,
                                        only_pvalue=True,
                                        pvalue_filter_limit=0.05,
                                        fdr_correction=True,
                                        fdr_filter_limit=0.05,
                                        return_only_significant=True)

*ANOVA Tukey*

.. code-block:: python

   df_pv = a2i.editing.anova_tukey_test(df,
                                        only_pvalue=True,
                                        pvalue_filter_limit_anova=0.05,
                                        pvalue_filter_limit_tukey=0.05,
                                        return_only_significant=True)

*Kruskal Dunn*

.. code-block:: python

   df_pv = a2i.editing.kruskal_dunn_test(df,
                                         only_pvalue=True,
                                         pvalue_filter_limit_kruskal=0.05,
                                         pvalue_filter_limit_dunn=0.05,
                                         return_only_significant=True)

**Statistics for PROPORTION**

*Pooling replicates*

May you need to pool the replicates to perform a chi-square or fisher tests. To do that you can use the ``pool_positions()`` to sum the coordinates by independecy G-test.

.. code-block:: python

   a, g = a2i.editing.pool_positions(df_a,
                                     df_g,
                                     pvalue_filter_limit=0.05,
                                     gtest_filter_limit=0,
                                     bh_correction=False)

*Chi-square test*

.. code-block:: python

   chi = a2i.editing.chi2_test(a, g,
                               only_pvalue=True,
                               pvalue_filter_limit=0.05)

*Fisher test*

.. code-block:: python

   fis = a2i.editing.fisher_test(a, g,
                                 only_pvalue=True,
                                 pvalue_filter_limit=.05)
