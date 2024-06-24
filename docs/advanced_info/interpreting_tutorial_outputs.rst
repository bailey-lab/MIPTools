.. _advanced_wrangler_interpretation:

Interpreting Wrangler Tutorial Output
=====================================

If you'd like to see the 'raw' outputs of the wrangler, the main output file is
called allInfo.tsv.gz. After unzipping this file (e.g. with gunzip) you can
open the file in a spreadsheet program. This file has a few essential columns:

- s_Sample: the name of the sample
- p_geneName: the gene being analyzed
- p_targetName: the MIP being used to target a section of the gene
- h_popUID: the name of the haplotype that was pulled for a particular MIP on a
  particular sample
- h_seq: the sequence of the haplotype
- c_barcodeCnt: the number of UMIs that support the haplotype pulled for a MIP
  for a sample

Every row represents a UMI count of a single haplotype of a single MIP of a
single sample. Although there are multiple rows associated with any given
sample, each row of a given sample is usually a different MIP. When there are
multiple rows that share the same sample and the same MIP, this indicates
multiple haplotypes.

In the tutorial dataset, the sample AG-00-01-PRX-00-1 has 108 rows associated
with it, but only 5 of the MIPs in these rows occur more than once
(crt_S0_Sub0_mip12, crt_S0_Sub0_mip20, dhps_S0_Sub0_mip14,
k13_S0_Sub0_mip8, and mdr1_S0_Sub0_mip15). For the MIPs that do occur more than
once, these represent different haplotypes. For example, for the rows that have
AG-00-01-PRX-00-1 as the sample, the subset of rows that have k13_S0_Sub0_mip8
are associated with two different haplotypes, k13_S0_Sub0_mip8.00 and
k13_S0_Sub0_mip8.08. The k13_S0_Sub0_mip8.00 row has a c_barcodeCnt (UMI count)
of 4, while k13_S0_Sub0_mip8.08 row has a c_barcodeCnt of 1. This means that in
the sample AG-00-01-PRX-00-1, 80 percent of the k13_S0_Sub0_mip8 data was
associated with haplotype 00, while 20 percent was associated with haplotype 08

An alignment of the two haplotypes (from the h_seq column) reveals that there
are only a few point mutations that separate these haplotypes:
haplotype 00: ATTATCAAGAATATAAAAATTTTGAGAATGATAAA
haplotype 08: ATTATGAACTAAATTAAAATTTTGCTCAATAAAAA

By aligning these haplotypes to the genome (not shown here) you might find that
the upper haplotype maps perfectly to the 3D7 genome, while the lower haplotype
has 11 mutations that affect 7 codons.

.. _advanced_umi_countvs_probe_coverage_interpretation:

Interpreting umi_counts
=======================

In the document umi_count_vs_probe_coverage.html, the x-axis is UMI count and
the y-axis is number of MIPs with >10 UMIs. Because DR23K has 121 UMIs, the
maximum value on the y-axis is 121. As some example samples:

  - AM-07-89-PRX-07-1 (zoom in as far as possible on the bottom left cluster of
    points and hover over the leftmost point within the cluster) has no MIPs
    (out of 121) with greater than 10 UMIs. If you hover over this point,
    you'll see a UMI count of 69, and a read count of 70. This means that every
    UMI in this sample was sequenced an average of just over one time. More
    reads might easily recover more UMIs, which could lead to more MIPs being
    sequenced with at least 10 UMIs each. It may be possible to recover these
    additional reads by skipping over samples that previously had high read
    counts, and repooling and re-sequencing samples like this one to boost read
    counts.
  - KB-04-06-PRX-04-1 (zoom in on the bottom left and hover over rightmost
    point) also has no MIPs (out of 121) with greater than 10 UMIs. However, if
    you hover over this point, you'll see a UMI count of 209, and a read count
    of 5,600. This means that every UMI in this sample was sequenced an average
    of 26.7 times (5,600/209). The number of reads per UMI is plenty, and more
    sequencing would probably not recover more UMIs, so re-doing the MIP
    capture reaction will probably work better.
  - KO-07-001-PRX-07-1 (top right) has 117 MIPs (out of 121) that have greater
    than 10 UMIs each. Hovering on the dot, the UMI count for this sample was
    622,789 and the read count for this sample was 1,310,326. This means there
    were an average of 2.104 reads per MIP (1,310,326/622,789), or an average of
    5,147 UMIs per MIP (622,789/121). This sample has plenty of UMIs for each
    MIP, so it's only the MIPs that perform poorly across the entire dataset
    that will be missing here.