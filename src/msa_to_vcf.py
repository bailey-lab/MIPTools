import mip_functions as mip
import pandas as pd


def msa_to_vcf(alignment_file, vcf_file, ref=None, snp_only=False):
    """ Take a multiple sequence alignment file and create a vcf file
    containing all variants.
    """
    # read in the alignment file
    aln = mip.fasta_parser(alignment_file)

    # convert each alignment string to list
    aln = {k: list(aln[k]) for k in aln}
    # create alignment dataframe
    aln_df = pd.DataFrame(aln)

    if not snp_only:
        # find indexes of all indels
        indel_index = aln_df.loc[aln_df.apply(
            lambda a: "-" in a.values, axis=1)].index
        if len(indel_index) == 0:
            variant_df = aln_df
        else:
            # merge all neighboring indel indexes
            indel_index = [[i, i] for i in indel_index]
            indel_index = mip.merge_overlap(indel_index, spacer=1)
            # include the prior base for indel calls
            indel_index = [[i[0] - 1, i[1]] if i[0] > 0 else i
                           for i in indel_index]
            # get indexes for snps (non-indel changes)
            snp_sets = [set(range(indel_index[i][1] + 1, indel_index[i+1][0]))
                        for i in range(len(indel_index) - 1)]
            # include the positions from last indel to the end of alignment
            snp_sets.append(set(range(indel_index[-1][1] + 1, len(aln_df))))
            # include the positions from the beginning of the alignment to the
            # first indel
            snp_sets.append(set(range(0, indel_index[0][0])))
            # create a set of snp indexes
            snp_index = set()
            for s in snp_sets:
                snp_index.update(s)

            # go through each indel and create a dataframe for each where the
            # each haplotype has the allele for the indel at the beginning
            # index position
            indel_list = []
            for ind in indel_index:
                indel_list.append(
                    # slice the alignment df for the indel location
                    aln_df.loc[ind[0]: ind[1]].apply(
                        # join the nucleotides (remove "-") for each haplotype
                        lambda a: "".join([v for v in a.values if v != "-"])))
            # concatanate the indel dataframes
            indel_df = pd.concat(indel_list, axis=1).T
            # set index using the first position for each indel
            indel_df.index = [ind[0] for ind in indel_index]

            # remove snp indexes where all haplotypes are same as reference
            snp_index = sorted([ind for ind in snp_index
                                if len(set(aln_df.loc[ind].values)) > 1])

            # merge indel and snp dfs
            variant_df = pd.concat([indel_df, aln_df.loc[snp_index]])
    else:
        snp_index = sorted([ind for ind in aln_df.index if len(set(
            aln_df.loc[ind].values).difference(["-"])) > 1])
        snp_df = aln_df.loc[snp_index]
        variant_df = snp_df.loc[snp_df["ref"] != "-"].replace("-", "*")
    variant_df.sort_index(inplace=True)

    # if reference sequence name is provided, sort alleles by that
    if ref is not None:
        allele_df = variant_df.apply(
            lambda a: [a["ref"]] + list(set(a).difference([a["ref"]])), axis=1)
    else:
        allele_df = variant_df.apply(
            lambda a: list(set(a)), axis=1)
    allele_dict = allele_df.to_dict()

    def get_genotype(row):
        allele_list = allele_dict[row.name]
        return pd.Series([allele_list.index(v) for v in row.values])

    genotypes = variant_df.apply(get_genotype, axis=1)
    genotypes.columns = variant_df.columns
    allele_df.name = "alleles"
    vcf = pd.concat([genotypes, allele_df], axis=1)
    vcf.index.name = "POS"
    vcf = vcf.reset_index()
    vcf["POS"] = vcf["POS"] + 1
    vcf["REF"] = vcf["alleles"].apply(lambda a: a[0])
    vcf["ALT"] = vcf["alleles"].apply(lambda a: ",".join(a[1:]))
    fields = ["#CHROM", "ID", "QUAL", "FILTER", "INFO", "FORMAT"]
    for f in fields:
        vcf[f] = "."
    fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
              "FORMAT"] + list(genotypes.columns)
    vcf = vcf[fields]
    vcf.to_csv(vcf_file, sep="\t", index=False)
