import json
import pandas as pd
import numpy as np

#these functions are taken from here:
#https://github.com/kathryn1995/pmo_converter/blob/main/src/transformer.py

def transform_mhap_info(df, bioinfo_id, field_mapping, additional_hap_detected_cols=None):
    """Reformat the DataFrame based on the provided field mapping."""
    # renamed_columns = {col: field_mapping[col]
    #                    for col in field_mapping if field_mapping[col] != "None"}
    # transformed_df = df.rename(columns=renamed_columns)
    transformed_df = microhaplotype_table_to_pmo_dict(
        df, bioinfo_id, sampleID_col=field_mapping["sampleID"], locus_col=field_mapping['locus'], mhap_col=field_mapping['asv'], reads_col=field_mapping['reads'], additional_hap_detected_cols=additional_hap_detected_cols)
    return transformed_df


def transform_panel_info(df, panel_id, field_mapping, target_genome_info, additional_target_info_cols=None):
    """Reformat the DataFrame based on the provided field mapping."""
    transformed_df = panel_info_table_to_pmo_dict(
        df,
        panel_id,
        target_genome_info,
        target_id_col=field_mapping["target_id"],
        forward_primers_seq_col=field_mapping["forward_primers"],
        reverse_primers_seq_col=field_mapping["reverse_primers"],
        additional_target_info_cols=additional_target_info_cols)
    return transformed_df


def microhaplotype_table_to_pmo_dict(
    contents: pd.DataFrame,
    bioinfo_id: str,
    sampleID_col: str = 'sampleID',
    locus_col: str = 'locus',
    mhap_col: str = 'asv',
    reads_col: str = 'reads',
    additional_hap_detected_cols: list | None = None
):
    """
    Convert a dataframe of a microhaplotype calls into a dictionary containing a dictionary for the haplotypes_detected and a dictionary for the representative_haplotype_sequences.

    :param contents: The dataframe containing microhaplotype calls
    :param bioinfo_id: the bioinformatics ID of the microhaplotype table
    :param sampleID_col: the name of the column containing the sample IDs
    :param locus_col: the name of the column containing the locus IDs
    :param mhap_col: the name of the column containing the microhaplotype sequence
    :param reads_col: the name of the column containing the reads counts
    :param additional_hap_detected_cols: optional additional columns to add to the microhaplotype detected dictionary, the key is the pandas column and the value is what to name it in the output
    :return: a dict of both the haplotypes_detected and representative_haplotype_sequences
    """

    representative_microhaplotype_dict = create_representative_microhaplotype_dict(
        contents, locus_col, mhap_col)

    detected_mhap_dict = create_detected_microhaplotype_dict(contents, sampleID_col, locus_col,
                                                             mhap_col, reads_col, representative_microhaplotype_dict,
                                                             additional_hap_detected_cols)

    output_data = {"microhaplotypes_detected": {bioinfo_id: {'experiment_samples': detected_mhap_dict}},
                   "representative_microhaplotype_sequences": {bioinfo_id: {"representative_microhaplotype_id": bioinfo_id, 'targets': representative_microhaplotype_dict}}
                   }
    output_data = json.dumps(output_data, indent=4)
    return output_data


def create_representative_microhaplotype_dict(
        microhaplotype_table: pd.DataFrame,
        locus_col: str,
        mhap_col: str
):
    """
    Convert the read-in microhaplotype calls table into a representative microhaplotype JSON-like dictionary.

    :param microhaplotype_table: The parsed microhaplotype calls table.
    :param locus_col: The name of the column containing the locus IDs.
    :param mhap_col: The name of the column containing the microhaplotype sequence.
    :return: A dictionary formatted for JSON output with representative microhaplotype sequences.
    """
    # Drop duplicates and reset index
    unique_table = microhaplotype_table[[
        locus_col, mhap_col]].drop_duplicates().reset_index(drop=True)

    json_data = {
        locus: {
            "seqs": {
                f"{locus}.{idx}": {
                    "microhaplotype_id": f"{locus}.{idx}",
                    "seq": seq
                }
                for idx, seq in enumerate(group[mhap_col])
            }
        }
        for locus, group in unique_table.groupby(locus_col)
    }

    return json_data


def create_detected_microhaplotype_dict(
    microhaplotype_table: pd.DataFrame,
    sampleID_col: str,
    locus_col: str,
    mhap_col: str,
    reads_col: str,
    representative_microhaplotype_dict: dict,
    additional_hap_detected_cols: list | None = None
):
    """
    Convert the read-in microhaplotype calls table into the detected microhaplotype dictionary.

    :param microhaplotype_table: Parsed microhaplotype calls table.
    :param sampleID_col: Column containing the sample IDs.
    :param locus_col: Column containing the locus IDs.
    :param mhap_col: Column containing the microhaplotype sequences.
    :param reads_col: Column containing the read counts.
    :param representative_microhaplotype_dict: Dictionary of representative microhaplotypes.
    :param additional_hap_detected_cols: Optional additional columns to add to the microhaplotypes detected, the key is the pandas column and the value is what to name it in the output.
    :return: A dictionary of detected microhaplotype results.
    """
    # Validate additional columns if provided
    if additional_hap_detected_cols:
        check_additional_columns_exist(
            microhaplotype_table, additional_hap_detected_cols)

    # Map sequences to representative haplotype IDs for fast lookup
    rep_hap_map = {
        (locus, seq["seq"]): seq["microhaplotype_id"]
        for locus, reps in representative_microhaplotype_dict.items()
        for seq in reps["seqs"].values()
    }

    def build_haplotype_info(row):
        """Helper to construct haplotype info for each row."""
        locus, seq = row[locus_col], row[mhap_col]
        matching_id = rep_hap_map.get((locus, seq))
        if not matching_id:
            raise ValueError(
                f"No representative haplotype ID found for {seq} at locus {locus}")

        haplotype_info = {
            "haplotype_id": matching_id,
            "read_count": row[reads_col],
        }

        if additional_hap_detected_cols:
            haplotype_info.update({input_col: row[input_col]
                                   for input_col in additional_hap_detected_cols})

        return matching_id, haplotype_info

    # Build the JSON-like structure
    json_data = (
        microhaplotype_table
        .groupby(sampleID_col)
        .apply(lambda sample_group: {
            "sample_id": sample_group[sampleID_col].iloc[0],
            "target_results": {
                locus: {
                    "microhaplotypes": {
                        hap_id: hap_info
                        for _, row in locus_group.iterrows()
                        for hap_id, hap_info in [build_haplotype_info(row)]
                    }
                }
                for locus, locus_group in sample_group.groupby(locus_col)
            }
        })
        .to_dict()
    )
    return json_data


def check_additional_columns_exist(df, additional_column_list):
    if additional_column_list:
        missing_cols = set(additional_column_list) - \
            set(df.columns)
        if missing_cols:
            raise ValueError(f"Missing additional columns: {missing_cols}")


def panel_info_table_to_pmo_dict(target_table: pd.DataFrame,
                                 panel_id: str,
                                 genome_info: dict,
                                 target_id_col: str = 'target_id',
                                 forward_primers_seq_col: str = 'fwd_primer',
                                 reverse_primers_seq_col: str = 'rev_primer',
                                 forward_primers_start_col: int | None = None,
                                 forward_primers_end_col: int | None = None,
                                 reverse_primers_start_col: int | None = None,
                                 reverse_primers_end_col: int | None = None,
                                 insert_start_col: int | None = None,
                                 insert_end_col: int | None = None,
                                 chrom_col: str | None = None,
                                 strand_col: str | None = None,
                                 gene_id_col: str | None = None,
                                 target_type_col: str | None = None,
                                 additional_target_info_cols: list | None = None,
                                 ):
    """
    Convert a dataframe containing panel information into dictionary of targets and reference information


    :param target_table: The dataframe containing the target information
    :param panel_id: the panel ID assigned to the panel
    :param genome_info: A dictionary containing the genome information
    :param target_id_col: the name of the column containing the target IDs
    :param forward_primers_seq_col: the name of the column containing the sequence of the forward primer
    :param reverse_primers_seq_col: the name of the column containing the sequence of the reverse primer
    :param forward_primers_start_col (Optional): the name of the column containing the 0-based start coordinate of the forward primer
    :param forward_primers_end_col (Optional): the name of the column containing the 0-based end coordinate of the forward primer
    :param reverse_primers_start_col (Optional): the name of the column containing the 0-based start coordinate of the reverse primer
    :param reverse_primers_end_col (Optional): the name of the column containing the 0-based end coordinate of the reverse primer
    :param insert_start_col (Optional): the name of the column containing the 0-based start coordinate of the insert
    :param insert_end_col (Optional): the name of the column containing the 0-based end coordinate of the insert
    :param chrom_col (Optional): the name of the column containing the chromosome for the target
    :param gene_id_col (Optional): the name of the column containing the gene id
    :param strand_col (Optional): the name of the column containing the strand for the target
    :param target_type_col (Optional): A classification type for the target
    :param additional_target_info_cols (Optional): dictionary of optional additional columns to add to the target information dictionary. Keys are column names and values are the type.
    :return: a dict of the panel information
    """

    if not isinstance(target_table, pd.DataFrame):
        raise ValueError("target_table must be a pandas DataFrame.")
    if not isinstance(genome_info, dict):
        raise ValueError("genome_info must be a dictionary.")

    # Check additional columns if any are added
    check_additional_columns_exist(target_table, additional_target_info_cols)

    # If one location column set, check all location columns are set
    location_cols = check_location_columns(forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col,
                                           reverse_primers_end_col, insert_start_col, insert_end_col, chrom_col, strand_col)

    # Create dictionary of targets
    targets_dict = create_targets_dict(target_table, target_id_col, forward_primers_seq_col, reverse_primers_seq_col,
                                       location_info_cols=location_cols, gene_id_col=gene_id_col, target_type_col=target_type_col,
                                       additional_target_info_cols=additional_target_info_cols)
    # Put together components
    panel_info_dict = {"panel_info": {panel_id: {"panel_id": panel_id,
                                                 "target_genome": genome_info, "targets": targets_dict}}}
    # Convert to json format
    panel_info_json = json.dumps(panel_info_dict, indent=4)
    return panel_info_json


def create_targets_dict(
    target_table: pd.DataFrame,
    target_id_col: str,
    forward_primers_seq_col: str,
    reverse_primers_seq_col: str,
    location_info_cols: list | None = None,
    gene_id_col: str = None,
    target_type_col: str = None,
    additional_target_info_cols: list = None
):
    """
    Convert the read in target information into the panel information dictionary


    :param target_table: The dataframe containing the target information
    :param target_id_col: the name of the column containing the target IDs
    :param forward_primers_seq_col: the name of the column containing the sequence of the forward primer
    :param reverse_primers_seq_col: the name of the column containing the sequence of the reverse primer
    :param location_info_cols: list of column names for location information
    :param gene_id_col: the name of the column containing the gene id
    :param chrom_col: the name of the column containing the chromosome the set of primers target
    :param strand_col: the name of the column containing the strand for the target
    :return: a dictionary of the target information
    """
    # Check targets before putting into JSON
    columns_to_check = []
    if location_info_cols:
        forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col, reverse_primers_end_col, insert_start_col, insert_end_col, chrom_col, strand_col = location_info_cols
        columns_to_check = columns_to_check + \
            [chrom_col, strand_col, insert_start_col, insert_end_col]
    if gene_id_col:
        columns_to_check.append(gene_id_col)
    if target_type_col:
        columns_to_check.append(target_type_col)
    if additional_target_info_cols:
        columns_to_check = columns_to_check + additional_target_info_cols

    # Check these columns are unique for each target
    check_columns_unique_for_target(
        target_table, target_id_col, columns_to_check)

    # Put targets together in dictionary
    targets_dict = {}
    for target_id in target_table[target_id_col].unique():
        target_df = target_table.loc[target_table[target_id_col] == target_id]
        target_dict = {
            "target_id": target_id,
            "forward_primers": [],
            "reverse_primers": [],
        }
        if gene_id_col:
            target_dict["gene_id"] = target_df[gene_id_col].iloc[0]
        if target_type_col:
            target_dict["target_type"] = target_df[target_type_col].iloc[0]
        if additional_target_info_cols:
            for col in additional_target_info_cols:
                value = target_df[col].iloc[0]
                # Convert numpy types to native Python types
                if isinstance(value, (np.integer, np.int64)):
                    value = int(value)
                elif isinstance(value, (np.floating, np.float64)):
                    value = float(value)
                elif pd.isna(value):
                    value = None
                target_dict[col] = value
        # Add insert location info if location_info_cols are provided
        if location_info_cols:
            target_dict["insert_location"] = {
                "chrom": target_df[chrom_col].iloc[0],
                "start": int(target_df[insert_start_col].iloc[0]),
                "end": int(target_df[insert_end_col].iloc[0]),
                "strand": target_df[strand_col].iloc[0]
            }
        # Extract primer information for each row
        for _, row in target_df.iterrows():
            fwd_primer_dict = {"seq": row[forward_primers_seq_col]}
            rev_primer_dict = {"seq": row[reverse_primers_seq_col]}
            if location_info_cols:
                fwd_primer_dict["location"] = {
                    "chrom": row[chrom_col],
                    "end": int(row[forward_primers_end_col]),
                    "start": int(row[forward_primers_start_col]),
                    "strand": row[strand_col]
                }
                rev_primer_dict["location"] = {
                    "chrom": row[chrom_col],
                    "end": int(row[reverse_primers_end_col]),
                    "start": int(row[reverse_primers_start_col]),
                    "strand": row[strand_col]
                }
            target_dict["forward_primers"].append(fwd_primer_dict)
            target_dict["reverse_primers"].append(rev_primer_dict)

        targets_dict[target_id] = target_dict
    return targets_dict


def check_location_columns(
    forward_primers_start_col,
    forward_primers_end_col,
    reverse_primers_start_col,
    reverse_primers_end_col,
    insert_start_col,
    insert_end_col,
    chrom_col,
    strand_col,
):
    location_cols = [
        forward_primers_start_col,
        forward_primers_end_col,
        reverse_primers_start_col,
        reverse_primers_end_col,
        insert_start_col,
        insert_end_col,
        chrom_col,
        strand_col,
    ]
    if any(location_cols):
        if not all(location_cols):
            raise ValueError(
                "If one of the location params (forward_primers_start_col, forward_primers_end_col, reverse_primers_start_col, reverse_primers_end_col, insert_start_col, insert_end_col, chrom_col, strand_col) is set then all must be set."
            )
        return location_cols
    return None


def check_columns_unique_for_target(df, target_id_col, columns_to_check):
    for col in columns_to_check:
        duplicates = df.groupby(target_id_col)[col].nunique()
        duplicates = duplicates[duplicates > 1]
        if not duplicates.empty:
            duplicate_targets = duplicates.index.tolist()
            raise ValueError(
                f"The following target_ids have multiple unique {col}: {duplicate_targets}")
