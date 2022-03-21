import mip_functions as mip
import argparse

# Read command line arguments
parser = argparse.ArgumentParser(
    description="""Call variants using Freebayes."""
)
parser.add_argument(
    "-b",
    "--bam-dir",
    help="bam directory to create or use bam files in.",
    default="/opt/analysis/padded_bams",
)
parser.add_argument(
    "-d",
    "--fastq-dir",
    help="fastq directory to create or use fastqs in.",
    default="/opt/analysis/padded_fastqs",
)
parser.add_argument(
    "-o",
    "--output-vcf",
    help="vcf.gz file name to create.",
    default="/opt/analysis/variants.vcf.gz",
)
parser.add_argument(
    "-g",
    "--targets-file",
    help="tab separated targets file for variant calling.",
)
parser.add_argument(
    "-k",
    "--skip-fastq",
    help=("set this flag to skip creating fastq files from MIP data."),
    action="store_false",
)
parser.add_argument(
    "-p",
    "--skip-align",
    help="set this flag to skip bwa alignment to genome.",
    action="store_false",
)
parser.add_argument(
    "-z",
    "--bam-files",
    nargs="*",
    help=(
        "bam file to use for variant calling."
        "Use multiple times for more than one file."
    ),
)
parser.add_argument(
    "-l",
    "--bam-list",
    help="file containing absolute paths to bam files to use.",
)


parser.add_argument(
    "-s",
    "--settings-file",
    help="MIPTools analysis settings file to use.",
    default="/opt/analysis/settings.txt",
)
parser.add_argument(
    "-f",
    "--fastq-padding",
    help="number of reference genome bases to flank " "haplotypes.",
    type=int,
    default=20,
)
parser.add_argument(
    "-q",
    "--min-base-quality",
    help="minimum base qual to consider an allele",
    type=int,
    default=1,
)
parser.add_argument(
    "-t",
    "--threads",
    help="number of CPU threads to use.",
    type=int,
    default=20,
)
parser.add_argument(
    "-e",
    "--extra-freebayes-options",
    nargs="*",
    help=(
        "additional Freebayes options to pass directly to Freebayes. Options "
        "must have + in place of -. For example, ++pooled+continuous if you "
        "want to pass --pooled-continuous."
    ),
)
parser.add_argument(
    "-w",
    "--extra-bwa-options",
    nargs="*",
    help=(
        "additional bwa options to pass directly to bwa during alignments. "
        "Options must have '+'' in place of '-'."
    ),
)

# Parse command line arguments
args = vars(parser.parse_args())

# Process extra Freebayes options
extra_freebayes_options = args["extra_freebayes_options"]
if extra_freebayes_options is not None:
    extra_freebayes_options = [
        e.replace("+", "-") for e in extra_freebayes_options
    ]
else:
    extra_freebayes_options = []

# Read settings from settings file
settings_file = args["settings_file"]
settings = mip.get_analysis_settings(settings_file)

# Parse and set extra bwa options
ebo = args["extra_bwa_options"]
if ebo is not None:
    ebo = [[e.replace("+", "-") for e in ebo]]
    settings["bwaOptions"].extend(ebo)

# Set processor number setting
settings["processorNumber"] = args["threads"]

# Call Freebayes
mip.freebayes_call(
    bam_dir=args["bam_dir"],
    fastq_dir=args["fastq_dir"],
    options=extra_freebayes_options,
    vcf_file=args["output_vcf"],
    targets_file=args["targets_file"],
    make_fastq=args["skip_fastq"],
    align=args["skip_align"],
    settings=settings,
    bam_files=args["bam_files"],
    bam_list=args["bam_list"],
    fastq_padding=args["fastq_padding"],
    min_base_quality=args["min_base_quality"],
)
