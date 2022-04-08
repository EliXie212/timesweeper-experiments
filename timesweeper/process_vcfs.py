import argparse
import logging
import multiprocessing as mp
import os
import subprocess
from glob import glob
from itertools import cycle
from find_sweeps import read_config

logging.basicConfig()
logger = logging.getLogger("vcf_processing")
logger.setLevel("INFO")


def read_multivcf(input_vcf):
    """Reads in file and returns as list of strings."""
    with open(input_vcf, "r") as input_file:
        raw_lines = [i.strip() for i in input_file.readlines()]

    return raw_lines


def split_multivcf(vcf_lines, header):
    """Splits the lines of multi-vcf file into list of vcf entries by <header> using itertools."""
    header_idxs = [i for i in range(len(vcf_lines)) if vcf_lines[i] == header]

    split_vcfs = []
    for idx in range(len(header_idxs[:-1])):
        split_vcfs.append(vcf_lines[header_idxs[idx] : header_idxs[idx + 1]])

    split_vcfs.append(vcf_lines[header_idxs[-1] :])

    return split_vcfs


def make_vcf_dir(input_vcf):
    """Creates directory named after vcf basename."""
    dirname = os.path.basename(input_vcf).split(".")[0]
    dirpath = os.path.dirname(input_vcf)
    vcf_dir = os.path.join(dirpath, dirname)
    os.makedirs(vcf_dir, exist_ok=True)

    return vcf_dir


def write_vcfs(vcf_lines, vcf_dir):
    """Writes list of vcf entries to numerically-sorted vcf files."""
    for idx, lines in enumerate(vcf_lines):
        with open(os.path.join(vcf_dir, f"{idx}.vcf"), "w") as outfile:
            outfile.writelines("\n".join(lines))


def index_vcf(vcf):
    cmd = f"""
    bgzip -f {vcf} > {vcf}.gz
    tabix -f -p vcf {vcf}.gz
    bcftools sort -Ov {vcf}.gz | bgzip -f > {vcf}.sorted.gz
    tabix -f -p vcf {vcf}.sorted.gz
    """
    subprocess.run(
        cmd, shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL
    )


def merge_vcfs(vcf_dir):
    num_files = len(glob(f"{vcf_dir}/*.vcf.sorted.gz"))
    cmd = f"""bcftools merge -Ov \
            --force-samples -0 \
            {" ".join([f"{vcf_dir}/{i}.vcf.sorted.gz" for i in range(num_files)])} > \
            {vcf_dir}/merged.vcf \
            """
    subprocess.run(cmd, shell=True)


def cleanup_intermed(vcf_dir):
    for ifile in glob(f"{vcf_dir}/*"):
        if "merged" not in ifile:
            os.remove(ifile)


def worker(input_vcf, num_tps, vcf_header):
    try:
        # Split into multiples after SLiM just concats to same file
        raw_lines = read_multivcf(input_vcf)
        split_lines = split_multivcf(raw_lines, vcf_header)

        split_lines = split_lines[len(split_lines) - num_tps :]

        # Creates subdir for each rep
        vcf_dir = make_vcf_dir(input_vcf)
        write_vcfs(split_lines, vcf_dir)

        # Now index and merge
        [index_vcf(vcf) for vcf in glob(f"{vcf_dir}/*.vcf")]
        merge_vcfs(vcf_dir)
        cleanup_intermed(vcf_dir)

    except Exception as e:
        print(e)
        pass


def get_ua():
    uap = argparse.ArgumentParser(
        description="Splits and re-merges VCF files to prepare for fast feature creation."
    )
    uap.add_argument(
        "--vcf-header",
        required=False,
        type=str,
        default="##fileformat=VCFv4.2",
        dest="vcf_header",
        help="String that tops VCF header, used to split entries to new files.",
    )
    uap.add_argument(
        "--threads",
        required=False,
        type=int,
        default=mp.cpu_count() - 1,
        dest="threads",
        help="Number of processes to parallelize across.",
    )
    subparsers = uap.add_subparsers(dest="config_format")
    subparsers.required = True
    yaml_parser = subparsers.add_parser("yaml")
    yaml_parser.add_argument(
        metavar="YAML CONFIG",
        dest="yaml_file",
        help="YAML config file with all cli options defined.",
    )

    cli_parser = subparsers.add_parser("cli")
    cli_parser.add_argument(
        "-w",
        "--work-dir",
        dest="work_dir",
        type=str,
        help="Directory used as work dir for simulate modules. Should contain simulated vcfs processed using process_vcf.py.",
        required=False,
        default=os.getcwd(),
    )
    cli_parser.add_argument(
        "--sample_sizes",
        required=True,
        type=int,
        nargs="+",
        dest="sample_sizes",
        help="Number of individuals to sample without replacement at each sampling point. Will be multiplied by ploidy to sample chroms from slim. Must match the number of entries in the -y flag.",
    )
    return uap.parse_args()


def main(ua):
    """Splits multivcf generated by SLiM using simulate.py and creates a merged/indexed vcf to create training data with."""
    if ua.config_format == "yaml":
        yaml_data = read_config(ua.yaml_file)
        work_dir, samp_sizes, threads, vcf_header = (
            yaml_data["work dir"],
            yaml_data["sample sizes"],
            ua.threads,
            ua.vcf_header,
        )

    elif ua.config_format == "cli":
        work_dir, samp_sizes, threads, vcf_header = (
            ua.work_dir,
            ua.sample_sizes,
            ua.threads,
            ua.vcf_header,
        )

    logger.info(f"Processing multiVCFs in {work_dir} using {threads} threads.")
    input_vcfs = glob(f"{work_dir}/vcfs/*/*.multivcf")
    pool = mp.Pool(threads)
    pool.starmap(
        worker,
        zip(input_vcfs, cycle([len(samp_sizes)]), cycle([vcf_header])),
        chunksize=5,
    )


if __name__ == "__main__":
    ua = get_ua()
    main(ua)