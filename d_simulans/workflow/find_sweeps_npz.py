import argparse as ap
import logging
import os

from itertools import cycle
import numpy as np
import pandas as pd
import yaml
from tensorflow.keras.models import load_model
from tqdm import tqdm

from timesweeper.utils.frequency_increment_test import fit

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"

logging.basicConfig()
logger = logging.getLogger("timesweeper")
logger.setLevel("INFO")

pd.options.display.float_format = "{:.2f}".format


def load_npz(npz_path):
    npz_obj = np.load(npz_path)
    aft = npz_obj["aftIn"]  # (snps, timepoints, window)
    locs = npz_obj["aftInPosition"]  # (snps, window_locs)

    return aft, locs


def parse_npz_name(npz_path):
    splitpath = npz_path.split("_")
    chrom = splitpath[2]
    rep = splitpath[-1].split(".")[0]

    return chrom, rep


def run_aft_windows(ts_aft, locs, chrom, model):
    """
    Iterates through windows of MAF time-series matrix and predicts using NN.

    Args:
        snps (list[tup(chrom, pos,  mut)]): Tuples of information for each SNP. Contains mut only if benchmarking == True.
        genos (allel.GenotypeArray): Genotypes of all samples. 
        samp_sizes (list[int]): Number of chromosomes sampled at each timepoint.
        win_size (int): Number of SNPs to use for each prediction. Needs to match how NN was trained.
        model (Keras.model): Keras model to use for prediction.

    Returns:
        dict: Prediction values in the form of dict[snps[center]]
        np.arr: the central-most window, either based on mutation type or closest to half size of chrom.
    """
    left_edges = list(locs[:, 0])
    right_edges = list(locs[:, -1])
    centers = list(locs[:, 25])
    probs = model.predict(ts_aft)

    results_list = zip(cycle([chrom]), centers, left_edges, right_edges, probs)

    return results_list


def run_fit_windows(ts_aft, locs, chrom):
    """
    Iterates through windows of MAF time-series matrix and predicts using NN.

    Args:
        snps (list[tup(chrom, pos,  mut)]): Tuples of information for each SNP. Contains mut only if benchmarking == True.
        genos (allel.GenotypeArray): Genotypes of all samples. 
        samp_sizes (list[int]): Number of chromosomes sampled at each timepoint.
        win_size (int): Number of SNPs to use for each prediction. Needs to match how NN was trained.
        gens (list[int]): List of generations that were sampled.

    Returns:
        dict: P values from FIT.
    """
    gens = [i * 10 for i in range(7)]

    results_list = []
    for snp in tqdm(range(len(locs)), desc="Calculating FIT values"):
        results_list.append(
            (chrom, locs[snp, 25], fit(list(ts_aft[snp, :, 25]), gens)[1])
        )  # tval, pval

    return results_list


def load_nn(model_path, summary=False):
    """
    Loads the trained Keras network.

    Args:
        model_path (str): Path to Keras model.
        summary (bool, optional): Whether to print out model summary or not. Defaults to False.

    Returns:
        Keras.model: Trained Keras model to use for prediction.
    """
    model = load_model(model_path)
    if summary:
        print(model.summary())

    return model


def write_fit(fit_list, outfile):
    """
    Writes FIT predictions to file.

    Args:
        fit_dict (dict): FIT p values and SNP information.
        outfile (str): File to write results to.
    """
    chrom, bps, pval = zip(*fit_list)
    inv_pval = [1 - i for i in pval]

    predictions = pd.DataFrame({"Chrom": chrom, "BP": bps, "Inv_pval": inv_pval})
    # predictions = predictions[predictions["Inv_pval"] > 0.9]

    predictions.sort_values(["Chrom", "BP"], inplace=True)
    predictions.to_csv(
        os.path.join(outfile), header=True, index=False, sep="\t", float_format="%.3f"
    )


def write_preds(results_list, outfile, benchmark):
    """
    Writes NN predictions to file.

    Args:
        results_dict (dict): SNP NN prediction scores and window edges.
        outfile (str): File to write results to.
    """
    lab_dict = {0: "Neut", 1: "Soft"}
    chrom, centers, left_edges, right_edges, probs = zip(*results_list)

    neut_scores = [i[0] for i in probs]
    soft_scores = [i[1] for i in probs]
    classes = [lab_dict[np.argmax(i)] for i in probs]

    predictions = pd.DataFrame(
        {
            "Chrom": chrom,
            "BP": centers,
            "Class": classes,
            "Neut_Score": neut_scores,
            "Soft_Score": soft_scores,
            "Win_Start": left_edges,
            "Win_End": right_edges,
        }
    )

    # predictions = predictions[predictions["Neut_Score"] < 0.5]

    predictions.sort_values(["Chrom", "BP"], inplace=True)

    predictions.to_csv(outfile, header=True, index=False, sep="\t", float_format="%.3f")

    bed_df = predictions[["Chrom", "Win_Start", "Win_End", "BP"]]
    bed_df.to_csv(outfile.replace(".csv", ".bed"), header=False, index=False, sep="\t")


def add_file_label(filename, label):
    """Injects a model identifier to the outfile name."""
    splitfile = filename.split(".")
    newname = f"{''.join(splitfile[:-1])}_{label}.{splitfile[-1]}"
    return newname


def parse_ua():
    uap = ap.ArgumentParser(
        description="Module for iterating across windows in a time-series vcf file and predicting whether a sweep is present at each snp-centralized window."
    )
    uap.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="NPZ file with allele frequencies already processed in proper shape.",
        required=True,
    )
    uap.add_argument(
        "-o",
        "--out-dir",
        dest="outdir",
        help="Directory to write output to.",
        required=True,
    )
    uap.add_argument(
        "--benchmark",
        dest="benchmark",
        action="store_true",
        help="If testing on simulated data and would like to report the mutation \
            type stored by SLiM during outputVCFSample as well as neutral predictions, use this flag. \
            Otherwise the mutation type will not be looked for in the VCF entry nor reported with results.",
        required=False,
    )
    uap.add_argument(
        "--aft-model",
        dest="aft_model",
        help="Path to Keras2-style saved model to load for aft prediction.",
        required=True,
    )
    subparsers = uap.add_subparsers(dest="config_format")
    subparsers.required = True
    yaml_parser = subparsers.add_parser("yaml")
    yaml_parser.add_argument(
        metavar="YAML CONFIG",
        dest="yaml_file",
        help="YAML config file with all cli options defined.",
    )

    return uap.parse_args()


def read_config(yaml_file):
    """Reads in the YAML config file."""
    with open(yaml_file, "r") as infile:
        yamldata = yaml.safe_load(infile)

    return yamldata


def main(ua):
    outdir, aft_model = (
        ua.outdir,
        load_nn(ua.aft_model),
    )

    try:
        if not os.path.exists(outdir):
            os.makedirs(outdir)
    except:
        #running in high-parallel sometimes it errors when trying to check/create simultaneously
        pass

    # Load in everything
    logger.info(f"Loading data from {ua.input_file}")
    chrom, rep = parse_npz_name(ua.input_file)
    ts_aft, locs = load_npz(ua.input_file)

    # aft
    logger.info("Predicting with AFT")
    aft_predictions = run_aft_windows(ts_aft, locs, chrom, aft_model)
    write_preds(aft_predictions, f"{outdir}/aft_{chrom}_{rep}_preds.csv", ua.benchmark)
    logger.info(f"Done, results written to {outdir}/aft_{chrom}_{rep}_preds.csv")

    # FIT
    fit_predictions = run_fit_windows(ts_aft, locs, chrom)
    write_fit(fit_predictions, f"{outdir}/fit_{chrom}_{rep}_preds.csv")


if __name__ == "__main__":
    ua = parse_ua()
    main(ua)
