import os, sys
from glob import glob
from tqdm import tqdm
import numpy as np
from collections import Counter
import argparse
from math import ceil
import random as rand
import multiprocessing as mp
from itertools import cycle
import pandas as pd

# This iteration of the haplotype module has had most of the arbitrary options removed
# Meaning that the shape and sampling frequency of simulations is dependent on the schema generated by inject_slim based on whatever data was used to generate the simulation
# PhysLen still needs to be supplied, and is found in the SLiM script used to generate data


class MsHandler:
    """Handles haplotype-tracked MS-style formatting from a standard SLiM output file.
    Runner function is parse_slim for easy tracking."""

    def __init__(self, mutfile):
        self.mutfile = mutfile

    def parse_slim(self, tol, physLen):
        """
        Runs all necessary steps to parse SLiM output and format it for hfs creation.

        Args:
            tol: Tolerance for window size
            physLen int: Length of chromosomes        
            
        Returns:
            list[str]: List of "lines" of a typical ms-format output. Used to be output as an intermediate file, but is now just passed as a list of str for parsing.
        """
        with open(self.mutfile, "r") as infile:
            lines = [i.strip() for i in infile.readlines()]

        cleaned_lines = self.remove_restarts(lines)
        mutations, genomes, samp_sizes, gens_sampled = self.readSampleOutFromSlimRun(
            cleaned_lines
        )
        newMutLocs = self.get_mutLocs(mutations, tol)
        unfilteredMuts = self.buildMutationPosMapping(newMutLocs, physLen)
        polyMuts = self.removeMonomorphic(unfilteredMuts, genomes)
        positionsStr = self.buildPositionsStr(polyMuts)

        # Iterate through timepoints, mutations is just a length indicator at this point
        segsitesStr = f"segsites: {len(polyMuts)}"
        haps = self.make_haps(polyMuts, genomes)
        out_ms = self.emitMsEntry(positionsStr, segsitesStr, haps)

        return out_ms, samp_sizes, gens_sampled

    def remove_restarts(self, lines):
        gens = []
        out_idxs = []
        for idx, line in enumerate(lines):
            if "#OUT" in line:
                gens.append(int(line.split(" ")[1]))
                out_idxs.append(idx)

        for gen in gens:
            gen_inds = [i for i, x in enumerate(gens) if x == gen]
            if (
                len(gen_inds) > 1
            ):  # Get indices of any duplicated gens - spans gens and out_idxs
                # Remove lines between the first and last occurence
                # Only want to keep the ones after the restart
                # Technically only restarts should happen at the dumpfile ggen, but this is flexible for anything I suppose
                del lines[out_idxs[gen_inds[0]] : out_idxs[gen_inds[-1]]]

        return lines

    def readSampleOutFromSlimRun(self, lines):
        """
        Adds genomes and mutations to object-wide dicts every time a new one is encountered.
        Scans through each line of SLiM output, at the end of a sample it will collate genomes and mutation info and store in dicts.
        """
        mode = 0
        mutations = {}
        genomes = []
        samp_sizes_list = []
        sampleText = []
        gens_sampled = []
        for idx, line in enumerate(lines):
            if mode == 0:
                if "#OUT" in line:
                    samp_sizes_list.append(int(line.split(" ")[4]))
                    gens_sampled.append(int(line.split(" ")[1]))
                    if idx != 0:
                        all_samp_genomes = self.addMutationsAndGenomesFromSample(
                            sampleText, mutations,
                        )
                        genomes.extend(all_samp_genomes)

                    sampleText = []

                else:
                    sampleText.append(line)

        # Last one
        all_samp_genomes = self.addMutationsAndGenomesFromSample(sampleText, mutations)
        genomes.extend(all_samp_genomes)

        return mutations, genomes, samp_sizes_list, gens_sampled

    def addMutationsAndGenomesFromSample(self, sampleText, mutations):
        """
        Maps mutation IDs to chromosomes that contain them, resulting in a genotype string that is added to the genomes list.

        Args:
            sampleText (list[str]): Lines of SLiM output relating to one timepoint sample from a series.
            mutations (dict[int]): Dict of mutation locations binned by their ID in the chromosome being sampled.
            genomes (list[set(int)]): List of genome IDs for each sample that have muts

        Mutations and Genomes are added in-scope.
        """
        mode = 0
        idMapping = {}
        genomes = []
        for line in sampleText:
            if mode == 0:
                if "Mutations" in line:
                    mode = 1
            elif mode == 1:
                if "Genomes" in line:
                    mode = 2
                elif len(line.strip().split()) == 9:
                    (
                        tempId,
                        permId,
                        mutType,
                        pos,
                        selCoeff,
                        domCoeff,
                        subpop,
                        gen,
                        numCopies,
                    ) = line.strip().split()
                    pos = int(pos)
                    if not pos in mutations:
                        mutations[pos] = {}
                    mutations[pos][permId] = 1
                    idMapping[tempId] = permId
            elif mode == 2:
                line = line.strip().split()
                gId, auto = line[:2]
                mutLs = line[2:]
                genomes.append(set([idMapping[x] for x in mutLs]))

        return genomes

    def get_mutLocs(self, mutations, tol):
        """
        Build new mutation map based on windows of mutations.

        Args:
            mutations (dict[int]): Dict of mutations binned by location.

        Returns:
            list[tuple(int, int)]: List of paired mutation (positions, IDs) in new format.
        """
        newMutLocs = []
        for mutPos in mutations:
            if len(mutations[mutPos]) == 1:
                mutId = list(mutations[mutPos].keys())[0]
                newMutLocs.append((mutPos, mutId))
            else:
                firstPos = mutPos - tol
                lastPos = mutPos + tol
                interval = (lastPos - firstPos) / (len(mutations[mutPos]) - 1)
                currPos = firstPos
                for mutId in mutations[mutPos]:
                    newMutLocs.append((currPos, mutId))
                    currPos += interval

        return newMutLocs

    def buildMutationPosMapping(self, mutLocs, physLen):
        """
        Creates new mapping relative to length of chromosome, adds to information tuple for mutation.

        Args:
            mutLocs list[tuple(int, int)]: List of paired mutation (positions, IDs) in new format.
            physLen int: Length of chromosomes
        Returns:
            list[tuple(int, int, float, int)]: Tuples of (newID, abs position, continuous position, permID).
        """
        mutMapping = []
        mutLocs.sort()
        for i in range(len(mutLocs)):
            pos, mutId = mutLocs[i]
            contPos = pos / physLen
            mutMapping.append((i, pos, contPos, mutId))

        return mutMapping

    def removeMonomorphic(self, allMuts, genomes):
        """
        Removes singletons by selecting only mutations that are polymorphic.

        Args:
            allMuts (list[tuple(int, int, float, int)]): Tuples of (newID, abs position, continuous position, permID)

        Returns:
            list[tuple(int, int, float, int)]: Tuples of (newID, abs position, continuous position, permID) for polymorphic mutations only.
        """
        newMuts = []
        newLocI = 0
        for locI, loc, contLoc, mutId in allMuts:
            freq = self.getFreq(mutId, genomes)
            if freq > 0 and freq < len(genomes):
                newMuts.append((newLocI, loc, contLoc, mutId))
                newLocI += 1

        return newMuts

    def getFreq(self, mut, genomes):
        """
        Calculate ocurrence of a mutation in each genome.

        Args:
            mut (tuple): Mutation information to query against the genome list.

        Returns:
            int: Number of times input mutation appears in all genomes.
        """
        mutCount = 0
        for genome in genomes:
            if mut in genome:
                mutCount += 1
        return mutCount

    def buildPositionsStr(self, muts):
        """
        Uses new mutation locations to build an MS-style mutation positions string.

        Args:
            muts (list[tuple]): Tuples of (newID, abs position, continuous position, permID) for each mutation.

        Returns:
            str: ms-style chromosome mutation positions string for use in downstream parsing.
        """
        positionsStr = []
        for _, locationDiscrete, _, mutId in muts:
            positionsStr.append(f"{locationDiscrete}.{mutId}")

        return "positions: " + " ".join(positionsStr)

    def make_haps(self, polyMuts, sampled_genomes):
        """
        Creates genotype 0/1 strings for each haplotype in ms-style format.

        Args:
            polyMuts (list[tuple]): Polymorphic mutations with ID, location, and permID fields.

        Returns:
            list[str]: All haplotype genotype strings for a given sample.
        """
        haps = []

        for i in range(len(sampled_genomes)):
            haps.append(["0"] * len(polyMuts))

        for i in range(len(sampled_genomes)):
            for locI, loc, contLoc, mutId in polyMuts:
                if mutId in sampled_genomes[i]:
                    haps[i][locI] = "1"

        return haps

    def emitMsEntry(self, positionsStr, segsitesStr, haps):
        """
        Writes a list of strings that is equivalent to the lines in an ms-formatted output.
        Can be edited to output to file instead easily.

        Args:
            positionsStr (str): Str of all positions with segsites
            segsitesStr (str): Str of number of segsites total in MS entry
            haps (list[str]]): All haplotype genotype strings for a given sample.

        Returns:
            List[str]: Expected MS output format for entire time series of sampled points and haps.
        """

        ms = []
        ms.append(f"slim {len(haps)} 1")
        ms.append("foo")
        ms.append("//")
        ms.append(segsitesStr)
        ms.append(positionsStr)
        for line in haps:
            ms.append("".join(line))

        return ms


class HapHandler:
    """
    Handles haplotype frequency spectrum generation, sorting, and output.
    Structure is very nested, user-facing function is readAndSplitMsData.
    """

    def __init__(self, hap_ms, maxSnps, samp_sizes, gens_sampled):
        self.hap_ms = hap_ms
        self.maxSnps = maxSnps
        self.samp_sizes = self.bin_samps(samp_sizes, gens_sampled)

    def bin_samps(self, samp_sizes, gens_sampled, gen_threshold=25, size_threshold=3):
        """
        Bins a list of ints into condensed bins where the minimum value is equal to <size_threshold>.
        Each bin must also not be larger than <gen_threshold>.

        Args:
            samp_sizes (list[int]): List of ints to bin.
            gens_sampled (list[int]): List of generations sampled.
            threshold (int, optional): Minimum value any given bin can be. Defaults to 3.

        Returns:
            list[int]: Binned values.
        """
        binned_sizes = []
        # print("Samp sizes:", len(samp_sizes))
        # print("Gens sampled:", len(gens_sampled))
        i = 0
        while i < len(samp_sizes):
            if samp_sizes[i] >= size_threshold:
                binned_sizes.append(samp_sizes[i])
                i += 1
            else:
                j = 0
                while sum(samp_sizes[i : i + j]) < size_threshold:
                    if gens_sampled[i + j] - gens_sampled[i] > gen_threshold:
                        # Good to go, append sample
                        break
                    elif (i + j) == len(gens_sampled) - 1:
                        # Hit the end before it's good, just take whatever's left
                        break
                    else:
                        # Need more samples, add the next timepoint
                        j += 1
                binned_sizes.append(sum(samp_sizes[i : i + j]))
                i += j

        return binned_sizes

    def readAndSplitMsData(self, inFileName):
        """Runner function that allows for broad exception catching from nested functions."""
        try:
            currTimeSeriesHFS = self.readMsData()
            X = np.array(currTimeSeriesHFS, dtype="float32")
            return X, "/".join([inFileName.split("/")[-3], inFileName.split("/")[-1]])

        except Exception as e:
            print(
                "couldn't make {} because of: {}".format(
                    inFileName.split("/")[-1].split(".")[0]
                ),
                e,
            )
            return None, None

    def readMsData(self):
        """
        Iterates through haplotype-tracked MS entry and creates haplotype matrices.


        Returns:
            list[list[float]]: Haplotype frequency spectrums for all timepoints; sorted by most common freq at any sampling point in series.
        """
        readMode = 0
        hapMats = []
        for line in self.hap_ms:
            if readMode == 0:
                if line.startswith("positions:"):
                    readMode = 1
                    currHaps = []
                elif line.startswith("segsites:"):
                    numSnps = int(line.strip().split()[-1])
                    if numSnps >= self.maxSnps:
                        start = int((numSnps - self.maxSnps) / 2)
                        end = start + self.maxSnps
                    else:
                        start, end = 0, numSnps
            elif readMode == 1:
                line = line.strip()
                if not line:
                    pass
                else:
                    if line[0] in ["0", "1"]:
                        currHaps.append(line[start:end])

        hapMats.append(self.getTimeSeriesHapFreqs(currHaps))

        return hapMats

    def getTimeSeriesHapFreqs(self, currHaps):
        """
        Build haplotype frequency spectrum for a single timepoint.

        Args:
            currHaps (list[str]): List of haplotypes read from MS entry.

        Returns:
            list[float]: Haplotype frequency spectrum for a single timepoint sorted by the most common hap in entire set.
        """
        winningHap = self.getMostCommonHapInEntireSeries(currHaps)
        hapBag = list(set(currHaps))
        hapBag.pop(hapBag.index(winningHap))

        hapToIndex = {}
        index = 0
        hapToIndex[winningHap] = index

        while len(hapBag) > 0:
            index += 1
            mostSimilarHapIndex = self.getMostSimilarHapIndex(hapBag, winningHap)
            mostSimilarHap = hapBag.pop(mostSimilarHapIndex)
            hapToIndex[mostSimilarHap] = index

        hapFreqMat = []
        i = 0
        for j in range(len(self.samp_sizes)):
            currHapFreqs = self.getHapFreqsForTimePoint(
                currHaps[i : i + self.samp_sizes[j]], hapToIndex, sum(self.samp_sizes),
            )
            hapFreqMat.append(currHapFreqs)
            i += self.samp_sizes[j]

        return hapFreqMat

    def getMostCommonHapInEntireSeries(self, haps):
        """
        Iterates through haplotype counts in entire series to find the most common hap.

        Args:
            haps (list[str]): List of haplotypes, each is a 1D genotype string

        Returns:
            str: Most frequent haplotype in entire sampling process
        """
        # Calculate haplotype frequency for each haplotype
        allFreqs = []
        i = 0
        j = 0
        for j in range(len(self.samp_sizes)):
            freqsInSamp = {}
            for hap in haps[i : i + self.samp_sizes[j]]:
                if not hap in freqsInSamp:
                    freqsInSamp[hap] = 0
                freqsInSamp[hap] += 1 / (i + self.samp_sizes[j])
            allFreqs.append(freqsInSamp)
            i += self.samp_sizes[j]

        # Calculate haplotype freq change from the start of sampling at each timestep and sort
        allHaps = []
        for timeStep in range(1, len(allFreqs)):
            for hap in allFreqs[timeStep]:
                if hap in allFreqs[0]:
                    freqChange = allFreqs[timeStep][hap] - allFreqs[0][hap]
                else:
                    freqChange = allFreqs[timeStep][hap]
                allHaps.append((allFreqs[timeStep][hap], freqChange, timeStep, hap))

        allHaps.sort()

        winningHapFreq, winningHapFreqChange, winningHapTime, winningHap = allHaps[-1]

        return winningHap

    def getMostSimilarHapIndex(self, haps, targHap):
        """
        Calculate distances between a current haplotype and all given haps in sample.

        Args:
            haps (list[str]): Haplotypes for a given sample point.
            targHap (str): Haplotype to calculate distance from.

        Returns:
            int: Index of the haplotype in the hapbag that has the min distance from targHap.
        """
        minDist = float("inf")
        for i in range(len(haps)):
            dist = self.seqDist(haps[i], targHap)
            if dist < minDist:
                minDist = dist
                minIndex = i

                return minIndex

    def getHapFreqsForTimePoint(self, currSample, hapToIndex, maxPossibleHaps):
        """
        Create haplotype freq spectrum for a given sample and haplotype.

        Args:
            currSample (list[str]): Set of haplotypes in current time-sample.
            hapToIndex (int): Index of hap from hap-bag to calculate with.
            maxPossibleHaps (int): Number of total possible haplotypes.

        Returns:
            list[float]: Haplotype frequency spectrum for a given set of haplotypes.
        """
        hfs = [0] * maxPossibleHaps
        for hap in currSample:
            hfs[hapToIndex[hap]] += 1

        hfs = [x / len(currSample) for x in hfs]

        return hfs

    def seqDist(self, hap1, hap2):
        """
        Calculates pairwise distance between two haplotypes

        Args:
            hap1 (list): Haplotype 1
            hap2 (list): Haplotype 2

        Returns:
            int: Number of pairwise differences (sequence distance) between haps
        """
        assert len(hap1) == len(hap2)
        numDiffs = 0
        for i in range(len(hap1)):
            if hap1[i] != hap2[i]:
                numDiffs += 1

        return numDiffs


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Haplotype frequency spectrum feature vector preparation.\
            Each run will result in an npz file named {samp_frequency}_{samp_size}.npz"
    )

    parser.add_argument(
        "-i",
        "--input-dir",
        metavar="INPUT_DIRECTORY",
        help="Base mutation type (hard/soft/etc) directory with *.pop files to create feature vectors from.",
        dest="in_dir",
        type=str,
        required=False,
    )

    parser.add_argument(
        "-p",
        "--phys-len",
        metavar="PHYS_LEN",
        help="Length of chromosome being simulated for training, will be specified in the stdpopsim call as well as in all SLiM scripts.",
        dest="physLen",
        default=12462531,
        type=int,
        required=False,
    )

    parser.add_argument(
        "--schema-name",
        metavar="SCHEMA-NAME",
        help="Name to use for output files. This is optional if running one instance, but necessary when doing multiple Snakemake runs at once.",
        dest="schema_name",
        type=str,
        required=True,
    )

    parser.add_argument(
        "--nthreads",
        metavar="NUM-PROCESSES",
        help="Number of threads available to multiprocessing module, more threads reduces runtime drastically. Defaults to all available - 1.",
        dest="nthreads",
        type=int,
        required=False,
        default=mp.cpu_count() - 1 or 1,
    )

    parser.add_argument(
        "-o",
        "--out-dir",
        metavar="OUT-DIR",
        help="Directory to write *.npz files to. Defaults to input dir.",
        dest="out_dir",
        type=str,
        required=False,
    )

    args = parser.parse_args()

    return args


def worker(args):
    mutfile, tol, physLen, maxSnps = args
    # try:
    # Handles MS parsing
    msh = MsHandler(mutfile, tol, physLen)
    hap_ms, samp_sizes, gens_sampled = msh.parse_slim(tol, physLen)
    # print("Sample sizes", samp_sizes)

    # Convert MS into haplotype freq spectrum and format output
    hh = HapHandler(hap_ms, maxSnps, samp_sizes, gens_sampled)
    X, id = hh.readAndSplitMsData(mutfile)
    #! (TPs * sampsize)
    X = np.squeeze(X)

    if X is not None and id is not None:
        return (id, X)
    # except Exception as e:
    #    print(e)
    #    pass


def main():
    argp = parse_arguments()

    print(f"Using {argp.nthreads} threads.")
    print("Data dir:", argp.in_dir)

    # Sanitize output dir
    if argp.out_dir is None:
        out_dir = argp.in_dir
    else:
        out_dir = argp.out_dir
    print("Output dir:", out_dir)

    filelist = glob(os.path.join(argp.in_dir, "*.pop"))
    physLen = argp.physLen
    tol = 0.5  # Allows for infinite sites model
    maxSnps = 50

    id_arrs = []

    args = zip(filelist, cycle([tol]), cycle([physLen]), cycle([maxSnps]))

    chunksize = 4
    pool = mp.Pool(processes=argp.nthreads)
    for proc_result in tqdm(
        pool.imap_unordered(worker, args, chunksize=chunksize),
        desc="Submitting processes...",
        total=len(filelist),
    ):
        id_arrs.append(proc_result)

    ids = []
    arrs = []
    # Have to do sanity check
    for i in id_arrs:
        if i:
            ids.append(i[0])
            arrs.append(i[1])

    print("Number of samples processed:", len(ids))
    print("Shape of single sample:", arrs[0].shape)

    np.savez(
        os.path.join(out_dir, f"hfs_{argp.schema_name}.npz"), **dict(zip(ids, arrs)),
    )
    print(
        "HFS data saved to:", os.path.join(out_dir, f"hfs_{argp.schema_name}.npz"),
    )


if __name__ == "__main__":
    main()
