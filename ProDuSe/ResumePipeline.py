#!/usr/bin/env python

"""
Resumes a previously terminated analysis. Only works for analysis generated using ProDuSe version 0.9 and above
"""

import argparse
import os
import sys
try:
    import ProdusePipeline
except ImportError:
    from ProDuSe import ProdusePipeline


def isValidDir(dir, parser):
    """
    Checks to ensure that the specified directory exists

    :param dir: A string containing a filepath to the directory
    :param parser: An argparse.ArgumentParser()object, or equivelent
    :return: dir: The same string
    :raises: parser.error: If the specified directory does not exist
    """

    if not os.path.isdir(dir):
        raise parser.error("ERROR: Unable to locate directory %s" % dir)
    else:
        return dir

def getArgs(stdin):

    parser = argparse.ArgumentParser(description="Resumes analysis of a previously terminated Pipeline")
    parser.add_argument("-d", "--produse_dir", required=True, type=lambda x: isValidDir(x, parser), help="An existing output directory for ProDuSe analysis")
    if stdin is None:
        return parser.parse_args()
    else:
        return parser.parse_args(stdin)


def main(args=None, sysStdin=None):
    if args is None:
        args = getArgs(sysStdin)

    args = vars(args)
    args["produse_dir"] = os.path.abspath(args["produse_dir"])
    # Assuming that the directory the user provided is actually a "produse_analysis_directory", we need to figure out what samples are in it
    # Since each directory corresponds to a different sample, obtain a list of all subdirectories
    samples = list(os.path.join(args["produse_dir"], x) for x in os.listdir(path = args["produse_dir"]) if os.path.isdir(os.path.join(args["produse_dir"], x)))

    # To actually check that this is a "produse_analysis_directory", check for the required config files for each sample
    validSamples = []
    requiredConfigs = ["trim_task.ini", "bwa_task.ini", "collapse_task.ini", "clipoverlap_task.ini"]
    for sample in samples:
        validSample = True
        for config in requiredConfigs:
            configPath = os.path.join(sample, "config", config)
            if not os.path.exists(configPath):
                validSample = False
                break

        if validSample:
            validSamples.append(sample)

    print(validSamples)

    # If no valid samples were found, then this is not a valid analysis directory
    if len(validSamples) == 0:
         sys.stderr.write("ERROR: Unable to find any valid sample directories in \"%s\". Check that this a directory created by \"ProDuSe Pipeline\", and that setup completed sucessfully" % args["produse_dir"] + os.linesep)

    # Obtain the sample name, and re-run this sample
    for sample in validSamples:
        sampleName = os.path.basename(sample)
        ProdusePipeline.runPipeline(sampleName, sample)


if __name__ == "__main__":
    main()
