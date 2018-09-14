#! /usr/bin/env python

import argparse
import os
import sys
import subprocess
import re
import time
import multiprocessing
from packaging import version
from configobj import ConfigObj

try: # If not installed, this works
    import Trim
    import Collapse
    import ClipOverlap
    import __version
    import AdapterPredict
    import Call
except ImportError:
    from ProDuSe import Trim, Collapse, ClipOverlap, __version, AdapterPredict, Call


def isValidFile(file, parser, default=None):
    """
    Checks to ensure the provided file is exists, and throws an error if it is not.
    If a default value is provided, an argument matching the default will not throw
    an error even if the default is not a valid file

    :param file: A string containing a filepath to the file of interest
    :param parser: An argparse.ArgumentParser() object.
    :param default: The default argument value

    :returns: The "file" variable, if the file is valid
    :raises parser.error: If the file is not valid
    """
    if default is not None and default == file:
        return file
    elif os.path.exists(file):
        return file
    else:
        raise parser.error("Unable to locate \'%s\': No such file or directory" % (file))


def makeConfig(configName, configPath, arguments):
    """
    Creates a config file with the specified arguments
    :param configName: A string containing the base name of the config file
    :param configPath: A string containing a directory in which the config file will be placed
    :param arguments: A dictionary containing {argument: parameter} pairings. To be written to the config file
    :return:
    """
    config = ConfigObj()
    config.filename = configPath + os.sep + configName + "_task.ini"
    config[configName] = arguments
    config.write()


def configureOutput(sampleName, sampleParameters, outDir, argsToScript, printPrefix):
    # Create a sample-specific output directory inside the specified output directory
    # Inside this directory, we are going to create seperate directories for the data,
    # intermediate files, and results

    # If a directory corresponding to this sample already exists, don't proccess this sample,
    # as it likely has already been processed (at least partially)
    samplePath = outDir + os.sep + sampleName
    if os.path.exists(samplePath):
        sys.stderr.write("WARNING: A folder corresponding to \'%s\' already exists inside \'%s\'.\n" % (sampleName, outDir))
        sys.stderr.write("We will attempt to finish analyzing this sample using the existing configuration.\n")
        return samplePath

    os.mkdir(samplePath)

    # Create a tmp directory, which will hold intermediate files
    tmpDir = samplePath + os.sep + "tmp" + os.sep

    # Create a results directory which will hold the finalized BAM file, as well as the variant calls
    resultsDir = samplePath + os.sep + "results" + os.sep

    # Create a config directory, which will hold the config file used to run each step of the pipeline
    configDir = samplePath + os.sep + "config" + os.sep
    os.mkdir(configDir)

    # Group the arguments for this sample by the script they will be run in
    # I am going to specify defaults for bwa here, as there is no script to run BWA anymore
    scriptToArgs = {}

    bwaR1In = tmpDir + sampleName + ".trim_R1.fastq.gz"
    bwaR2In = tmpDir + sampleName + ".trim_R2.fastq.gz"
    bwaOut = tmpDir + sampleName + ".trim.bam"
    collapseOut = tmpDir + sampleName + ".collapse.bam"
    collapseSortedOut = resultsDir + sampleName + ".collapse.sort.bam"
    callPassedOut = resultsDir + sampleName + ".call.passed.vcf"
    callAllOut = resultsDir + sampleName + ".call.all.vcf"

    # Read group for BWA
    # Since the input and output of each script will be unique, specify them here
    scriptToArgs["bwa"] = {"input": [bwaR1In, bwaR2In], "output": bwaOut, "reference": sampleParameters["reference"]}
    scriptToArgs["trim"] = {"input": sampleParameters["fastqs"], "output": [bwaR1In, bwaR2In]}
    scriptToArgs["collapse"] = {"input": bwaOut, "output": collapseOut}
    scriptToArgs["call"] = {"input": collapseSortedOut, "output": callPassedOut, "unfiltered": callAllOut}

    for argument, scripts in argsToScript.items():

        for script in scripts:
            # Skip pipeline args, since those are no longer relevent
            if script == "pipeline":
                continue

            if script not in scriptToArgs:
                # Add the input and output mappings for this script
                scriptToArgs[script] = {}
            scriptToArgs[script][argument] = sampleParameters[argument]

    # Next, write a config file for each script
    for script, args in scriptToArgs.items():
        makeConfig(script, configDir, args)

    # Finally, lets create tmp directories (for the intermediate files), figure directories, and output directories
    os.mkdir(samplePath + os.sep + "tmp")
    os.mkdir(samplePath + os.sep + "results")
    os.mkdir(samplePath + os.sep + "figures")
    return samplePath


def combineArgs(confArgs, args):
    """
    Merges arguments specified in the config file with those parsed from the command line

    :param confArgs: A dictionary listing arguments parsed from the configuration file
    :param args: A Namespace object containing command line paramters
    """

    # Loop through each command line parameter, and if an argument was not provided, parse it from the
    # config file (assuming one was provided there)
    for parameter, arg in args.items():

        try:
            if arg is None and parameter in confArgs["Pipeline"]:
                # If this parameter is specified in the config file, add it to the argument groupings
                args[parameter] = confArgs["Pipeline"][parameter]
        except KeyError as e:
            raise TypeError("No section named \'Pipeline\' was found in the specified config file") from e

    return args


def checkArgs(rawArgs):
    """
    Re-parses arguments, and checks to ensure that they are both the required type, and that they exist if they
    are required

    This is mainly used to validate arguments from the configuration file

    :param rawArgs: A dictionary containing {argument: parameter}
    :returns: A namespace object coresponding to the validated parameters
    """

    # Convert the dictionary of arguments to a list, to allow for parsing
    listArgs = []
    for arg, parameter in rawArgs.items():
        if parameter is None or parameter is False or parameter == "None" or parameter == "False":
            continue
        listArgs.append("--" + arg)
        # If the parameter is a boolean, ignore it, as this will be reset once the arguments are re-parsed
        if parameter is True or parameter == "True":
            continue

        # Convert paramters that are lists into strings
        if isinstance(parameter, list):
            for p in parameter:
                listArgs.append(str(p))

        else:
            listArgs.append(str(parameter))

    # Point to the default filter included in the distribution
    try:  # If ProDuSe is installed
        produsePath = ProDuSe.__path__[0]
    except NameError:  # i.e. ProDuSe is not installed
        produsePath = os.path.dirname(os.path.realpath(__file__))
    defaultFilt = os.path.join(produsePath, "default_filter.pkl")

    # Did the user specify no_barcodes? if so, we don't need a lot of the arguments that the pipeline would normally
    # require, as we will set some defaults later
    if rawArgs["no_barcodes"] is True or rawArgs["no_barcodes"] == "True":
        barParamReq = False
    else:
        barParamReq = True

    # To validate the argument type, recreate the parser
    # BUT HERE, INDICATE IF AN ARGUMENT IS REQUIRED OR NOT
    parser = argparse.ArgumentParser(description="Runs all stages of the ProDuSe pipeline on the designated samples")
    parser.add_argument("-c", "--config", metavar="INI", default=None, type=lambda x: isValidFile(x, parser),
                        help="A configuration file, specifying one or more arguments. Overriden by command line parameters")
    parser.add_argument("-d", "--outdir", metavar="DIR", default="." + os.sep,
                        help="Output directory for analysis directory")
    parser.add_argument("-r", "--reference", metavar="FASTA", required=True,
                        help="Reference genome, in FASTA format. BWA indexes should be present in the same directory")
    parser.add_argument("--no_barcodes", action="store_true", help="Does not use semi-degenerate barcoded adapters")
    inputFiles = parser.add_mutually_exclusive_group(required=True)
    inputFiles.add_argument("--fastqs", metavar="FASTQ", default=None, nargs=2,
                            type=lambda x: isValidFile(x, parser), help="Two paired end FASTQ files")
    inputFiles.add_argument("-sc", "--sample_config", metavar="INI", default=None,
                            type=lambda x: isValidFile(x, parser),
                            help="A sample cofiguration file, specifying one or more samples")

    trimArgs = parser.add_argument_group("Arguments used when Trimming barcodes")
    trimArgs.add_argument("-b", "--barcode_sequence", metavar="NNNWSMRWSYWKMWWT", type=str,
                          help="The sequence of the degenerate barcode, represented in IUPAC bases")
    trimArgs.add_argument("-p", "--barcode_position", metavar="0001111111111110",  type=str, required=barParamReq,
                          help="Barcode positions to use when comparing expected and actual barcode sequences (1=Yes, 0=No)")
    trimArgs.add_argument("-mm", "--max_mismatch", metavar="INT", type=int, required=barParamReq,
                          help="The maximum number of mismatches permitted between the expected and actual barcode sequences [Default: 3]")
    trimArgs.add_argument("--trim_other_end", action="store_true",
                        help="In addition, examine the other end of the read for a barcode. Will not remove partial barcodes")

    collapseArgs = parser.add_argument_group("Arguments used when Collapsing families into a consensus")
    collapseArgs.add_argument("-fm", "--family_mask", metavar="0001111111111110", type=str, required=barParamReq,
                              help="Positions to consider when identifying reads are members of the same family. Usually the same as \'-b\'")
    collapseArgs.add_argument("-fmm", "--family_mismatch", metavar="INT", type=int, required=barParamReq,
                              help="Maximum number of mismatches allowed between two barcodes before they are considered as members of different families")
    collapseArgs.add_argument("-dm", "--duplex_mask", metavar="0000000001111110", type=str, required=barParamReq,
                              help="Positions to consider when determining if two families are in duplex")
    collapseArgs.add_argument("-dmm", "--duplex_mismatch", metavar="INT", type=int, required=barParamReq,
                              help="Maximum number of mismatches allowed between two barcodes before they are classified as not in duplex")
    collapseArgs.add_argument("-t", "--targets", metavar="BED", type=lambda x: isValidFile(x, parser),
                              help="A BED file containing capture regions of interest. Read pairs that do not overlap these regions will be filtered out")
    collapseArgs.add_argument("--tag_family_members", action="store_true",
                              help="Store the names of all reads used to generate a consensus in the tag 'Zm'")

    callArgs = parser.add_argument_group("Arguments used when calling variants")
    callArgs.add_argument("-f", "--filter", metavar="PICKLE", type=lambda x: isValidFile(x, parser), default=defaultFilt,
                          help="A python pickle containing a trained Random Forest variant filter")
    callArgs.add_argument("--threshold", metavar="FLOAT", type=float, default=0.33,
                          help="Classifier threshold to use when filtering variants. Decrease to be more lenient [Default: 0.33]")

    miscArgs = parser.add_argument_group("Miscellaneous Args")
    miscArgs.add_argument("-j", "--jobs", metavar="INT", type=int, default=1,
                        help="If multiple samples are specified (using \'-sc\'), how many samples will be processed in parallel")
    miscArgs.add_argument("--bwa", default="bwa", type=lambda x: isValidFile(x, parser, default="bwa"),
                        help="Path to bwa executable")
    miscArgs.add_argument("--samtools", default="samtools", type=lambda x: isValidFile(x, parser, default="samtools"),
                        help="Path to samtools executable")
    miscArgs.add_argument("--directory_name", default="produse_analysis_directory",
                        help="Default output directory name. The results of running the pipeline will be placed in this directory [Default: \'produse_analysis_directory\']")
    miscArgs.add_argument("--append_to_directory", action="store_true",
                        help="If \'--directory_name\' already exists in the specified output directory, simply append new results to that directory [Default: False]")
    miscArgs.add_argument("--cleanup", action="store_true", help="Remove intermediate files")

    validatedArgs = parser.parse_args(args=listArgs)

    return vars(validatedArgs)


def checkCommand(command, path, versionStr=None, minVer=None):
    """
    Ensures the command is installed on the system, and returns the version if it does exist

    If a command is not found, python will throw an OS error. Catch this, and inform the user.
    Otherwise, save and return the version number

    :param command: Literal name of the command
    :param path: A string containing a filepath to the command executable
    :param versionStr: The argument passed to the command to print out the version
    :param minVer: A minimum version number. If specified, the version number of the specified command will be checked
    :returns: The installed version of the specified command
    """

    try:
        runCom = [path]
        if versionStr:
            runCom.append(versionStr)
        runCheck = subprocess.Popen(runCom, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        versionRaw = []
        for line in runCheck.stdout:
            line = line.decode("utf-8")
            versionRaw.append(line)

        # Here, I am assuming that the above command 1) Does not include strings in the format of n.n before the
        # real version number, and 2) The version number (roughly) coresponds with standards (i.e. at least
        # major.minor )
        regex = re.compile('[0-9]+[.][.0-9]+')
        matches = re.findall(regex, "".join(versionRaw))

        # If there are no matches, then likely a version number is not generated by the above command.
        # In this case, simply warn the user, and continue
        if len(matches) == 0:
            sys.stderr.write("WARNING: Unable identify a version number for %s. Proceeding anyways...\n" % (command))
            currentVer = "0"
        else:
            currentVer = matches[0]
            # Check version number, if specified
            if minVer:

                if version.parse(minVer) > version.parse(currentVer):
                    sys.stderr.write(
                        "ERROR: The minimum version of %s required is %s, but the version provided is %s\n" % (
                            command, minVer, currentVer))
                    exit(1)

    # Prints out an error if the command is not installed
    except OSError:
        sys.stderr.write("ERROR: Unable to run \'%s\'\n" % (command))
        sys.stderr.write("Please ensure the executable exists, and try again\n")
        sys.exit(1)

    return currentVer


def createLogFile(filePath, args, **toolVer):
    """
    Creates a log file, specifying both the command line parameters used, as well as the version number of any tools

    :param filePath: A string listing the path of the log file. Any existing file with this name will be overwritten
    :param args: A dictionary storing argument: parameter values
    :param toolVer: A dictionary storing toolName: toolVersion

    :return: None
    """

    with open(filePath, "w") as o:

        # Add the version number of ProDuSe
        o.write("# ProDuSe Version " + __version.__version__ + os.linesep)
        o.write(os.linesep)

        # Add a status message explaining that this log file can actually be used as a configuration file
        o.write("# This log file can be used as a ProDuSe configuration file to repeat this analysis" + os.linesep)
        o.write("[Pipeline]")
        o.write(os.linesep)
        o.write(os.linesep)

        # Add command line parameters
        for argument, parameter in args.items():
            o.write(str(argument) + "=" + str(parameter) + os.linesep)
        o.write(os.linesep)

        # Add tool version numbers
        for tool, version in toolVer.items():
            o.write("# " + str(tool) + ": " + str(version.replace("\n", " ")) + os.linesep)
        o.write(os.linesep)


def checkIndexes(refFasta, bwaExec, samtoolsExec, outDir):
    """

    :param refFasta: A string containing a filepath to the reference genome
    :param bwaExec: A string containing a filepath to the bwa executable
    :param samtoolsExec: A string containing a filepath to the samtools executable
    :param outDir: A string containing a filepath to the base output directory
    :return: refPath: A string containing a filepath to the reference genome
    """

    hasBWAIndexes = True
    hasFastaIndex = True

    # Check to see if BWA indexes are present
    # These indexes are current as of 0.7.13-r1126
    indexSuffixes = (".amb", ".ann", ".bwt", ".sa")
    for suffix in indexSuffixes:
        if not os.path.exists(refFasta + suffix):
            hasBWAIndexes = False
            break  # If one index is missing, we will need to regenerate them all

    # Check to see if there is a fasta index
    if not os.path.exists(refFasta + ".fai"):
        hasFastaIndex = False

    # Now, check to see if the required indexes are present. If not, generate them
    if hasBWAIndexes and hasFastaIndex:  # All the required indexes are present
        return refFasta

    sys.stderr.write(
        "\t".join(["PRODUSE-MAIN\t", time.strftime('%X'), "Could not locate all indexes for \'%s\'\n" % refFasta]))
    sys.stderr.write("\t".join(["PRODUSE-MAIN\t", time.strftime('%X'), "Generating Indexes...\n"]))

    # Create a folder to hold the reference data
    refDir = os.path.join(outDir, "RefGenome")
    newRefFasta = os.path.join(refDir, os.path.basename(refFasta))
    os.mkdir(refDir)

    # Symlink the reference genome over, so we don't waste a bunch of disk space
    os.symlink(refFasta, newRefFasta)

    # If the BWA indexes exist, just symlink them over as well
    if hasBWAIndexes:
        for suffix in indexSuffixes:
            os.symlink(refFasta + suffix, refDir)
    else:
        # Generate the indexes
        bwaCommand = [bwaExec, "index", newRefFasta]
        bwaRun = subprocess.Popen(bwaCommand, stderr=subprocess.PIPE)
        for line in bwaRun.stderr:
            line = line.decode("utf-8")
            if "[BWTIncConstructFromPacked]" in line:
                continue
            elif "[main]" in line:  # Just cleaning up the terminal output a bit
                continue
            else:
                sys.stderr.write(line)

        bwaRun.wait()
        if bwaRun.returncode != 0:
            raise subprocess.CalledProcessError(" ".join(bwaCommand), bwaRun.returncode)

    if hasFastaIndex:
        # Some versions of BWA create this
        try:
            os.symlink(refFasta + ".fai", refDir)
        except FileExistsError:
            pass
    else:
        samtoolsCom = [samtoolsExec, "faidx", newRefFasta]
        subprocess.check_call(samtoolsCom)

    sys.stderr.write("\t".join(["PRODUSE-MAIN\t", time.strftime('%X'), "Indexes Generated. Using \'%s\' as the reference genome\n" % newRefFasta]))

    return newRefFasta


def runPipelineMultithread(args):
    """
    This function does basically nothing, except unpack the arguments provided and pass them to runPipeline

    This is required since "map" only passes each function one argument

    :param args: An iterable containing the arguments to be passed to runPipeline()
    :return:
    """

    runPipeline(args[0], args[1], args[2])

def runPipeline(sampleName, sampleDir, cleanup=False):
    """
    Run all scripts in the ProDuSe pipeline on the specified sample
    :param sampleName: A string containing the sample name, for status message updates
    :param sampleDir: A string containg the filepath to the base sample directory
    :param cleanup: A boolean indicating if temporary files should be deleted
    :return:
    """

    def runBWA(configPath):
        """
        Aligns the reads in the specified FASTQ files using the Burrows-Wheeler aligner
        In addition, a read group is added, and the resulting BAM file is sorted

        :param configPath: A string containing a filepath to a ini file listing bwa's parameters
        :return: None
        """

        sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Running BWA...\n"]))
        # Read the arguments from the config file
        try:
            bwaConfig = ConfigObj(configPath)["bwa"]
        except KeyError:  # Thrown if the section is not labelled "bwa"
            sys.stderr.write(
                "ERROR: The config file \'%s\' does not appear to be a bwa config, as no section is labelled \'bwa\'\n" % (
                configPath))
            exit(1)
        # Parse the arguments from the config file in the required order
        try:
            bwaCommand = [bwaConfig["bwa"], "mem", "-C",
                          bwaConfig["reference"],
                          bwaConfig["input"][0],
                          bwaConfig["input"][1],
                          ]
            sortCommand = [bwaConfig["samtools"],
                           "sort", "-O", "BAM", "-o",
                           bwaConfig["output"]]

            # To supress BWA's status messages, we are going to buffer the stderr stream of every process into a variable
            # If BWA or a samtools task crashes (exit code != 0), we will print out everything that is buffered
            bwaStderr = []
            samtoolsStderr = []
            bwaCounter = 0

            bwaCom = subprocess.Popen(bwaCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            sortCom = subprocess.Popen(sortCommand, stdin=bwaCom.stdout, stderr=subprocess.PIPE)

            # Parse through the stderr lines of BWA and samtools, and buffer them as necessary
            for bwaLine in bwaCom.stderr:
                # If this line indicates the progress of BWA, print it out
                bwaLine = bwaLine.decode("utf-8")
                if bwaLine.startswith("[M::mem_process_seqs]"):
                    bwaCounter += int(bwaLine.split(" ")[2])
                    sys.stderr.write(
                        "\t".join([printPrefix, time.strftime('%X'), "Reads Processed:" + str(bwaCounter) + "\n"]))
                bwaStderr.append(bwaLine)

            for samtoolsLine in sortCom.stderr:
                samtoolsStderr.append(samtoolsLine.decode("utf-8"))

            bwaCom.stdout.close()
            bwaCom.wait()
            sortCom.wait()

            if bwaCom.returncode != 0 or sortCom.returncode != 0:  # i.e. Something crashed
                sys.stderr.write("ERROR: BWA and Samtools encountered an unexpected error and were terminated\n")
                sys.stderr.write("BWA Standard Error Stream:\n")
                sys.stderr.write(os.linesep.join(bwaStderr))
                sys.stderr.write("Samtools Sort Standard Error Stream:\n")
                sys.stderr.write(os.linesep.join(samtoolsStderr))
                exit(1)

            sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Mapping Complete\n"]))

        except KeyError:  # i.e. A required argument is missing from the config file
            sys.stderr.write(
                "ERROR: Unable to locate a required argument in the bwa config file \'%s\'\n" % (configPath))
            exit(1)


    def sortAndRetag(inFile, outFile, bwaConfigPath):
        """
        Recalculare the MD and NM tags of the secified SAM file, and sort it

        :param inFile: A string containing a filepath to an input BAM file. Usually generated by clipOverlap
        :param outFile: A string containing an output filepath
        :param bwaConfigPath: A string containing a filepath to the BWA config file
        :return:
        """

        # Parse the reference genome location from the BWA config file
        bwaConfig = ConfigObj(bwaConfigPath)["bwa"]
        refGenome = bwaConfig["reference"]

        calmdCom = ["samtools", "calmd", inFile, refGenome, "-b"]  # Recalculate MD and NM tags
        sortCom = ["samtools", "sort", "-o", outFile]

        # To cleanup the terminal, we are going to buffer the stderr stream of every process into a variable
        # If a samtools task crashes (exit code != 0), we will print out everything that is buffered
        calmdStderr = []
        sortStderr = []

        calmdTask = subprocess.Popen(calmdCom, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        sortTask = subprocess.Popen(sortCom, stdin=calmdTask.stdout, stderr=subprocess.PIPE)

        # Parse through the stderr lines of samtools, and buffer them as necessary
        for calmdLine in calmdTask.stderr:
            calmdStderr.append(calmdLine.decode("utf-8"))

        for sortLine in sortTask.stderr:
            sortStderr.append(sortLine.decode("utf-8"))

        calmdTask.stdout.close()
        calmdTask.wait()
        sortTask.wait()

        if calmdTask.returncode != 0 or sortTask.returncode != 0:  # i.e. Something crashed
            sys.stderr.write("ERROR: Samtools encountered an unexpected error and was terminated\n")
            sys.stderr.write("Samtools calmd Standard Error Stream:\n")
            sys.stderr.write("\n".join(calmdStderr))
            sys.stderr.write("Samtools Sort Standard Error Stream:\n")
            sys.stderr.write("\n".join(sortStderr))
            exit(1)

        # Finally, index the BAM file
        subprocess.check_call(["samtools", "index", outFile])

    printPrefix = "PRODUSE-MAIN\t\t"+ sampleName
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Processing Sample \'%s\'\n" % sampleName.rstrip()]))

    # Run Trim
    trimDone = os.path.join(sampleDir, "config", "Trim_Complete")  # Similar to Make's "TASK_COMPLETE" file
    if not os.path.exists(trimDone):  # i.e. Did Trim already complete for this sample? If so, do not re-run it
        trimConfig = os.path.join(sampleDir, "config", "trim_task.ini")  # Where is Trim's config file?
        trimPrintPrefix = "PRODUSE-TRIM\t\t" + sampleName
        Trim.main(sysStdin=["--config", trimConfig], printPrefix=trimPrintPrefix)  # Actually run Trim
        open(trimDone,
             "w").close()  # After Trim completes, it will create this file, signifying to the end user that this task completed

    # Run bwa
    bwaDone = os.path.join(sampleDir, "config", "BWA_Complete")
    if not os.path.exists(bwaDone):
        bwaConfig = os.path.join(sampleDir, "config", "bwa_task.ini")
        runBWA(bwaConfig)
        open(bwaDone, "w").close()

    # Run Collapse
    collapseDone = os.path.join(sampleDir, "config", "Collapse_Complete")
    if not os.path.exists(collapseDone):
        collapseConfig = os.path.join(sampleDir, "config", "collapse_task.ini")
        collapsePrintPrefix = "PRODUSE-COLLAPSE\t" + sampleName
        Collapse.main(sysStdin=["--config", collapseConfig], printPrefix=collapsePrintPrefix)

        # Sort the collapse BAM output output
        # Parse the config file for the output file name
        collapseConfArgs = ConfigObj(collapseConfig)
        sortInput = collapseConfArgs["collapse"]["output"]
        # Append "sort" as the output file name
        sortOutput = sortInput.replace(".bam", ".sort.bam")
        tmpDir = os.sep + "tmp" + os.sep
        resultsDir = os.sep + "results" + os.sep
        sortOutput = sortOutput.replace(tmpDir, resultsDir)
        sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Sorting final BAM file, and recalculating tags...\n"]))
        bwaConfig = os.path.join(sampleDir, "config", "bwa_task.ini")
        sortAndRetag(sortInput, sortOutput, bwaConfig)

        open(collapseDone, "w").close()

    # Run call (variant calling)
    callDone = os.path.join(sampleDir, "config", "Call_Complete")
    pipelineDone = os.path.join(sampleDir, "config", "Pipeline_Complete")
    if not os.path.exists(callDone):
        callConfig = os.path.join(sampleDir, "config", "call_task.ini")
        callPrintPrefix = "PRODUSE-CALL\t\t" + sampleName
        Call.main(sysStdin=["--config", callConfig], printPrefix=callPrintPrefix)
        open(callDone, "w").close()

        # Mark this sample as fully processed
        open(pipelineDone, "w").close()

    # Cleanup intermediate files (if specified)
    if cleanup:
        tmpDir = os.path.join(sampleDir, "tmp")
        tmpFiles = os.listdir(tmpDir)
        for tmpFile in tmpFiles:
            os.remove(os.path.join(tmpDir, tmpFile))

    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "%s: Pipeline Complete\n" % sampleName.rstrip()]))


parser = argparse.ArgumentParser(description="Runs all stages of the ProDuSe pipeline on the designated samples")
parser.add_argument("-c", "--config", metavar="INI", default=None, type=lambda x: isValidFile(x, parser),
                    help="A configuration file, specifying one or more arguments. Overriden by command line parameters")
parser.add_argument("--fastqs", metavar="FASTQ", default=None, nargs=2, type=lambda x: isValidFile(x, parser),
                    help="Two paired end FASTQ files")
parser.add_argument("-sc", "--sample_config", metavar="INI", default=None, type=lambda x: isValidFile(x, parser),
                    help="A sample cofiguration file, specifying one or more samples")
parser.add_argument("-r", "--reference", metavar="FASTA",
                    help="Reference genome, in FASTA format. BWA indexes should present in the same directory")
parser.add_argument("-d", "--outdir", metavar="DIR", help="Output directory for analysis directory")
parser.add_argument("--no_barcodes", action="store_true", help="Does not use semi-degenerate barcoded adapters")

trimArgs = parser.add_argument_group("Arguments used when Trimming barcodes")
trimArgs.add_argument("-b", "--barcode_sequence", metavar="NNNWSMRWSYWKMWWT", type=str,
                      help="The sequence of the degenerate barcode, represented in IUPAC bases")
trimArgs.add_argument("-p", "--barcode_position", metavar="0001111111111110", type=str,
                      help="Barcode positions to use when comparing expected and actual barcode sequences (1=Yes, 0=No)")
trimArgs.add_argument("-mm", "--max_mismatch", metavar="INT", type=int,
                      help="The maximum number of mismatches permitted between the expected and actual barcode sequences [Default: 3]")
trimArgs.add_argument("--trim_other_end", action="store_true", help="In addition, examine the other end of the read for a barcode. Will not remove partial barcodes")

collapseArgs = parser.add_argument_group("Arguments used when Collapsing families into a consensus")
collapseArgs.add_argument("-fm", "--family_mask", metavar="0001111111111110", type=str,
                          help="Positions to consider when identifying reads are members of the same family. Usually the same as \'-b\'")
collapseArgs.add_argument("-fmm", "--family_mismatch", metavar="INT", type=int,
                          help="Maximum number of mismatches allowed between two barcodes before they are considered as members of different families")
collapseArgs.add_argument("-dm", "--duplex_mask", metavar="0000000001111110", type=str,
                          help="Positions to consider when determining if two families are in duplex")
collapseArgs.add_argument("-dmm", "--duplex_mismatch", metavar="INT", type=int,
                          help="Maximum number of mismatches allowed between two barcodes before they are classified as not in duplex")
collapseArgs.add_argument("-t", "--targets", metavar="BED", type=lambda x: isValidFile(x, parser),
                          help="A BED file containing capture regions of interest. Read pairs that do not overlap these regions will be filtered out")
collapseArgs.add_argument("--tag_family_members", action="store_true",
                          help="Store the names of all reads used to generate a consensus in the tag 'Zm'")

callArgs = parser.add_argument_group("Arguments used when calling variants")
callArgs.add_argument("--threshold", metavar="FLOAT", type=float, help="Classifier threshold to use when filtering variants. Decrease to be more lenient [Default: 0.33]")
callArgs.add_argument("-f", "--filter", metavar="PICKLE", type=lambda x:isValidFile(x, parser), help="A python pickle containing a trained Random Forest variant filter")

miscArgs = parser.add_argument_group("Miscellaneous Arguments")
miscArgs.add_argument("-j", "--jobs", metavar="INT", type=int, help="If multiple samples are specified (using \'-sc\'), how many samples will be processed in parallel")
miscArgs.add_argument("--bwa", help="Path to bwa executable [Default: \'bwa\']")
miscArgs.add_argument("--samtools", help="Path to samtools executable [Default: \'samtools\']")
miscArgs.add_argument("--directory_name",
                    help="Default output directory name. The results of running the pipeline will be placed in this directory [Default: \'produse_analysis_directory\']")
miscArgs.add_argument("--append_to_directory", action="store_true",
                    help="If \'--directory_name\' already exists in the specified output directory, simply append new results to that directory [Default: False]")
miscArgs.add_argument("--cleanup", action="store_true", help="Remove intermediate files")

# For config parsing purposes, assign each parameter to the pipeline component from which it originates
argsToPipelineComponent = {
    "fastqs": ["pipeline"],
    "barcode_sequence": ["trim"],
    "barcode_position": ["trim"],
    "max_mismatch": ["trim"],
    "sample_config": ["pipeline"],
    "reference": ["pipeline", "collapse", "call"],
    "bwa": ["bwa"],
    "samtools": ["bwa"],
    "family_mask": ["collapse"],
    "family_mismatch": ["collapse"],
    "duplex_mask": ["collapse"],
    "duplex_mismatch": ["collapse"],
    "targets": ["collapse", "call"],
    "tag_family_members": ["collapse"],
    "filter" : ["call"],
    "threshold": ["call"]
}


def main(args=None, sysStdin=None):
    if args is None:
        if sysStdin is None:
            args = parser.parse_args()
            args = vars(args)
        else:
            args = parser.parse_args(sysStdin)
            args = vars(args)

    # If a config file was specified, parse the arguments out of it
    confArgs = None
    if args["config"] is not None:
        confArgs = ConfigObj(args["config"])
        args = combineArgs(confArgs, args)

    # Next, ensure that required parameters were provided, and they are of the correct type
    # This is done here so that parameters passed from the config file can also be checked
    args = checkArgs(args)

    # Next, lets make sure that ProDuSe's dependencies exist, and they meet the minumum version
    # requirements
    bwaVer = checkCommand("bwa", args["bwa"], minVer="0.7.0")
    samtoolsVer = checkCommand("samtools", args["samtools"], minVer="1.3.1")
    pythonVer = sys.version

    printPrefix = "PRODUSE-MAIN\t"
    sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Starting..." + "\n"]))

    # Next, let's configure a base output directory for the analysis
    # First, check to see if this output directory already exists in the specified path
    args["outdir"] = os.path.abspath(args["outdir"])
    baseOutDir = args["outdir"] + os.sep + args["directory_name"]
    if os.path.exists(baseOutDir):
        if not args["append_to_directory"]:
            sys.stderr.write("ERROR: \'%s\' already exists in \'%s\'\n" % (args["directory_name"], args["outdir"]))
            sys.stderr.write("Use \'--append_to_directory\' to append results to the existing directory\n")
            sys.exit(1)
    else:
        os.mkdir(baseOutDir)

    # Create a log file specifying the command line parameters
    logFileName = baseOutDir + os.sep + "ProDuSe_task.log"
    createLogFile(logFileName, args, bwa=bwaVer, samtools=samtoolsVer, python=pythonVer)

    # Finally, organize samples and prepare to run the entire pipeline on each sample
    # If only a single sample was specified at the command (i.e. a single pair of FASTQ files, using -f)
    samples = {}
    if args["fastqs"]:
        # Generate a sample name
        baseName = os.path.basename(args["fastqs"][0]).split(".")[0]
        samples[baseName] = {}
    # Otherwise, we need to parse each sample from the sample config file
    # along with sample-specific arguments
    else:
        sConfig = ConfigObj(args["sample_config"])
        for sample, arguments in sConfig.items():
            samples[sample] = arguments

    # Check for the required BWA indexes. If they do not exist, we will generate them
    args["reference"] = checkIndexes(args["reference"], args["bwa"], args["samtools"],
                                     os.path.join(args["outdir"], args["directory_name"]))

    # Setup each sample
    samplesToProcess = {}
    barcodeWarned = False
    for sample, sampleArgs in samples.items():
        # Override the existing arguments with any sample-specific arguments
        runArgs = args.copy()
        if len(sampleArgs) != 0:
            for argument, parameter in sampleArgs.items():
                runArgs[argument] = parameter

        # Double check to ensure that --fastqs were provided, since we have assumed that they would be specified in the
        # sample config file
        if "fastqs" not in runArgs:
            sys.stderr.write("\t".join(
                [printPrefix, time.strftime('%X'), "WARNING: No FASTQ files were provided for \'%s\'\n" % (sample)]))
            sys.stderr.write(
                "Ensure \'fastqs = <fq1.fastq.gz> <fq2.fastq.gz>\' is specified in the sample config file. Skipping...\n")
            continue
        elif len(runArgs["fastqs"]) != 2:
            sys.stderr.write("\t".join(
                [printPrefix, time.strftime('%X'),
                 "WARNING: Two FASTQ files must be provided for \'%s\', Skipping...\n" % (sample)]))
            continue

        # If a barcode was not specified, estimate it using adapter_predict
        # Note that this will end catastrophically if the sample is multiplexed
        if runArgs["barcode_sequence"] is None:
            if not barcodeWarned:
                sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Predicting Barcode Sequences...\n"]))
                barcodeWarned = True
            adapterPredictArgs = ["--max_barcode_length", str(len(runArgs["barcode_position"])), "--input"]
            adapterPredictArgs.extend(runArgs["fastqs"])
            runArgs["barcode_sequence"] = AdapterPredict.main(sysStdin=adapterPredictArgs, supressOutput=True)
            # Check if the resulting barcode is garbage
            if set(runArgs["barcode_sequence"]) == {"N"}:
                sys.stderr.write("WARNING: Unable to predict barcode sequence for \'%s\'. Skipping..." % runArgs["sample"])
                continue
            elif "N" in runArgs["barcode_sequence"]:  # TODO: Improve this estimate
                sys.stderr.write("WARNING: The barcode sequence predicted for \'%s\' is quite broad\n" % sample)
                sys.stderr.write("We'll continue anyways, but if these FASTQs are multiplexed, the resulting analysis will merge those samples\n")

        # Configure the output directories and config files for this sample
        sampleDir = configureOutput(sample, runArgs, baseOutDir, argsToPipelineComponent, printPrefix)
        # Create a sample-specific log file
        sLogName = os.path.join(sampleDir, "config", sample + "_Task.log")
        createLogFile(sLogName, runArgs, bwa=bwaVer, samtools=samtoolsVer, python=pythonVer)
        samplesToProcess[sample] = sampleDir

    # If multiple samples were provided, and more than 1 job was specified, run these samples in parallel
    if args["jobs"] == 0:  # i.e. use as many threads as possible
        threads = None  # The multiprocessing module will chose the number of threads
    else:
        threads = min(len(samplesToProcess), args["jobs"], os.cpu_count())  # Just as a safety measure, don't let the user spawn more threads then there are CPU cores

    if threads > 1 or threads is None:
        processPool = multiprocessing.Pool(processes=threads)

        # Convert the dictionary which contains the sample arguments into a list of tuples, for multithreading use
        # In addition, to keep the command line status messages semi-reasonable, normalize for sample name length
        maxLength = max(len(x) for x in samplesToProcess.keys())
        multithreadArgs = list((x.ljust(maxLength, " "), y, args["cleanup"]) for x, y in samplesToProcess.items())
        try:
            # Run the jobs
            processPool.map_async(runPipelineMultithread, multithreadArgs)
            processPool.close()
            processPool.join()
        except KeyboardInterrupt as e:
            processPool.terminate()
            processPool.join()
            raise e

    else:
        # Run each sample in series
        for sample, sampleDir in samplesToProcess.items():
            runPipeline(sample, sampleDir, args["cleanup"])


if __name__ == "__main__":
    main()
