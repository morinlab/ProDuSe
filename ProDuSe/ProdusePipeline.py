#! /usr/bin/env python

# USAGE:
# 	See ProdusePipeline -h for details
#
# DESCRIPTION
# 	Runs the entire ProDuSe variant calling pipeline on the provided samples
# 	A single sample can be run by specifying -f <fastq1> <fastq2>
# 	To run multiple samples, specify a sample_config file
# 	Note that paramters specified in the sample_config file override other paramters for that sample
# 	Arguments specified in the produse_config file override command line arguments
#
# AUTHORS
# 	Christopher Rushton (ckrushto@sfu.ca)

# Import standard python modules
import configargparse
import os
import subprocess
import sys
import time

# Import ProDuSe modules
# If not installed or running in python2, this works fine
try:
	import configure_produse
	import bwa
	import trim
	import collapse
	import snv
	import filter_produse
	import SplitMerge
	import __version as ProDuSeVer

# If installed and running in python3
except ImportError:
	from ProDuSe import configure_produse
	from ProDuSe import bwa
	from ProDuSe import trim
	from ProDuSe import collapse
	from ProDuSe import snv
	from ProDuSe import filter_produse
	from ProDuSe import SplitMerge
	from ProDuSe import __version as ProDuSeVer
"""
Processes command line arguments

"""

# Look, we both know this is terrible. But resolve_conflicts is also terrible. Thus, I don't have any choice but to copy-paste arguments over
parser = configargparse.ArgumentParser(description="Runs the entire ProDuSe pipeline on the supplied samples. ")

# Universal Args
pipelineArgs = parser.add_argument_group("Pipeline Arguments")

pipelineArgs.add_argument("-c", "--config", required=False, is_config_file=True, help="ProDuSe config file, listing adapter sequences, positions, and other parameters to be passed to pipeline scripts. See (See \'etc/produse_config.ini\' for an example)")
pipelineArgs.add_argument("-x", "--stitcherpath", required=True, type=lambda x: isValidFile(x, parser), help="Path to Illumina's Stitcher.exe (Can be obtained from \'https://github.com/Illumina/Pisces\')")

# configure_produse args
confArgs = parser.add_argument_group("configure_produse Arguments")
inputArgs = confArgs.add_mutually_exclusive_group(required=True)
inputArgs.add_argument(
    "-f", "--fastqs",
    metavar="FASTA",
    nargs=2,
    type=lambda x: isValidFile(x, parser),
    help="Forward and reverse fastq files from paired end sequencing. Overwritten by --sample_config"
    )
confArgs.add_argument(
    "-d", "--output_directory", "-o",
    metavar="DIR",
    default="." + os.sep,
    help="Output directory for ProDuSe analysis [Default: %(default)s]"
    )
confArgs.add_argument(
    "-r", "--reference",
    metavar="FASTA",
    required=True,
    type=lambda x: isValidFile(x, parser),
    help="Reference genome, in FASTA format. A BWA Index should be located in the same directory"
    )

inputArgs.add_argument(
    "-sc", "--sample_config",
    metavar="INI",
    required=False,
    type=lambda x: isValidFile(x, parser),
    help="Config file listing sample names and FASTQ locations (See \'etc/sample_config.ini\' for an example)"
    )

confArgs.add_argument(
	"-n", "--normal",
	type=lambda x:isValidFile(x, parser),
	metavar="BAM/FASTQS",
	nargs="+",
	help="A BAM or paired FASTQ files coresponding to a matched normal"
	)
confArgs.add_argument(
	"-ns",
	"--normal_adapter_sequence",
	metavar="NNNWSMRWSYWKMWWT",
	help="The matched normal adapter sequence, if barcoded adapters were used"
	)
confArgs.add_argument(
	"-np",
	"--normal_adapter_position",
	metavar="0001111111111111",
	help="The matched normal adapter position, if barcoded adapters were used"
	)
confArgs.add_argument(
	"--append_to_directory",
	action="store_true",
	help="Place results into an existing \'produse_analysis_directory\'. Note that samples with conflicting names will be skipped"
)
# Trim args
trimArgs = parser.add_argument_group("Trim Arguments")
trimArgs.add(
    "--adapter_sequence",
    type=str,
    required=True,
    help="The randomized adapter sequence flanked in input fastq files described using IUPAC bases"
    )
trimArgs.add(
    "--adapter_position",
    type=str,
    required=True,
    help="The positions in the adapter sequence to include in distance calculations, 0 for no, 1 for yes"
    )
trimArgs.add(
    "--max_mismatch",
    type=int,
    required=True,
    help="The maximum number of mismatches allowed between the expected and actual adapter sequences",
    )
trimArgs.add(
    "-v",
    action="store_true",
    help="Instead, output entries that are distant from the adapter sequence"
    )
trimArgs.add(
    "-u",
    action="store_true",
    help="Instead, output entries without trimming the adapter sequence"
    )

# bwa args
bwaArgs = parser.add_argument_group("bwa Arguments")
bwaArgs.add(
    "-t", "--threads",
    required=False,
    default=1,
    help="Number of threads to use while running BWA [Default: %(default)s]"
    )

bwaArgs.add(
    "-R", "--readgroup",
    required=False,
    type=str,
    help="Fastq read group"
    )

# Collapse arguments
collapseArgs = parser.add_argument_group("Collapse Arguments")
collapseArgs.add(
    "-sp", "--strand_position",
    metavar="STR",
    type=str,
    required=True,
    help="The positions in the adapter sequence to include in distance calculations, 0 for no, 1 for yes."
    )

collapseArgs.add(
    "-dp", "--duplex_position",
    metavar="STR",
    type=str,
    required=True,
    help="The positions in the adapter sequence to include in distance calculations between forward and reverse reads, 0 for no, 1 for yes"
    )
collapseArgs.add(
    "-fp", "--family_plot",
    type=str,
    required=False,
    default=None,
    help="A histogram to plot molecule counts per read family (i.e. each consensus read)"
    )

# Used to maintain backwards compatibility with the poorly-named strand mis-match
adapterMismatch = collapseArgs.add_mutually_exclusive_group(required=True)
adapterMismatch.add(
    "-amm", "--adapter_max_mismatch",
    type=int,
    help="The maximum number of mismatches allowed between the expected and actual adapter sequences",
    )
adapterMismatch.add(
    "--strand_max_mismatch",
    type=int,
    help=configargparse.SUPPRESS,
)

collapseArgs.add(
    "-dmm", "--duplex_max_mismatch",
    type=int,
    required=True,
    help="The maximum number of mismatches allowed between the expected and actual duplexed adapter sequences",
    )

collapseArgs.add(
    "-smm", "--sequence_max_mismatch",
    type=int,
    required=False,
    default=20,
    help="The maximum number of mismatches allowed in an alignment"
    )

collapseArgs.add(
    "-oo", "--original_output",
    type=str,
    required=False,
    action="append",
    default=None,
    help="A pair of empty fastq files to rewrite original fastqs with duplex information"
    )

# SNV arguments
snvArgs = parser.add_argument_group("SNV Argumentss")

snvArgs.add_argument(
    "-tb", "--target_bed",
    required=False,
    help="A tab-delinated file listing regions on which variant calling will be restricted to"
    )
snvArgs.add_argument(
    "-vaft", "--variant_allele_fraction_threshold",
    default=0.01,
    type=float,
    help="Minimum variant frequency threshold for each strand [Default: %(default)s]"
    )
snvArgs.add_argument(
    "-mo", "--min_molecules",
    default=40,
    type=int,
    help="Number of total molecules (supporting or otherwise) required to call a variant at this position. Reduce this if you are running only on positions you expect to be mutated [Default: %(default)s]"
    )
snvArgs.add_argument(
    "-mum", "--mutant_molecules",
    default=3,
    required=False,
    type=int,
    help="Number of TOTAL molecules (i.e. on the forward and reverse strand) required to call a variant as real (set to 0 if you are forcing variant calling at known sites) [Default: %(default)s]"
    )
snvArgs.add_argument(
    "-mrpu", "--min_reads_per_uid",
    default=2,
    type=int,
    help="Bases with support between MRPU and SSBT will be classified as a weak supported base [Default: %(default)s]"
    )
snvArgs.add_argument(
    "-ssbt", "--strong_supported_base_threshold",
    default=3,
    type=int,
    help="Bases with support equal to or greater then SSBT, will be classified as a strong supported base [Default: %(default)s]"
    )

snvArgs.add_argument(
    "-eds", "--enforce_dual_strand",
    action='store_true',
    help="require at least one molecule to be read in each of the two directions"
    )
snvArgs.add_argument(
    "-mq", "--min_qual",
    default=30,
    type=int,
    help="Minimum base quality threshold, below which a WEAK base will be ignored'")

# Filter Args
filterArgs = parser.add_argument_group("Filter Arguments")
filterArgs.add_argument("-sv", "--allow_single_stranded", action="store_true", default=False, help="Allow variants with only single stranded support [Default: %(default)s]")
filterArgs.add_argument("-sb", "--strand_bias_threshold", default=0.05, type=float, help="Strand bias p-value threshold, below which vairants will be discarded [Default: %(default)s]")
filterArgs.add_argument("-ss", "--strong_singleton_threshold", default=1, type=int, help="Base threshold for strong singleton bases (SP, SN) [Default: %(default)s]")
filterArgs.add_argument("-sd", "--strong_duplex_threshold", default=1, type=int, help="Base threshold for strong duplex bases (DPN, DPn, DpN) [Default: %(default)s]")
filterArgs.add_argument("-wt", "--weak_base_threshold", default=2, type=int, help="Base threshold for weak supporting bases (Sn, Sp, Dpn) [Default: %(default)s]")
filterArgs.add_argument("-md", "--min_depth", type=int, default=2, help="Minimum depth threshold [Default: %(default)s]")
filterArgs.add_argument("-fl", "--filter_log", metavar="FILE", help="A log file to explain the thresholds used for each variant, and why variants failed filters")
filterArgs.add_argument("-nv", "--normal_vaf", default=0.05, type=float, metavar="FLOAT", help="VAF threshold for the normal sample, above which variants will be called as germline [Default: %(default)s]")
filterArgs.add_argument("-g", "--germline_output", metavar="FILE", help="If a matched normal was supplied, an output file for germline variants")


def isValidFile(file, parser):
	"""
	Checks to ensure the specified file exists, and throws an error if it does not

	Args:
		file: A filepath
		parser: An argparse.ArgumentParser() object. Used to throw an exception if the file does not exist

	Returns:
		type: The file itself

	Raises:
		parser.error(): An ArgumentParser.error() object, thrown if the file does not exist
	"""
	if os.path.isfile(file):
		return file
	else:
		parser.error("Unable to locate %s" % (file))


def createLogFile(args, logFile="ProDuSe_Task.log", *versionInfo):
	"""
	Creates a log file in the output directory specifying the arguments that ProDuSe was run with

	If a log file already exists, it is simply appended with the new log

	Input:
		args: A namespace object listing ProDuSe arguments
		logFile: Name of the output file
		versionInfo: A list of lists containing version information
	"""

	log = []

	# Timestamp this
	log.append("ProDuSe initialized at %s %s" % (time.strftime('%x'), time.strftime('%X')))

	log.append("python " + sys.argv[0])

	# Add each argument to the logString
	for argument, parameter in vars(args).items():
		log.append(" --" + argument + " " + str(parameter))

	logString = "\n".join(log) + "\n"

	# Add the version information for each subprocess
	for program in versionInfo:
		logString += "\n" + "".join(program)

	with open(logFile, "a") as o:
		o.write(logString)
		o.write("\n")


def check_command(command, versionStr=None):
	"""
	Ensures the command is installed on the system, and returns the version if it does exist

	If a command is not found, python will throw an OS error. Catch this, and inform the user.
	Otherwise, save and return the version number

	Args:
		command: Literal name of the command
	        versionStr: The argument paseed to the command to print out the version
	Returns:
		version: An array listing the output of running '--version' or 'version'
	"""

	try:
		runCom = [command]
		if versionStr:
			runCom.append(versionStr)
		runCheck = subprocess.Popen(runCom, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		version = []
		for line in runCheck.stdout:
			line = line.decode("utf-8")
			version.append(line)

		return version

	# Prints out an error if the command is not installed
	except OSError:
		sys.stderr.write("ERROR: Unable to run %s\n" % (command))
		sys.stderr.write("Please ensure it is installed, and try again\n")
		sys.exit(1)


def getConfig(sampleDir, task):
	"""
		Returns the path to the config file for the designated task

		Args:
			sampleDir: Directory for produse analysis
			task: Name of the task
	"""
	configFile = os.path.join(sampleDir, "config", task + "_task.ini")
	if not os.path.exists(configFile):
		sys.stderr.write("ERROR: Unable to locate %s\n" % (configFile))
		sys.stderr.write("Check to ensure configure_produse completed sucessfully\n")
		sys.exit(1)
	return configFile


def runStitcher(inputBAM, stitcherPath):
	"""
	Runs Illumina's stitcher on the supplied BAM file

	Args:
		inputBAM: Path to BAM file to be stitched
		stitcherPath: Path to stitcher.exe
	Returns:
		outFile: Path to stitched BAM
	"""

	# Output the stitched BAM file in the same directory
	outDir = os.path.dirname(inputBAM)

	stitcherArgs = ["mono", stitcherPath, "--Bam", inputBAM, "--OutFolder=" + outDir]
	subprocess.check_call(stitcherArgs)

	return inputBAM.replace(".bam", ".stitched.bam")


def runSort(bamFile, byName=False):
	"""
	Sorts the supplied BAM file

	Sorts the designated BAM file using Samtools sort. The sorted BAM file is placed in the same
	directory as the input BAM file

	Args:
		bamFile: Path to BAM file to be sorted
		byName: Whether or not to sort the BAM file by read name
	Returns:
		outFile: Path to sorted BAM file
	"""
	outFile = bamFile.replace(".bam", ".sorted.bam")
	samtoolsArgs = ["samtools", "sort", "-o", outFile, bamFile]
	if byName:
		samtoolsArgs.append("-n")
	subprocess.check_call(samtoolsArgs)
	return outFile


def runPipeline(args, sampleName, sampleDir):
	"""
		Runs the main ProDuSe analysis stages on the provided sample

		Args:
			args: A namespace object listing command line parameters to be passed to subscripts
			sampleName: Name of the sample currently being processed
			sampleDir: Output directory
	"""

	printPrefix = "PRODUSE-MAIN\t"

	# Run Trim
	args.config = getConfig(sampleDir, "trim")
	trim.main(args)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Trimming Complete\n"]))

	# If we need to trim the adapters off of a matched normal, then run it
	if os.path.exists(sampleDir + os.sep + "config" + os.sep + "trim_normal_task.ini"):
		sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Trimming Normal\n"]))
		args.config = getConfig(sampleDir, "trim_normal")
		trim.main(args)
		sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Normal Trimming Complete\n"]))

	# Run bwa on the trimmed fastqs
	args.config = getConfig(sampleDir, "trim_bwa")
	bwa.main(args)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Alignment Complete\n"]))

        # If fastq files were supplied for the matched normal, run BWA and align it
	if os.path.exists(sampleDir + os.sep + "config" + os.sep + "bwa_normal_task.ini"):
		sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Aligning Normal\n"]))
		args.config = getConfig(sampleDir, "bwa_normal")
		bwa.main(args)
		sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Normal Alignment Complete\n"]))

	# Run collapse on the trimmed BAM file
	args.config = getConfig(sampleDir, "collapse")
	collapse.main(args)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Collapse Complete\n"]))

	# Run bwa on the collapsed
	args.config = getConfig(sampleDir, "collapse_bwa")
	bwa.main(args)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Alignment Complete\n"]))

	# Run stitcher
	collapsedBamFile = os.path.abspath(os.path.join(sampleDir, "tmp", sampleName + ".collapse.bam"))
	stitchedBam = runStitcher(collapsedBamFile, args.stitcherpath)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Stitching Complete\n"]))

	# Sort files prior to splitmerge
	runSort(stitchedBam, byName=True)
	runSort(collapsedBamFile, byName=True)

	args.config = getConfig(sampleDir, "splitmerge")
	splitMergeBam = os.path.join(sampleDir, "results", sampleName + ".SplitMerge.bam")
	SplitMerge.main(args)
	runSort(splitMergeBam)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": SplitMerge Complete\n"]))

	# Time for SNV calling, what everyone has been waiting for
	args.config = getConfig(sampleDir, "snv")
	snv.main(args)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": SNV Calling Complete\n"]))

	# Filter variants
	args.config = getConfig(sampleDir, "filter_produse")
	filter_produse.main(args)
	# runFilter(args.vaf, vcfFile, scriptDir + os.sep + "filter_produse.pl")
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Variant Filtering Complete\n"]))

	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": ProDuSe analysis Complete\n"]))


def main(args=None):

	"""
	Runs ALL steps of the ProDuSe pipeline on the supplied samples

	Args:
		args: A namespace object containing command line parameters

	"""

	if not args:
		args = parser.parse_args()

	# Checks command line arguments

	# First things first, lets make sure that the programs required to run ProDuSe are installed,
	# and pull out the version number of each
	bwaVer = check_command("bwa")
	samtoolsVer = check_command("samtools", "--version")
	monoVer = check_command("mono", "--version")
	pythonVer = check_command("python", "--version")

	printPrefix = "PRODUSE-MAIN\t"
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Starting ProDuSe\n"]))

	# Setup ProDuSe using configure_produse
	configure_produse.main(args)
	sampleList = configure_produse.samples
	outDir = configure_produse.output_directory
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Configuration Complete\n"]))

	createLogFile(args, os.path.join(outDir, "ProDuSe_Task.log"),["ProDuSe Version " + ProDuSeVer.__version__ + "\n"], samtoolsVer, monoVer, pythonVer, bwaVer)

	# Run the pipeline on each sample
	for sample in sampleList:

		sys.stderr.write("\t".join(["\n" + printPrefix, time.strftime('%X'), "Processing sample " + sample + "\n"]))
		sampleDir = os.path.join(outDir, sample)
		runPipeline(args, sample, sampleDir)
	sys.stderr.write("\t".join(["\n" + printPrefix, time.strftime('%X'), "All Samples Processed" + "\n"]))


if __name__ == "__main__":
	print("")
	main()
