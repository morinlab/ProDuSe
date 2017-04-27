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

# Import stanard python modules
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

# If installed and running in python3
except ImportError:
	from ProDuSe import configure_produse
	from ProDuSe import bwa
	from ProDuSe import trim
	from ProDuSe import collapse
	from ProDuSe import snv
	from ProDuSe import filter_produse
	from ProDuSe import SplitMerge

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

# Used to maintain backwards compatibility with the poorly-named strand mis-match
adapterMismatch = collapseArgs.add_mutually_exclusive_group(required=True)
adapterMismatch.add(
    "-amm", "--adapter_max_mismatch",
    type=int,
    help="The maximum number of mismatches allowed between the expected and actual adapter sequences [Default: %(default)s]",
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

collapseArgs.add(
    "-sf", "--stats_file",
    type=str,
    required=False,
    default=None,
    help="An optional output file to list stats generated during collapse"
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
filterArgs.add_argument("-ss", "--allow_single_stranded", action="store_true", default=False, help="Allow variants with only single stranded support [Default: %(default)s]")
filterArgs.add_argument("-sb", "--strand_bias_threshold", default=0.05, type=float, help="Strand bias p-value threshold, below which vairants will be discarded [Default: %(default)s]")
filterArgs.add_argument("-st", "--strong_base_threshold", default=1, type=int, help="Strong supported base count threshold [Default: %(default)s]")
filterArgs.add_argument("-wt", "--weak_base_threshold", default=2, type=int, help="Weak supported base count theshold [Default: %(default)s]")
"""
filterArgs.add_argument("-dsv", "--totalvaf", type=float, default=0.05, help="Dual-strand VAF threshold [Default: %(default)s]")
filterArgs.add_argument("-md", "--min_duplex", type=int, default=3, help="Minimum duplex support required to call a variant [Default: %(default)s]")
filterArgs.add_argument("-mp", "--min_pos_strand", type=int, default=1, help="Minimum positive strand support required to call a varaint [Default: %(default)s]")
filterArgs.add_argument("-mn", "--min_neg_strand", type=int, default=1, help="Minimum negative strand support required to call a variant [Default: %(default)s]")
filterArgs.add_argument("-ms", "--min_singleton", type=int, default=3, help="Minimum singleton support required to call a variant [Default: %(default)s]")
filterArgs.add_argument("-ww", "--weak_base_weight", type=lambda x: isValidWeight(x, parser), default=0.1, help="Weight of weak bases, relative to strong supported bases [Default: %(default)s]")
filterArgs.add_argument("-sb", "--strand_bias", type=float, default=0.05, help="Strand bias p-value threshold, below which variants will be ignored")
"""


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
	Creates a log file in the cwd specifying the arguments that ProDuSe was run with

	Input:
		args: A namespace object listing ProDuSe arguments
		logFile: Name of the output file
		versionInfo: A list of lists containing version information
	"""

	logString = ["python ", sys.argv[0]]

	# Add each argument to the logString
	for argument, parameter in vars(args).items():
		logString.extend([" --" + argument + " ", str(parameter)])

	# Add the version information for each subprocess
	for program in versionInfo:
		logString.append("\n")
		logString.extend(program)

	with open(logFile, "w") as o:
		for item in logString:
			o.write(item)
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

def getSampleList(outDir):
	"""
	List all samples to be processed by ProDuSe

	Parse the names of the sample subdirectories created by configure_produse


	Args:
		outdir: ProDuSe output directory

	Returns:
		sampleList: A list of strings, containing sample names
	"""

	# Sanity Check: Ensure produse_analysis_directory actually exists
	if not os.path.exists(os.path.join(outDir, "produse_analysis_directory")):
		sys.stderr.write("ERROR: Unable to find " + os.path.join([outDir, "produse_analysis_directory"]) + "\n")
		sys.stderr.write("Check to ensure configure_produse completed sucessfully")
		exit(1)

	# Generates a list of samples by listing all directories created in 'produce_analysis_directory'
	sampleList = next(os.walk(os.path.join(outDir, "produse_analysis_directory")))[1]

	# If configure_produse.py created a reference folder, ignore that directory
	try:
		sampleList.remove("Reference_Genome")
	except ValueError:
		pass

	return sampleList


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
		sys.stderr.write("Check to ensure configure_produse completed sucessfully")
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


def runSplitMerge(inputBam, outputBam, scriptPath):
	"""
	Runs splitmerge.pl on the supplied BAM file

	Args:
		inputBam: Path to BAM file to be split
		outputBam: Path to output BAM file
		sciptPath: Path to splitmerge.pl
	"""

	# Sanity check. Ensure the script exists in that directory
	if not os.path.exists(scriptPath):
		sys.stderr.write("ERROR: Unable to locate \'splitmerge.pl\' in %s\n" % scriptDir)
		sys.stderr.write("Please ensure the script exists\n")

	preSMViewArgs = ["samtools", "view", "-h", inputBam]
	postSMViewArgs = ["samtools", "view", "-b"]

	# Runs splitmerge
	with open(outputBam, "w") as o:
		preSMView = subprocess.Popen(preSMViewArgs, stdout=subprocess.PIPE)
		exSM = subprocess.Popen(scriptPath, stdin=preSMView.stdout, stdout=subprocess.PIPE)

		subprocess.check_call(postSMViewArgs, stdin=exSM.stdout, stdout=o)


def runFilter(vaf, inputFile, scriptPath):
	"""
	Runs filter_produse.pl

	Filters variants called by snv.py using the specified VAF threshold

	Args:
		vaf: Minimum VAF threshold. Variants below this threshold will be discarded
		vcfFile: Path to input VCF file
		scriptPath: Path to filter_produse.pl
	"""
	filterArgs = ["perl", scriptPath, str(vaf)]
	outFile = inputFile.replace(".vcf", ".filtered.vcf")

	with open(inputFile) as f, open(outFile, "w") as o:
		subprocess.Popen(filterArgs, stdin=f, stdout=o)


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

	# Run bwa on the trimmed fastqs
	args.config = getConfig(sampleDir, "trim_bwa")
	bwa.main(args)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Alignment Complete\n"]))

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
	args.config = getConfig(sampleDir, "filter")
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
		# Clear these, since they are different for each script
		args.input = None
		args.output = None

	# Checks command line arguments

	# First things first, lets make sure that the programs required to run ProDuSe are installed,
	# and pull out the version number of each
	check_command("bwa")
	samtoolsVer = check_command("samtools", "--version")
	monoVer = check_command("mono", "--version")
	pythonVer = check_command("python", "--version")

	printPrefix = "PRODUSE-MAIN\t"
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Starting ProDuSe\n"]))

	# Setup ProDuSe using configure_produse
	configure_produse.main(args)
	sampleList = getSampleList(args.output_directory)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Configuration Complete\n"]))

	createLogFile(args, os.path.join(os.path.abspath(args.output_directory), "produse_analysis_directory", "ProDuSe_Task.log"), samtoolsVer, monoVer, pythonVer)

	# Run the pipeline on each sample
	for sample in sampleList:

		sys.stderr.write("\t".join(["\n" + printPrefix, time.strftime('%X'), "Processing sample " + sample + "\n"]))
		sampleDir = os.path.join(os.path.abspath(args.output_directory), "produse_analysis_directory", sample)
		runPipeline(args, sample, sampleDir)
	sys.stderr.write("\t".join(["\n" + printPrefix, time.strftime('%X'), "All Samples Processed" + "\n"]))


if __name__ == "__main__":
	print("")
	main()
