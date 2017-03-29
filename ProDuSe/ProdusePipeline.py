#! /usr/bin/env python

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

# If installed and running in python3
except ImportError:
	from ProDuSe import configure_produse
	from ProDuSe import bwa
	from ProDuSe import trim
	from ProDuSe import collapse
	from ProDuSe import snv
	from ProDuSe import filter_produse

"""
Processes command line arguments

"""

# Used to locate perl scripts required for running ProDuSe
global scriptDir
scriptDir = os.path.abspath(os.path.dirname(__file__))

# Pulls command like arguments from the sub-scripts
parser = configargparse.ArgumentParser(description="Runs the entire ProDuSe pipeline on the supplied samples. ", parents=[configure_produse.parser, trim.parser, bwa.parser, collapse.parser, snv.parser, filter_produse.parser], conflict_handler="resolve")
parser.add_argument("-c", "--config", required=False, is_config_file=True, help="ProDuSe config file, listing adapter sequences, positions, and other parameters to be passed to pipeline scripts. See (See \'etc/produse_config.ini\' for an example)")
parser.add_argument("-x", "--stitcherpath", required=True, type=lambda x: isValidFile(x), help="Path to Illumina's Stitcher.exe (Can be obtained from \'https://github.com/Illumina/Pisces\')")
parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads to use while running BWA [Default: %(default)s]")

# Supress command line arguments for these options, as there are irrelevent for this wrapper
parser.add_argument("-i", "--input", help=configargparse.SUPPRESS)
parser.add_argument("-o", "--output", help=configargparse.SUPPRESS)


def isValidFile(file):
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


def runSort(bamFile):
	"""
	Sorts the supplied BAM file

	Sorts the designated BAM file using Samtools sort. The sorted BAM file is placed in the same
	directory as the input BAM file

	Args:
		bamFile: Path to BAM file to be sorted
	Returns:
		outFile: Path to sorted BAM file
	"""
	outFile = bamFile.replace(".bam", ".sorted.bam")
	samtoolsArgs = ["samtools", "sort", "-o", outFile, bamFile]
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
	collapsedBamFile = os.path.abspath(os.path.join(sampleDir, "data", sampleName + ".collapse.bam"))
	stitchedBam = runStitcher(collapsedBamFile, args.stitcherpath)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": Stitching Complete\n"]))

	# Run splitmerge
	sortedStitchedBam = runSort(stitchedBam)
	global scriptDir
	splitMergeBam = os.path.join(sampleDir, "results", sampleName + ".SplitMerge.bam")
	runSplitMerge(sortedStitchedBam, splitMergeBam, scriptDir + os.sep + "splitmerge.pl")
	sortedSplitmergeBam = runSort(splitMergeBam)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": SplitMerge Complete\n"]))

	# Time for SNV calling, what everyone has been waiting for
	args.config = getConfig(sampleDir, "snv")
	snv.main(args)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), sampleName + ": SNV Calling Complete\n"]))

	# Filter variants
	vcfFile = os.path.join(sampleDir, "results", sampleName + ".variants.vcf")
	outFile = vcfFile.replace(".vcf", ".filtered.vcf")
	args.input=vcfFile
	args.output=outFile
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
	perlVer = check_command("perl", "--version")
	pythonVer = check_command("python", "--version")

	printPrefix = "PRODUSE-MAIN\t"
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Starting ProDuSe\n"]))

	# Setup ProDuSe using configure_produse
	configure_produse.main(args)
	sampleList = getSampleList(args.output_directory)
	sys.stderr.write("\t".join([printPrefix, time.strftime('%X'), "Configuration Complete\n"]))

	createLogFile(args, os.path.join(os.path.abspath(args.output_directory), "produse_analysis_directory", "ProDuSe_Task.log"), samtoolsVer, monoVer, perlVer, pythonVer)

	# Run the pipeline on each sample
	for sample in sampleList:

		sys.stderr.write("\t".join(["\n" + printPrefix, time.strftime('%X'), "Processing sample " + sample + "\n"]))
		sampleDir = os.path.join(os.path.abspath(args.output_directory), "produse_analysis_directory", sample)
		runPipeline(args, sample, sampleDir)
	sys.stderr.write("\t".join(["\n" + printPrefix, time.strftime('%X'), "All Samples Processed" + "\n"]))

if __name__ == "__main__":
	print("")
	main()
