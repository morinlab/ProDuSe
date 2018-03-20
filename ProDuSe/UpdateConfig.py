#!/usr/bin/env python

# DESCRIPTION:
#   Converts the old-style ProDuSe config file to an updated version
#
# AUTHOR
#   Christopher Rushton (ckrushto@sfu.ca)


import os
import argparse
import sys
from collections import OrderedDict

def isValidFile(file, parser):
    """
    Checks to ensure the specified file path is valid

    :param file: A string containing a filepath
    :param parser: An argparse.ArgumentParser() object
    :return: file, if the specified file exists
    :raises parser.error: If the specified file does not exist
    """
    if os.path.exists(file):
        return file
    parser.error("Unable to locate \'%s\': No such file or directory." % file)


parser = argparse.ArgumentParser(description="Converts older config files to a format compatible with this version of ProDuSe")
parser.add_argument("-i", "--input", metavar="INI", type=lambda x:isValidFile(x, parser), required=True, help="Input old-style configuration file")
parser.add_argument("-o", "--output", metavar="INI", required=True, help="Output reformatted configuration file")

def main(args=None, sysStdin=None):
    if args is None:
        args = parser.parse_args(sysStdin)

    outLines = OrderedDict()
    unformattedLines = OrderedDict()
    i = 0  # Keep the arguments in the same order in the output file
    # Reformat the INI file
    with open(args.input) as f, open(args.output, "w") as o:
        # Since the new config file only requires one section header, place it here
        o.write("[Pipeline]" + os.linesep)
        for line in f:
            # Ignore the old config sections. They are no longer required
            if line[0] == "[":
                configSection = line.rstrip()
                continue
            # Ignore empty lines.
            if line == os.linesep:
                continue
            # Write out comment lines unaltered
            if line[0] == "#":
                unformattedLines[i] = line
                i += 1
                continue
            # Collapse args
            oLine = line.replace("adapter_max_mismatch", "family_mismatch")
            oLine = oLine.replace("duplex_max_mismatch", "duplex_mismatch")
            oLine = oLine.replace("strand_position", "family_mask")
            oLine = oLine.replace("duplex_position", "duplex_mask")
            # Rename "adapter" to "barcode"
            oLine = oLine.replace("adapter_", "barcode_")
            try:
                argument, parameter = oLine.split("=")
            except ValueError:  # This doesn't look like a parameter. Ignore it, but warn the user
                sys.stderr.write("WARNING: Unrecognized entry \'%s\'. Ignoring..." % line.strip(os.linesep) + os.linesep)
                continue
            if argument in outLines:  # This argument has already been encountered
                if parameter == outLines[argument]:  # These entries are identical. One can safely be ignored
                    continue
                else:
                    # The argument name is the same, but the parameters are different. Since I don't know which one
                    # to keep, keep them both, but inform the user that they MUST remove one of the duplicates
                    sys.stderr.write("WARNING: Multiple parameters for \'%s\' were provided" % argument + os.linesep)
                    sys.stderr.write("Only one parameter should be specified for each argument."
                                     " Please remove one of the entires from the output config file" + os.linesep)
                    unformattedLines[i] = line
            else:
                outLines[argument] = parameter
            i += 1

        # Write the output config file
        i = 0
        for argument, parameter in outLines.items():
            while i in unformattedLines:
                o.write(unformattedLines[i])
                del unformattedLines[i]
                i += 1
            o.write("=".join([argument, parameter]))
            i += 1
        for line in unformattedLines.values():
            o.write(line)

if __name__ == "__main__":
    main()