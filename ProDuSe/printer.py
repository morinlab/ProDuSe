import time
import sys
import inspect
import os

function_to_prefix = {}
function_to_prefix["collapse.py"] = "PRODUSE-COLLAPSE   "
function_to_prefix["trim.py"] =     "PRODUSE-TRIM       "
function_to_prefix["bwa.py"] =      "PRODUSE-BWA        "
function_to_prefix["snv.py"] =      "PRODUSE-SNV        "

def message(message, error=False):
    function = os.path.basename(inspect.stack()[2][1])
    if not error:
        joined_message = ''.join([
            function_to_prefix[function],
            time.strftime("%X"),
            "    ",
            message,
            "\n"
            ])
        sys.stdout.write(joined_message)
    else:
        joined_message = ''.join([
            function_to_prefix[function],
            time.strftime("%X"),
            "    Error: ",
            message,
            "\n"
            ])
        sys.stdout.write(joined_message)
        sys.exit(1)      