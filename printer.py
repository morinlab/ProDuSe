import time
import sys

program = 'PROCESS-FASTQ'
version = "1.0"
authour = 'Marco Albuquerque'

def general_print( what ):
    print '    '.join([program, time.strftime('%X'), global_task] + what)

def issue( issue ):
    print '    '.join([program, time.strftime('%X'), 'ISSUE '] + [issue])
    sys.exit()

def general( task ):
    global global_task
    if task == 'START':
        global_task = 'START '
        general_print(['BARCODE ADAPTER PIPELINE'])
        general_print(['VERSION', version])
        general_print(['AUTHOUR', authour])
    elif task == 'ARGS':
        global_task = 'ARGS  '
        general_print(['FETCHING COMMAND ARGUMENTS'])
    elif task == 'SETUP':
        global_task = 'SETUP '
        general_print([' '])
    elif task == 'TRIM':
        global_task = 'TRIM  '
        general_print(['TRIMMING ADAPTER SEQUENCE'])
    elif task == 'BWA':
        global_task = 'BWA   '
        general_print(['ALIGNING TO REFERENCE'])
    elif task == 'PICARD':
        global_task = task
        general_print(['MARKING DUPLICATES'])
    elif task == 'REMARK':
        global_task = task
        general_print(['REMARKING ADAPTER DUPLICATES'])
    else:
        general_print([task])

def fastq(count, mode):
    if mode == 'r':
        print '    '.join([program, time.strftime('%X'), 'FASTQ '] + [': '.join(['Fastq Reads', str(count)])])
    elif mode == 'w':
        print '    '.join([program, time.strftime('%X'), 'FASTQ '] + [':'.join(['Fastq Writes', str(count)])])

def trim(count, discard):
    print '    '.join([program, time.strftime('%X'), 'TRIM  '] + [':'.join(['Reads Processed', str(count)])] + [':'.join(['Reads Discarded', str(discard)])])

def mark(count, discard):
    print '    '.join([program, time.strftime('%X'), 'MARK  '] + [':'.joing(['Adapter Classes Processed', str(count)])]) 
