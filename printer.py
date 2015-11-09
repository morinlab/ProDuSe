import time
import sys

program = 'PRODUSE'
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
        global_task = 'START   '
        general_print(['BARCODE ADAPTER PIPELINE'])
        general_print(['VERSION', version])
        general_print(['AUTHOUR', authour])
    elif task == 'ARGS':
        global_task = 'ARGS    '
        general_print(['FETCHING COMMAND ARGUMENTS'])
    elif task == 'SETUP':
        global_task = 'SETUP   '
        general_print([' '])
    elif task == 'TRIM':
        global_task = 'TRIM    '
        general_print(['TRIMMING ADAPTER SEQUENCE'])
    elif task == 'BWA':
        global_task = 'BWA     '
        general_print(['ALIGNING TO REFERENCE'])
    elif task == 'COLLAPSE':
        global_task = task
        general_print(['COLLAPSING ADAPTER DUPLICATES'])
    else:
        general_print([task])

def command(commands):
    general_print([' '.join(commands)])

def trim(count, discard):
    print '    '.join([program, time.strftime('%X'), 'TRIM  '] + [':'.join(['Reads Processed', str(count)])] + [':'.join(['Reads Discarded', str(discard)])])


class open:

    def __init__(self, pipeline, pipeline_stage, stage_task, list_of_counters, verbose):
        self.pipeline = pipeline
        self.stage = pipeline_stage
        self.task = stage_task
        self.counters = {}
        
    def pipeline(self, pipeline):
        self.pipeline = pipeline

    def stage(self, stage):
        self.stage = stage

    def task(self, task):
        self.task = tas
        if self.verbose == 3:
            print '    '.join([self.pipeline, time.strftime('%X'), self.stage, self.task,] + [':'.join([x, self.counter[x]]) for x in self.counter])

    def count(self, counter):
        self.counter[counter] =+ 1
        self.base =+ 1

    def update(self):       
        
        if self.verbose == 0:
            pass

        elif self.verbose == 1:
            if self.base % 100000:
                print '    '.join([self.pipeline, time.strftime('%X'), self.stage, self.task,] + [':'.join([x, self.counter[x]]) for x in self.counter])     


