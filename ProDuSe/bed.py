import itertools

class BedRead(file):

    def __init__(self, file):
        self.fh = open(file, 'r');
        self.data = list(itertools.islice(self.fh, 20000));
        self.data = [ x[:-1] for x in self.data ];
        self.chroms = [ str(x.split('\t')[0]) for x in self.data ]
        self.starts = [ int(x.split('\t')[1]) for x in self.data ]
        self.ends = [ int(x.split('\t')[2]) for x in self.data ]
        self.regions = [ ''.join( [ str(x.split('\t')[0]), ':', str(x.split('\t')[1]), '-', str(x.split('\t')[2])  ] ) for x in self.data ];

    def has(self, chrom, start, end):
        for i in range(len(self.chroms)):
            if not self.chroms[i] == chrom:
                continue
            else:
                if self.starts[i] >= start and self.ends[i] <= end:
                    return True
        return False 
            

def BedOpen(file, mode):
    if mode == 'r':
        return BedRead(file)
    else:
        sys.exit(':'.join(['This mode is currently not supported by BedOpen', mode]))




