import itertools

class BedRead(file):

    def __init__(self, file):
        self.fh = open(file, 'r');
        self.data = list(itertools.islice(self.fh, 20000));
        self.data = [ x[:-1] for x in self.data ];
        self.data = [ ':'.join( [ str(x.split('\t')[0]) , str(x.split('\t')[1]) ] ) for x in self.data ];
        print self.data

    def __contains__(self, other):
    	if type(other) == str:
    		return other in self.data;
    	else:
    		return False;


def BedOpen(file, mode):
    if mode == 'r':
        return BedRead(file)
    else:
        sys.exit(':'.join(['This mode is currently not supported by BedOpen', mode]))




