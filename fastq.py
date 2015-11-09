import itertools
import printer

class Fastq:
    
    def __init__(self, id, seq, strand, qual):
        self.id = id
        self.seq = seq
        self.strand = strand
        self.qual = qual

    def __str__(self):
        return '\n'.join([
            ':\t'.join(['ID', self.id]),
            ':\t'.join(['SEQ', self.seq]),
            ':\t'.join(['STRAND', self.strand]),
            ':\t'.join(['QUAL', self.qual])
            ])

    def trim(self, length):
        
        tmp = Fastq('','','','')
        tmp.id = self.id
        tmp.strand = self.strand

        if length > 0:
            
            # 
            tmp.seq = (self.seq[:length])
            tmp.qual = (self.qual[:length])
            
            #
            self.seq = self.seq[length:]
            self.qual = self.qual[length:]

        elif length < 0:
            
            #
            tmp.seq = self.seq[length:]
            tmp.append = self.qual[length:]

            #
            self.seq = self.seq[:length]
            self.qual = self.qual[:length]      

        else:
            pass
        return tmp

    def glue(self, fastq):

        # Probably can do more high level stuff using the strand
        self.seq = ''.join(self.seq, fastq.seq)
        self.qual = ''.join(self.qual, fastq.qual)


class FastqRead:

    def __init__(self, file):
        self.fh = open(file, 'r')
        self.data = list(itertools.islice(self.fh, 20000))
        self.end = len(self.data) < 20000
        self.count = len(self.data)

    def __iter__(self):
        return self

    def next(self):

        if self.end and len(self.data) == 0:
            raise StopIteration()

        elif len(self.data) == 0:
            self.data = list(itertools.islice(self.fh, 20000))
            self.end = len(self.data) < 20000
            self.count = self.count + len(self.data)
            return self.next()

        else:
            return Fastq(
                self.data.pop(0)[:-1], 
                self.data.pop(0)[:-1], 
                self.data.pop(0)[:-1], 
                self.data.pop(0)[:-1]
                )

    def __del__(self):
        self.close()

    def __exit__(self):
        self.close()

    def close(self):
        self.fh.close()


class FastqWrite:

    def __init__(self,file):
        self.file = file
        self.data = []
        self.data_len = len(self.data)
        self.count = len(self.data)

    def force_write(self):
        tmp = open(self.file, 'a')
        tmp.write(''.join(self.data))
        self.data = []
        tmp.close()

    def next(self, fastq):

        self.data.append(
            ''.join([fastq.id,'\n'])
            )
        self.data.append(
            ''.join([fastq.seq,'\n'])
            )
        self.data.append(
            ''.join([fastq.strand,'\n'])
            )
        self.data.append(
            ''.join([fastq.qual,'\n'])
            )
        if len(self.data) >= 20000:
            self.count = self.count + len(self.data)
            self.force_write()

    def __exit__(self):
        self.force_write()

    def __del__(self):
        self.force_write()

    def close(self):
        self.force_write()


def FastqOpen(file, mode):
    if mode == 'r':
        return FastqRead(file)

    elif mode == 'w':
        return FastqWrite(file)

    else:
        sys.exit(':'.join(['This mode is currently not supported by FastqOpen', mode]))
