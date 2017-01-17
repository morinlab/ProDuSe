import requests
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest='input')
    parser.add_argument('-o', dest='output')
    args = parser.parse_args()
    
    url = 'http://www.broadinstitute.org/oncotator/mutation/'
    
    # open files
    try:
        mafFile = open(args.input, 'r')
    except IOError:
        print 'Error: maf input does not exist.'    
    annoFile = open(args.output, 'w')

    # headers
    annoHeader = 'gene\ttranscript\tvariant_class\tdbsnp\taa_change\tpph2_class'
    mafHeader = mafFile.readline()
    annoFile.write(mafHeader)
    
    # read maf file for chr, start, end, ref, var
    i = 0
    for row in mafFile:
        i += 1
        fields = row.rstrip('\n').split('\t')
        query = url + fields[4] + '_' + fields[5] + '_' + fields[6] + '_' + fields[10] + '_' + fields[12]
        
        query = requests.get(query)
        result = query.json()
        AF = result["1000Genome_AF"]
        print(len(fields))
        try:
            if float(AF) > 0.01:
                fields[108] = "common_variant"
        except ValueError:
            print "what's up with %s" % AF
        fields[107] = AF
        out_s = '\t'.join([str(x) for x in fields])
        annoFile.write(out_s)
        annoFile.write("\n")
        
    
    mafFile.close()
    annoFile.close()

if __name__ == '__main__':
    main()

