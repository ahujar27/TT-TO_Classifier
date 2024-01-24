#import pysam
import json
import sys

sourceBam = sys.argv[1]
outName = sys.argv[3]

jsonDict = {'BamToUnmappedBams.input_bam':sourceBam}

with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
