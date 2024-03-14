import pysam
import json
import sys

sourceFQ = sys.argv[1]
fqRoot = sys.argv[2]
outName = sys.argv[3]

bam = pysam.AlignmentFile(sourceFQ, mode='rb')
bam_header = dict((k, v) for k, v in bam.header.items())
if type(bam_header['RG'])==list:
    dictRG = bam_header['RG'][0]
else:
    dictRG = bam_header['RG']

tagDict = {}
tags = ['SM','ID','LB','PU','DT','PL','CN']

for tag in tags:
    if tag in dictRG:
        tagDict[tag] = dictRG[tag]
    else:
        tagDict[tag] = 'NA'

jsonDict = {'ConvertPairedFastQsToUnmappedBamWf.sample_name':tagDict['SM'],
            'ConvertPairedFastQsToUnmappedBamWf.readgroup_name':tagDict['ID'],
            'ConvertPairedFastQsToUnmappedBamWf.fastq_1':fqRoot+'.1.fq.gz',
            'ConvertPairedFastQsToUnmappedBamWf.fastq_2':fqRoot+'.2.fq.gz',
            'ConvertPairedFastQsToUnmappedBamWf.library_name':tagDict['LB'],
            'ConvertPairedFastQsToUnmappedBamWf.platform_unit':tagDict['PU'],
            'ConvertPairedFastQsToUnmappedBamWf.run_date':tagDict['DT'],
            'ConvertPairedFastQsToUnmappedBamWf.platform_name':tagDict['PL'],
            'ConvertPairedFastQsToUnmappedBamWf.sequencing_center':tagDict['CN'],
            'ConvertPairedFastQsToUnmappedBamWf.make_fofn':'true'}

with open(outName, 'w') as fout:
    json.dump(jsonDict, fout)
