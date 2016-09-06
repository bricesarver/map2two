from collections import defaultdict
import pysam
import os
import sys

#one = "1_motile.merged.CZII.dedup.realigned.bam"
#two = "1_motile.merged.PWK.dedup.realigned.bam"

def main(one, two):
    first = defaultdict(list)
    errors = defaultdict(list)
    
    print("processing first file...")
    
    with pysam.AlignmentFile(one, "rb") as bamfile:
        for i, read in enumerate(bamfile):
            #print(read)
            if i % 100000 == 0:
            	print("first file, read: " + str(i))
            readname = read.qname
            try:
            	first[readname].extend((read.mapq, read.qlen, read.reference_name, read.pos, read.aend))
            except:
                errors[readname].extend(readname)
                pass
    
    first_copy = dict(first)
    
    print("processing second file...")
    
    with pysam.AlignmentFile(two, "rb") as bamfile:
        for i, read in enumerate(bamfile):
            #print(read)
            if i % 100000 == 0:
            	print("second file, read: " + str(i))
            readname = read.qname
            if readname in first_copy:
            	try:
                    first[readname].extend((read.mapq, read.qlen, read.reference_name, read.pos, read.aend))
                except:
                    errors[readname].extend(readname)
                    pass
    
    print("calculating PE keys in both...")
    twenty = [key for key in first if len(first[key]) == 20]
    
    print("calculating SE keys in both...")
    ten = [key for key in first if len(first[key]) == 10]
    
    print("writing PE keys in both...")
    with open(one.split(".")[0] + ".results.20.txt", "w") as results:
        results.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("readID", "mapq.1", "qlen.1", "chr.1", "start.1", "end.1", "mapq.2", "qlen.2", "chr.2", "start.2", "end.2", "mapq.3", "qlen.3", "chr.3", "start.3", "end.3", "mapq.4", "qlen.4", "chr.4", "start.4", "end.4"))
        for key in twenty:
            results.write(key + "\t" + "\t".join([str(x) for x in first[key]]) + "\n")
    
    print("writing SE keys in both...")
    with open(one.split(".")[0] + ".results.10.txt", "w") as results:
        results.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("readID", "mapq.1", "qlen.1", "chr.1", "start.1", "end.1", "mapq.2", "qlen.2", "chr.2", "start.2", "end.2"))
        for key in ten:
            results.write(key + "\t" + "\t".join([str(x) for x in first[key]]) + "\n")

    print("writing errors...")
    with open(one.split(".")[0] + ".errors.txt", "w") as errorfile:
        for key in errors:
        	errorfile.write(key + "\n")

if __name__ == "__main__":
    one = sys.argv[1]
    two = sys.argv[2]
    #if not os.path.isfile(one + ".bai"):
        #pysam.index(one)
    #if not os.path.isfile(two + ".bai"):
        #pysam.index(two)
    main(one, two)

