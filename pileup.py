import pysam
import argparse

#parse opt
parser = argparse.ArgumentParser()
parser.add_argument("snv")
parser.add_argument("bam")
args = parser.parse_args()

ref = "/home/kchang3/references/hsap_36.1_hg18.fa"

# open bam
samfile = pysam.Samfile(args.bam, "rb")

# read snv file
print args.snv
fin = open(args.snv, "r")
for line in fin:
  if line.startswith("#") | line.startswith("contig"):
    continue
  arr = line.strip().split("\t")
  chr = arr[0]
  pos = int(arr[1])
  ref = arr[3]
  var = arr[4]
  ref_num = 0
  var_num = 0

  # parse pileup
  for column in samfile.pileup(chr, pos-1, pos): 
    for read in column.pileups:
      snpcaller = pysam.IteratorSNPSNPCalls(column)
      sc = snpcaller.call(chr, pos-1) 
      print str(sc)
   


