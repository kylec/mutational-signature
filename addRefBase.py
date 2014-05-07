#!/bin/python
#Author: Kyle Chang
#Count number of mutations per mutation spectrum

import string
import pysam
import argparse
from collections import OrderedDict

# read snp
# find ref  base around snp
# tally each base change type

# parse opt
parser = argparse.ArgumentParser()
parser.add_argument("input")
args = parser.parse_args()
#print "input=", args.input

# global var
ref = "/home/kchang3/references/hsap_36.1_hg18.fa"
subtypes = OrderedDict() # count seq context and var combinations

# reverse complement
def revcomp(seq):
  rc_seq = seq.translate(string.maketrans("ATGCN","TACGN"))[::-1]
  return rc_seq

# load 96 subtypes
fin = open("/home/kchang3/src/mutation_signature/input/96_mut_types.txt")
for line in fin:
  subtypes[line.strip()] = 0

# read snp and write output
fin = open(args.input, 'r')
fsubtypes = open(args.input + ".subtypes", 'w')

for line in fin:
  if line.startswith("#") | line.startswith("contig"):
    continue
  val = line.strip().split("\t")
  # skip rejected calls
  if val[34] == "REJECT":
    continue
  # replace x in CCCxAAC with ref(T) then substring to just CTA, set as key in dict and count occurences
  seqContext = val[2].replace("x",val[3])[2:5]
  key = seqContext + val[4]
  if key in subtypes:
    subtypes[key] += 1
  else:
    #print "lalala" + revcomp(seqContext) + "," + revcomp(val[4])
    rc_key = revcomp(seqContext) + revcomp(val[4])
    #a=rc(seqContext)
   
    subtypes[rc_key] += 1
  
# print seq context count
for seq, num in subtypes.iteritems(): 
  fsubtypes.write(str(num) + "\n")
  #fsubtypes.write(seq + "\t" + str(num) + "\n")
  
  #print val[2], seqContext, val[34]
#fa = pysam.Fastafile(ref)
#print fa.fetch('chr1', 57, 64)
