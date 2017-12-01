#!/usr/bin/python

'''
    Test Lara 1 versus Lara 2 on Bralibase 2 data sets.
    
    Author: Joerg Winkler <j.winkler@fu-berlin.de>
'''

import re
import os
import sys
import time
import numpy
import filecmp
import subprocess
from Bio import AlignIO
from Bio.Statistics.lowess import lowess
from shutil import copyfile
from matplotlib import use
use('PDF')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from datetime import datetime

######################
##  SET FILE PATHS  ##
######################

l2args = [[], ["-my", "0.5"], ["-my", "2"]]
'''
    ["-g", "2", "-lbm", "1", "-tb", "0.1", "-ssc", "10", "-stsc", "1", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-tb", "0.1", "-ssc", "5", "-stsc", "1", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-tb", "0.1", "-ssc", "1", "-stsc", "1", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-tb", "0.1", "-ssc", "1", "-stsc", "0", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-ssc", "10", "-stsc", "1", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-ssc", "5", "-stsc", "1", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-ssc", "1", "-stsc", "0", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-ssc", "1", "-stsc", "0", "-tcm", "0"],
          ["-lbm", "1", "-tb", "0.1", "-ssc", "10", "-stsc", "1", "-tcm", "0"],
          ["-lbm", "1", "-tb", "0.1", "-ssc", "5", "-stsc", "1", "-tcm", "0"],
          ["-lbm", "1", "-tb", "0.1", "-ssc", "1", "-stsc", "1", "-tcm", "0"],
          ["-lbm", "1", "-tb", "0.1", "-ssc", "1", "-stsc", "0", "-tcm", "0"],
          ["-lbm", "1", "-tb", "0.1", "-ssc", "10", "-stsc", "1", "-tcm", "2"],
          ["-lbm", "1", "-tb", "0.1", "-ssc", "5", "-stsc", "1", "-tcm", "2"],
          ["-lbm", "1", "-tb", "0.1", "-ssc", "1", "-stsc", "1", "-tcm", "2"],
          ["-lbm", "1", "-tb", "0.1", "-ssc", "1", "-stsc", "0", "-tcm", "2"],
          ["-g", "2", "-lbm", "1", "-tb", "0.1", "-ssc", "10", "-stsc", "1", "-lsm", "../ribosum_matrices/RIBOSUM45_N.txt", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-tb", "0.1", "-ssc", "1", "-stsc", "1", "-lsm", "../ribosum_matrices/RIBOSUM45_N.txt", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-tb", "0.1", "-ssc", "1", "-stsc", "0", "-lsm", "../ribosum_matrices/RIBOSUM45_N.txt", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-tb", "0.1", "-ssc", "10", "-stsc", "1", "-lsm", "../ribosum_matrices/RIBOSUM95_N.txt", "-tcm", "0"],
          ["-g", "2", "-lbm", "1", "-tb", "0.1", "-ssc", "1", "-stsc", "1", "-lsm", "../ribosum_matrices/RIBOSUM95_N.txt", "-tcm", "0"]]'''

# define macros
SCORES = (SPS, SCI, MPI) = (0, 1, 2)
SCORELBL = ('Sum of Pairs Score (compalignp)', 'Structure Conservation Index (RNAz)', 'Mean Pairwise Identity (RNAz)')
PROGLBL = ['Reference', 'LaRA1 + T-Coffee', 'SeqAn::T-Coffee', 'MAFFT', 'SeqAn::Lara + T-Coffee', 'SeqAn::LaRA + SeqAn::T-Coffee', 'SeqAn::LaRA + SeqAn::T-Coffee my=0.5', 'SeqAn::LaRA + SeqAn::T-Coffee my=0.5']
#PROGLBL.extend(['SeqAn::LaRA + SeqAn::T-Coffee' + str(i) for i in range(len(l2args) - 1)])
PROGRAMS = range(len(PROGLBL))
(REF, LA1, STC, MAF, L2O, L2N) = PROGRAMS[:6]

work_dir = os.getcwd()
results_dir  = os.path.join(work_dir, "results")
# old lara and old tcoffee
oldlara_dir = os.path.join(work_dir, "lara-1.3.2")
oldlara_bin = os.path.join(oldlara_dir, "lara")
oldtcof_bin = os.path.join(work_dir, "..", "..", "..", "..", "..", "apps", "tcoffee", "bin", "t_coffee")

# new lara and new tcoffee
seqan_dir = os.path.join(work_dir, "..", "..", "..", "..", "build", "bin") # adapt this to your system!
newlara_bin = os.path.join(seqan_dir, "laragu")
newtcof_bin = os.path.join(seqan_dir, "seqan_tcoffee")
tc_tempfile = os.path.join(results_dir, "tcoffeLara.lib")

# compiled MAFFT program
mafft_bin = os.path.join(work_dir, "..", "..", "..", "..", "..", "apps", "tcoffee", "plugins", "linux", "mafft")
# compiled RNAz program ( http://www.tbi.univie.ac.at/~wash/RNAz/ )
rnaz_bin = os.path.join(work_dir, "..", "..", "..", "..", "..", "apps", "rnaz", "bin", "RNAz") # adapt
# compiled compalignp program ( http://www.biophys.uni-duesseldorf.de/bralibase/compalignp.tgz )
compali_bin = os.path.join(work_dir, "..", "..", "..", "..", "..", "apps", "compalignp", "compalignp") # adapt

########################
##  CREATE FILE LIST  ##
########################

files = []
for d in "g2intron", "rRNA", "tRNA", "U5":
  p = os.path.join(work_dir, "benchmarks", "bralibase2", "data-set1", d, "unaligned")
  files.extend([(os.path.join(p, f), os.path.join(results_dir, d + "-" + f)) for f in os.listdir(p)])

if not os.path.isdir(results_dir):
  os.mkdir(results_dir)
  
print(str(len(files)) + " alignments to compute.")

############################
##  CALCULATE ALIGNMENTS  ##
############################

stats = {}
progtime = [0.0] * len(PROGRAMS)

for (infile, outfile) in files:
  basename = os.path.basename(outfile)
  print >>sys.stderr, "processing", basename

  # Reference
  refalignment = infile.replace("unaligned", "structural-no_str")
  copyfile(refalignment, outfile + str(REF) + ".fasta")

  # run Lara1
  if not os.path.isfile(outfile + str(LA1) + ".fasta"):
    print >>sys.stderr, "  --> Lara 1"
    t = time.time()
    os.chdir(oldlara_dir)
    exe = [oldlara_bin, "-i", infile, "-w", outfile + str(LA1) + ".fasta"]
    os.system(' '.join(exe) + " >/dev/null")
    os.chdir(work_dir)
    progtime[LA1] += time.time() - t

  # run SeqAn::TCoffee without Lara
  if not os.path.isfile(outfile + str(STC) + ".fasta"):
    print >>sys.stderr, "  --> SeqAn::TCoffee"
    t = time.time()
    exe = [newtcof_bin, "-s", infile, "-a", "iupac", "-b", "wavg", "-o", outfile + str(STC) + ".fasta"]
    os.system(' '.join(exe))
    progtime[STC] += time.time() - t

  # run MAFFT
  if not os.path.isfile(outfile + str(MAF) + ".fasta"):
    print >>sys.stderr, "  --> MAFFT"
    t = time.time()
    os.system(mafft_bin + " " + infile + " >" + outfile + str(MAF) + ".fasta 2>/dev/null")
    progtime[MAF] += time.time() - t

  # run Lara2
  sys.stderr.write('  --> Lara 2 ')
  for i in range(len(l2args)):
    if os.path.isfile(outfile + str(L2N + i) + ".fasta") and os.path.isfile(outfile + str(L2O) + ".fasta"):
      continue
    
    sys.stderr.write('.')
    sys.stderr.flush()
    libfile = outfile + str(L2N + i) + ".lib"
    exe = [newlara_bin, "-i", infile, "-td", results_dir, "-t", "1"]
    exe.extend(l2args[i])
    t = time.time()
    os.system(' '.join(exe))
    lt = time.time() - t
    
    # SeqAn::TCoffee
    exe = [newtcof_bin, "-s", infile, "-l", tc_tempfile, "-o", outfile + str(L2N + i) + ".fasta", "-a", "iupac", "-m", "global"]
    t = time.time()
    os.system(' '.join(exe))
    progtime[L2N + i] += time.time() - t + lt
    
    # old tcoffee
    if i == 0:
      exe = [oldtcof_bin, "-lib", tc_tempfile, "-output", "fasta", "-outfile", outfile + str(L2O) + ".fasta", "-newtree", outfile + str(L2O) + ".dnd"]
      t = time.time()
      os.system(' '.join(exe) + " >/dev/null 2>/dev/null")
      progtime[L2O] += time.time() - t + lt
  sys.stderr.write('\n')
  sys.stderr.flush()

  # transform into Clustal format and run RNAz analysis
  print >>sys.stderr, "  --> RNAz"
  for file in [outfile + str(x) for x in PROGRAMS]:
    if not os.path.isfile(file + ".aln") or not os.path.isfile(file + ".stat"):
      ali = AlignIO.read(file + ".fasta", "fasta")
      AlignIO.write([ali], file + ".aln", "clustal")
      os.system(rnaz_bin + " -o " + file + ".stat -n " + file + ".aln")

#  for file in [outfile + str(x) for x in PROGRAMS]:
#    ali = AlignIO.read(file + ".aln", "clustal")
#    AlignIO.write([ali], file + ".test", "fasta")
#    if not filecmp.cmp(file + ".fasta", file + ".test", False):
#      print >>sys.stderr, "      File " + file + " not converted properly."

##########################
##  ALIGNMENT ANALYSIS  ##
##########################

print ("\nAnalyze alignments")
for (infile, outfile) in files:
  sys.stderr.write('.')
  sys.stderr.flush()
  basename = os.path.basename(outfile)
  stats[basename] = tuple([[] for _ in range(len(PROGRAMS))])
  
  # run compalignp
  refalignment = infile.replace("unaligned", "structural-no_str")
  for (i, file) in [(x, outfile + str(x) + ".aln") for x in PROGRAMS]:
    proc = subprocess.Popen([compali_bin, "-t", file, "-r", refalignment],\
           bufsize=-1, executable=compali_bin, stdout=subprocess.PIPE, shell=False)
    stats[basename][i].append(float(proc.communicate()[0]))
  
  # extract statistics from files
  for (i, file) in [(x, outfile + str(x) + ".stat") for x in PROGRAMS]:
    (sci, mpi) = (None, None)
    for line in open(file, 'r'):
      if line.startswith(" Mean pairwise identity:"):
        mpi = float(re.findall(r"(-?\d+\.\d*)", line)[0])
      if line.startswith(" Structure conservation index:"):
        sci = float(re.findall(r"(-?\d+\.\d*)", line)[0])
    
    if sci == None or mpi == None:
      print >>sys.stderr, "Could not parse file " + file
      exit(1)
      
    stats[basename][i].extend([sci, mpi])
sys.stderr.write('\n')

###############
##  RESULTS  ##
###############

# provide data sorted by MPI
data = [zip(*x) for x in zip(*stats.values())]
for i in PROGRAMS:
  (data[i][MPI], data[i][SPS], data[i][SCI]) = (list(x) for x in zip(*sorted(zip(data[i][MPI], data[i][SPS], data[i][SCI]))))
  
# print cumulative values
def _make_stat_str(values):
  return " \tMean "   + str(tuple([round(numpy.mean(x),2) for x in values]))\
       + " \tMedian " + str(tuple([round(numpy.median(x),2) for x in values]))\
       + " \tStdDev " + str(tuple([round(numpy.std(x),2) for x in values]))
       
# write to file
tfile = open(os.path.join(results_dir, "time.txt"), "a")
tfile.write("######## ")
tfile.write("RUN on " + datetime.now().strftime("%A, %d. %B %Y %H:%M") + '\n\n')
print "Values for (Sum of Pairs Score, Structure Conservation Index, Mean Pairwise Identity):"
for (pname,pval) in zip(PROGLBL,data):
  print "  " + pname + _make_stat_str(pval)
  tfile.write("  " + pname + _make_stat_str(pval) + '\n')

tfile.write('\nTotal run times:\n')
for (pname, tm) in zip(PROGLBL, progtime)[1:]:
  tfile.write('  ' + pname + ' \t{} seconds.'.format(tm) + '\n')
tfile.write('\n\n')
tfile.close()

# set view area and file for plots
view = ([30, 100, 0.2, 1], [30, 100, 0, 1.4])
pdf = PdfPages(os.path.join(results_dir, "benchmark.pdf"))

for score in (SPS, SCI):
  # calculate lowess function
  x = numpy.array(data[REF][MPI], numpy.float)
  y = [prog[score] for prog in data]
  f = [lowess(x, numpy.array(y[i], numpy.float)) for i in PROGRAMS]
  c = ['green','black','blue','red','magenta','cyan','yellow','cyan','magenta']
   
  # plot data
  fig = plt.figure()
  for prog in (0,5,4,1,2,3,6,7):
    plt.plot(map(int, x), y[prog], ".", ms=3, color=c[prog])
    plt.plot(x, f[prog], label=PROGLBL[prog], color=c[prog])
    
  plt.axis(view[score])
  plt.xlabel(SCORELBL[MPI])
  plt.ylabel(SCORELBL[score])
  plt.legend(loc=4, fontsize='x-small', borderaxespad=0.) # bbox_to_anchor=(1.05, 1), 
  pdf.savefig(fig)
  plt.clf()

print "Plot saved to results/benchmark.pdf"
pdf.close()

