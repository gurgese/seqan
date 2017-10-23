# laragu
The new version of Lara tool for the structural-alignment of RNA sequences 

To test the T-Coffee library preparation computed with the FIXEDINTER mode
-i ./apps/lara2/demo/rna.fasta -t 2 -a -lgo -2 -lge -1 -lbm 1 -ssc 10 -tcm 3 -tb 0.2 -v 1 -td ./apps/lara2/demo/tmp -iter 100 

To test the tcoffe lib file production (A problem with the multithread execution is present and should be fixed)
-i ./apps/lara2/demo/input.fa -t 2 -a -lsm ./apps/lara2/ribosum_matrices/RIBOSUM45-30_N.txt -lgo -2 -lge -1 -lbm 1 -ssc 10 -tcm 2 -tb 0.2 -v 1 -t 1 -td /home/vitrusky8/git/seqan/apps/lara2/demo/tmp

To test 200 iterations over a small dataset with lower bound MWM and a scaling factor applied on the ribosum matrices
-i ./apps/lara2/demo/input_small.fa -t 2 -v 2 -a -lsm ./apps/lara2/ribosum_matrices/RIBOSUM45-30_N.txt -lgo -2 -lge -1 -lbm 1 -ssc 10

To test 200 iterations over a small dataset with lower bound MWM
-i ./apps/lara2/demo/U5.fa -v 2 -a -lsm ./apps/lara2/ribosum_matrices/RIBOSUM45-30_N.txt -lgo -2 -lge -1 -lbm 0 -iter 200 -tb 1e-15

To test 200 iterations over a small dataset with lower bound approximation
-i ./apps/lara2/demo/U5.fa -v 2 -a -lsm ./apps/lara2/ribosum_matrices/RIBOSUM45-30_N.txt -lgo -2 -lge -1 -lbm 1 -iter 200 -tb 1e-15

To test 100 iterations over a small dataset
-i ./apps/lara2/demo/input_small.fa -t 2 -v 3 -a -lsm ./apps/lara2/ribosum_matrices/RIBOSUM45-30_N.txt -lgo -2 -lge -1 -lbm 1 -iter 100

To test the usage of a small alignment useful for the lower and upper bound checking (-lbm 1 set the algorithm to compute an approximation of the MWM that should be deeply tested)
-i ./apps/lara2/demo/input_small.fa -t 2 -v 3 -a -lsm ./apps/lara2/ribosum_matrices/RIBOSUM45-30_N.txt -lgo -2 -lge -1 -lbm 1

To test the creation of the lamb string that will store the aligned interections that will be updated at each iteration
-i ./apps/lara2/demo/input.fa -ir ./apps/lara2/demo/input.fa -t 2 -v 3 -a -lsm ./apps/lara2/ribosum_matrices/RIBOSUM45-30_N.txt -lgo -2 -lge -1

To test the first local alignment
-i ./apps/lara2/demo/input.fa -ir ./apps/lara2/demo/input.fa -t 2 -v 3 -a

To run with the acquisition of ribosum matrices
-i ./apps/lara2/demo/input.fa -ir ./apps/lara2/demo/input.fa -t 2 -v 1 -lsm ./apps/lara2/ribosum_matrices/RIBOSUM45-30_N.txt

To run with a single fasta file run the following command
-i ./apps/lara2/demo/input.fa -t 2 -v 2

To run with two fasta files run the following command
-i ./apps/lara2/demo/input.fa -ir ./apps/lara2/demo/input.fa -t 2 -v 2

To test the dbn files run the following command
-i ./apps/lara2/demo/input.fa -ir ./apps/lara2/demo/rnafold.dbn -t 2 -v 2

To test the gap open/extension that produce the desired configuration
-i ./apps/lara2/demo/input_jorg.fasta -t 1 -lbm 0 -v 2 -td ./apps/lara2/demo/tmp -iter 5  -stsc 3 -ggo -0.3 -gge 0.1


Configuration of g++ version to be used in the INSTALL instructions

These commands are based on a askubuntu answer http://askubuntu.com/a/581497
To install gcc-6 (gcc-6.2.0), I had to do more stuff as shown below.
USE THOSE COMMANDS AT YOUR OWN RISK. I SHALL NOT BE RESPONSIBLE FOR ANYTHING.
ABSOLUTELY NO WARRANTY.

sudo apt-get update && \
sudo apt-get install build-essential software-properties-common -y && \
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y && \
sudo apt-get update && \
sudo apt-get install gcc-snapshot -y && \
sudo apt-get update && \
sudo apt-get install gcc-6 g++-6 -y && \
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6 60 --slave /usr/bin/g++ g++ /usr/bin/g++-6 && \
sudo apt-get install gcc-4.8 g++-4.8 -y && \
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 60 --slave /usr/bin/g++ g++ /usr/bin/g++-4.8;

When completed, you must change to the gcc you want to work with by default
sudo update-alternatives --config gcc

To verify if it worked. Just type in your terminal
gcc -v

If everything went fine you should see gcc 6.2.0
