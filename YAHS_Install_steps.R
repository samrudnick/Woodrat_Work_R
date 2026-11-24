# This is how I created an executable wrapper for YAHS:

#Step 1: Clone the repository from github:
git clone https://github.com/c-zhou/yahs

#Step 2: In the directory downloaded called "yahs" run the following dependencies:
module purge
module load gcc/11.2.0
module load zlib/1.2.11

make clean
make CC=gcc

#Step 3: Test it
./yahs --help
ldd ./yahs

#You can move the yahs file to other directories to run jobs mannually but you must give it permission to run using this:
chmod u+x yahs