# fcgr

### Installation
```
git clone https://github.com/pg-space/fcgr.git
cd fcgr
mkdir build ; cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
./fcgr -h
```

### Usage

#### Single mode
`fcgr` takes as input a list of KMC databases (prefix used while running KMC).
For each of them, it outputs the FCGR in `.npy` format.

```
./fcgr single [-h] <in-kmc-list>
```

`MASK` is a bit string indicating which bases of a kmer to keep. For instance, if `MASK=10011` and kmer is `AACGT`, kmer `AGT` will be stored in the FCGR. By default it's a k-long string of 1s.

**Note:** if input kmers are in canonical form, masked kmers won't be in canonical form anymore.

##### Example
```
cd example
KMC=../build/kmc-prefix/src/kmc/bin/kmc
for fa in 1.fa 2.fa 3.fa ; do $KMC -v -m4 -sm -b -ci0 -cs65535 -k6 -fm $fa $fa.db . ; done
echo $PWD/1.fa.db $PWD/2.fa.db $PWD/3.fa.db | tr " " "\n" > kmc.list
../fcgr single -o single-output kmc.list
ls single-output/*.npy # <- these are the FCGR
```

#### Multi mode
`fcgr` takes as input a FASTA file and, for each record, outputs the FCGR in `.npy` format.

```
./fcgr multi [-h] <fasta>
```

##### Example
```
cd example
../fcgr multi -o multi-output 123.fa
ls multi-output/*.npy # <- these are the FCGR
```
