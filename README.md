# fcgr

### Installation
```
git clone https://github.com/pg-space/fcgr.git
cd fcgr
mkdir build ; cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..
```

### Usage
`fcgr` takes as input the list of KMC databases (prefix used while running KMC).
For each of them, it outputs the FCGR in `.npy` format.

```
./fcgr [-m MASK] <in-kmc-list>
```

`MASK` is a bit string indicating which bases of a kmer to keep. For instance, if `MASK=10011` and kmer is `AACGT`, kmer `AGT` will be stored in the FCGR. By default it's a k-long string of 1s.

**Note:** if input kmers are in canonical form, masked kmers won't be in canonical form anymore.

### Example
```
cd example
KMC=../build/kmc-prefix/src/kmc/bin/kmc
for fa in *.fa ; do $KMC -v -m4 -sm -b -ci0 -cs65535 -k6 -fm $fa $fa . ; done
ls $PWD/*.fa > list
../fcgr list
ls *.npy # <- these are the FCGR
```
