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
./fcgr <in-kmc-list>
```

### Example
```
cd example
KMC=../build/kmc-prefix/src/kmc/bin/kmc
for fa in *.fa ; do $KMC -v -m4 -sm -b -ci0 -cs65535 -k6 -fm $fa $fa . ; done
ls $PWD/*.fa > list
../fcgr list
ls *.npy # <- these are the FCGR
```