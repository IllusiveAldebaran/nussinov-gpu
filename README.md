# Build Instructions

Compile and Run with different parameters: SEQ_N, LOOP, BLOCK, and GRID
(SEQ_N = 30, BLOCK = 64, LOOP = 4, GRID = 1 by default)

```bash
make SEQ_N=50 BLOCK=512 GRID=8
```

Run on CPU (Run on GPU after CPU)
```bash
make CPU=1
```

Enable debug flags.
```bash
make DEBUG=1 CPU=1
```


### Instructions

Run as:
```bash
./nussinov <seq_file_list.txt>
```


Or on DSMLP:
```
/opt/launch-sh/bin/launch.sh -v a30 -c 8 -g 1 -m 8 -i yatisht/ece213-wi26:latest -f ./nussinov-gpu/run-commands.sh
```


# Create and Use Synthetic RNA Sequences


First generate a list of rna sequences. Below we generate 128 sequences of 256 length.
```
cd synthetic
./synthetic.py -l 256 -n 128
cd ..
```

```
make clean
make SEQ_N=256
./nussinov synthetic/synthetic_sequences.txt
```
