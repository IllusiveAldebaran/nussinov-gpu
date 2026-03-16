### Create and Use Synthetic RNA Sequences

First generate a list of rna sequences. Below we generate 128 sequences of 256 length.
```
cd synthetic
./synthetic.py -l 256 -n 128
cd ..
```

### Instructions
1. If running it on DSMLP:

```
/opt/launch-sh/bin/launch.sh -v a30 -c 8 -g 1 -m 8 -i yatisht/ece213-wi26:latest -f ./nussinov-gpu/run-commands.sh
```
Make sure to match SEQ_N parameter in run-commands.sh and the length of the synthetic sequence you just generated

2. If using your own GPU:
```bash
./nussinov <seq_file_list.txt>
```

### Parameters

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

