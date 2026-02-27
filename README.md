# Build Instructions

Compile and Run with only 50 nucleotides of input file

```bash
make SEQ_N=50
```


Build and run entire sequence:
```bash
make SEQ_N=$(wc -c inputs/ec16s.seq | cut -d ' ' -f 1)
```


### Instructions

Run as:
```bash
./nussinov <SEQ>
```


Or on DSMLP:
```
/opt/launch-sh/bin/launch.sh -v a30 -c 8 -g 1 -m 8 -i yatisht/ece213-wi26:latest -f ./nussinov-gpu/run-commands.sh
```
