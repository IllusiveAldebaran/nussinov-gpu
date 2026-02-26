
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
