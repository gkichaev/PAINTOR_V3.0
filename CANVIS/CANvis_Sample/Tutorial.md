# INSTRUCTIONS

1. Install additional libraries by running `install_mac.sh` or `install_linux.sh`. (see README for other dependencies)

```
$ ./install_mac.sh
```
2.
Run the following command inside of the CANVIS_sample folder 
```
$ python CANVIS.py -l chr4.3473139.rs6831256.post.filt.300 -z tg.Zscore -r chr4.3473139.rs6831256.ld.filt.300 -a chr4.3473139.rs6831256.annot.filt.300 -s E066.H3K27ac.narrowPeak.Adult_Liver E066.H3K4me1.narrowPeak.Adult_Liver -t 99 -i 3381705 3507346
```
3.
The image produced should match the `example.pdf` file given in the sample folder. 