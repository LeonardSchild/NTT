# Entity
C++ and CUDA header only implementation of the number theoretic transform and
associated polynomial product over the cyclotomic ring <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{Z}_Q[X]/(X^N&space;&plus;&space;1),\,\,\,\,\,\,\&space;N=2^k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{Z}_Q[X]/(X^N&space;&plus;&space;1),\,\,\,\,\,\,\&space;N=2^k" title="\mathbb{Z}_Q[X]/(X^N + 1),\,\,\,\,\,\,\ N=2^k" /></a>

Source https://eprint.iacr.org/2014/646.pdf

# Some Timings

All timings are in *microseconds*

### Fast Engine (32 bit)

| N | Q | T | 
|---| ----|----|
| 1024 | 134215681 | uint32_t |

|  | g++ | clang++ | nvcc & g++ | nvcc & clang++
|---| ----|---------|------------|----------------
| Forward transform | 15 | 17 | TBA | TBA
| Backward transform | 17| 16 | TBA | TBA
| Multiplication | 54 | 49| TBA | TBA