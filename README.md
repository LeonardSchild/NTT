# NTT
C++ and CUDA header only implementation of the number theoretic transform and
associated polynomial product over the cyclotomic ring <a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{Z}_Q[X]/(X^N&space;&plus;&space;1),\,\,\,\,\,\,\&space;N=2^k" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{Z}_Q[X]/(X^N&space;&plus;&space;1),\,\,\,\,\,\,\&space;N=2^k" title="\mathbb{Z}_Q[X]/(X^N + 1),\,\,\,\,\,\,\ N=2^k" /></a>

Source https://eprint.iacr.org/2014/646.pdf

# Timings

N = 1024, Q = 134215681, T = uint32_t


g++ | clang++ | nvcc & g++ | nvcc & clang++
----|---------|------------|----------------
55us | 52us | TBA | TBA
