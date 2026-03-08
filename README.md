# CompPhysHack 2026

Benchmarks for the numerical approaches of one-dimensional quantum mechanical problems

## Hamiltonian

```math
\hat{H} = - \frac{1}{2} \frac{\mathrm{d}^2}{\mathrm{d}x^2} + \frac{1}{2}x^2
```

## FDM v.s. QTT + TCI

- FDM: Descritize with Finite diffential approximation
- QTT + TCI: 

Trial wave function:
```math
\psi(x) = \exp \left( -\frac{1}{2}x^2 \right)
```

Exact energy:
```math
E = \frac{\langle \psi_0 | \hat{H} | \psi_0 \rangle}{\langle \psi_0 | \psi_0 \rangle} = \frac{1}{2}
```

Numerical energy:

|  c |        n |               FDM |           QTT+TCI |
| -- | -------- | ----------------- | ----------------- |
|  6 |       64 | 0.499131263832253 | 0.696303556584334 |
|  7 |      128 | 0.481589657154101 | 0.481589657154243 |
|  8 |      256 | 0.495255158330897 | 0.495255158330209 |
|  9 |      512 | 0.498807047001768 | 0.498807047003291 |
| 10 |     1024 | 0.499701631478844 | 0.499701631481162 |
| 11 |     2048 | 0.499925436210840 | 0.499925436208489 |
| 12 |     4096 | 0.499981365376262 | 0.499981365378351 |
| 13 |     8192 | 0.499995342307949 | 0.499995342303344 |
| 14 |    16384 | 0.499998835708073 | 0.499998835706305 |
| 15 |    32768 | 0.499999708943988 | 0.499999708961086 |
| 16 |    65536 | 0.499999927238317 | 0.499999927238436 |
| 17 |   131072 | 0.499999981805869 | 0.499999981850532 |
| 18 |   262144 | 0.499999995457983 | 0.499999995434375 |
| 19 |   524288 | 0.499999998810651 | 0.499999999145209 |
| 20 |  1048576 | 0.499999999724325 | 0.499999995042597 |

Calculation Time:

|  c |        n |             FDM |         QTT+TCI |
| -- | -------- | --------------- | --------------- |
|  6 |       64 |       22.541 μs |        1.645 ms |
|  7 |      128 |       34.369 μs |        2.415 ms |
|  8 |      256 |       74.828 μs |        3.593 ms |
|  9 |      512 |       93.436 μs |        5.042 ms |
| 10 |     1024 |      180.841 μs |        6.342 ms |
| 11 |     2048 |      361.011 μs |       11.272 ms |
| 12 |     4096 |      725.843 μs |       11.710 ms |
| 13 |     8192 |        1.261 ms |       10.414 ms |
| 14 |    16384 |        3.093 ms |       13.150 ms |
| 15 |    32768 |        6.205 ms |       13.704 ms |
| 16 |    65536 |       13.965 ms |       13.348 ms |
| 17 |   131072 |       31.262 ms |       14.062 ms |
| 18 |   262144 |       59.420 ms |       24.675 ms |
| 19 |   524288 |      178.922 ms |       16.491 ms |
| 20 |  1048576 |      400.168 ms |       15.956 ms |

## How to Run

```sh
git clone https://github.com/ohno/CompPhysHack2026Ohno.git
cd CompPhysHack2026Ohno
julia --project=./julia --startup-file=no -e 'import Pkg; Pkg.instantiate()'
julia --project=./julia --startup-file=no benchmark.jl
```