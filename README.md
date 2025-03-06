# Gross lattices of supersingular elliptic curves

This repository contains the implementation of the algorithms and the results of the computations described in the paper.

**arXiv:** [2503.03478v1](https://arxiv.org/abs/2503.03478)

**Authors:** Chenfeng He, Gaurish Korpal, Ha Tran, and Christelle Vincent

## Repository Organization

### A. Functions

This directory contains all the functions that can be called to run various experiments.

```
functions/
├── __init__.py
├── fp/
│   ├── __init__.py
│   ├── gross/
│   │   ├── __init__.py
│   │   ├── numeric.py
│   │   ├── symbolic.py
│   │   └── type.py
│   └── ibukiyama/
│       ├── __init__.py
│       └── maximal.py
└── fp2/
    ├── __init__.py
    ├── deuring/
    │   ├── __init__.py
    │   ├── correspondence.py
    │   ├── cost.py
    │   ├── klpt.py
    │   └── xonly.py
    └── eichler/
        ├── __init__.py
        ├── bounds.py
        └── maximal.py
```

The `__init__.py` (empty) files are required to make Python treat directories containing the file as packages. The rest of the files in each sub-directory are described below.

#### A.1 `fp`
Contains all functions that work only for supersingular elliptic curves with j-invariant lying in $\mathbb{F}_p$.
- `gross`: Contains the algorithms from our paper.
    - `numeric.py` and `symbolic.py` contain SageMath implementations of our algorithm that let one compute the Gram matrix of the minimal Gross lattice given the first successive minima $D_1$ value.
    - `type.py` contains a function to classify and extract parameters from the Gram matrix of the minimal Gross lattice into one of the four types described in the paper for the 13 CM curves over $\mathbb{Q}$.
- `ibukiyama`: `maximal.py` contains a SageMath implementation of the algorithm from [Li-Ouyang-Xu](http://staff.ustc.edu.cn/~yiouyang/curve_over_Fp_1101.pdf) that lets us compute Ibukiyama's maximal orders $O(q,r)$ and $O'(q,r)$. It also contains other helper functions to compute the Gross lattice and Gram matrix.

#### A.2 `fp2`
Contains all functions that work for all supersingular elliptic curves with j-invariant lying in $\mathbb{F}_{p^2}$.
- `deuring`: Contains a copy of the code by [Eriksen-Panny-Sotáková-Veroni](https://eprint.iacr.org/2023/106) and the necessary intra-package references. This lets us find the j-invariant corresponding to a given maximal order in quaternion algebra $B_{p}$.
- `eichler`: An add-on to SageMath's Quaternion Algebra functionality.
    - `maximal.py` utilizes SageMath's interface to Magma so that we can use the algorithm from [Kirschmer-Voight](https://arxiv.org/abs/0808.3833) to compute all the maximal orders of $B_{p}$.
    - `bounds.py` contains an implementation of the bounds on the third successive minima $D_3$ of the Gross lattice.

### B. Scripts

This directory contains scripts that implement various experiments.

```
scripts/
├── CM_Gross/
│   ├── cm.sage
│   └── ref/
├── CM_Gross_finite_cases/
│   ├── cm_finite.sage
│   └── ref/
├── family_20/
│   ├── fam20.sage
│   └── ref/
├── family_20_2mod3/
│   ├── fam2mod3.sage
│   └── ref/
├── finite_cases/
│   ├── finite.sage
│   └── ref/
├── NE_values/
│   ├── NEvalue.sage
│   └── ref/
├── toy_Gross/
│   ├── csidh419gross.sage
│   └── ref/
├── toy_Ibukiyama/
│    ├── csidh419.sage
│    └── ref/
└── toy_Kaneko/
     ├── csidh419kaneko.sage
     └── ref/
```

The sub-directory `ref/` contains the reference outputs of these scripts. The rest of the files in each sub-directory are described below.

#### B.1 `CM_Gross`
Contains the script that symbolically computes the Gram matrix of the minimal Gross lattice of all 13 CM curves over $\mathbb{Q}$.
- One can run the script by navigating to this directory and then executing `$ sage cm.sage <d>` where `<d>` is the absolute value of the CM discriminant of one of the 13 CM curves over $\mathbb{Q}$, namely $\{3,4,7,8,11,12,16,19,27,28,43,67,163\}$.
- The `ref/` directory contains outputs for all 13 values of `d`. For example, `$ sage cm.sage 163` will generate `cm_163.txt` that contains all 81 possible Gram matrices.

#### B.2 `CM_Gross_finite_cases`
Contains the script that computes the LLL-reduced Gram matrix of the Gross lattice for the CM curves over $\mathbb{Q}$ whose $N_E$ is smaller than 37.
- One can run the script by navigating to this directory and then executing `$ sage cm_finite.sage <d>` where `<d>` is the absolute value of the CM discriminant of one of the 5 CM curves over $\mathbb{Q}$, namely $\{3,4,7,8,11\}$.
- The `ref/` directory contains outputs for all 5 values of `d`. For example, `$ sage cm_finite.sage 8` will generate `cases_8.txt` that contains all 3 Gram matrices of interest.

#### B.3 `family_20`
Contains the script that computes LLL-reduced Gram matrices of all maximal orders with the first successive minima equal to 20 and lying in $B_{p}$ for $p \equiv c \pmod {20}$ lying between M and N.
- One can run the script by navigating to this directory and then executing `$ sage fam20.sage <c> <M> <N>` where `<c>` is $\{1,3,5,\ldots,17,19\}$ and $0<M<N$.
- The `ref/` directory contains outputs for all 10 values of `c`. For example, `$ sage fam20.sage 19 1 500` will generate `fam20_19_1_500.txt` that contains all 10 Gram matrices of interest. In particular, from `fam20_11_1_500.txt` and `fam20_19_1_500.txt`, one can observe that for the prime families like $p \equiv 11 \pmod{20}$ and $p \equiv 19 \pmod{20}$, all curves lie over $\mathbb{F}_p$ and the Gram matrices are one of the four types.

#### B.4 `family_20_2mod3`
Contains the script that computes the j-invariants of all maximal orders with the first successive minima equal to 20 and lying in $B_{p}$ for $p \equiv c \pmod {20}$ and $p \equiv 2 \pmod 3$ lying between M and N.
- One can run the script by navigating to this directory and then executing `$ sage fam2mod3.sage <c> <M> <N>` where `<c>` is $\{13,17\}$ and $0<M<N$.
- The `ref/` directory contains outputs for all 10 values of `c`. For example, `$ sage fam2mod3.sage 13 20 1000` will generate `fam2mod3_13_20_1000.txt` that contains all 10 Gram matrices of interest.

#### B.5 `finite_cases`
Contains the script that computes LLL-reduced Gram matrices of all maximal orders with the first successive minima equal to $D_1\leq 2*p^{2/3}$ and lying in $B_{p}$ for $p$ lying between M and N.
- One can run the script by navigating to this directory and then executing `$ sage finite.sage <M> <N>` where $0<M<N$.
- The `ref/` directory contains output for finite cases to verify $D_3$ bound. For example, `$ sage finite.sage 1 100` will generate `cases_1_100.txt` that contains all Gram matrices of interest. In particular, from `cases_1_100.txt` one can observe that for the primes up to 100, $D_3 < \frac{3}{5}p+5$ for curves not defined over $\mathbb{F}_p$.

#### B.6 `NE_values`
Contains the script to estimate $N_E$, the smallest supersingular prime (greater than 3) such that for $p \geq N_E$ we get $D_1 = d$ by looking for a continuous sequence of 10 LLL-reduced Gram matrices.
- One can run the script by navigating to this directory and then executing `$ sage NEvalue.sage <d>` where `<d>` is the absolute value of the CM discriminant of one of the 13 CM curves over $\mathbb{Q}$, namely $\{3,4,7,8,11,12,16,19,27,28,43,67,163\}$.
- The `ref/` directory contains outputs for all 13 values of `d`. For example, `$ sage NEvalue.sage 27` will generate `NE_27.txt` that contains the $N_E$ value.

#### B.7 `toy_Gross`
Contains the script to numerically compute the Gram matrix of the minimal Gross lattice of all elliptic curves given $D_1$ values $\{3, 4, 7, 12, 15, 16, 20, 23, 27, 28, 35, 36, 39\}$ obtained from `toy_Ibukiyama/` and test if they belong to the isogeny graph of CSIDH for $p=419>37$.
- One can run the script by navigating to this directory and then executing `$ sage csidh419gross.sage`.
- The `ref/` directory contains the output `csidh419gross.txt` with all the Gram matrices. From this, one can observe that for the fixed prime $p=419$, if the curve lies over $\mathbb{F}_p$, then the Gram matrices are one of the four types.

#### B.8 `toy_Ibukiyama`
Contains the script that computes LLL-reduced Gram matrices of all elliptic curves belonging to the isogeny graph of CSIDH for $p=419$.
- One can run the script by navigating to this directory and then executing `$ sage csidh419.sage`.
- The `ref/` directory contains output `csidh419.txt` with all the maximal orders and Gram matrices. From this, one can observe the $D_1$ values for all curves over $\mathbb{F}_p$ and with endomorphism ring of type $O(q,r)$.

#### B.9 `toy_Kaneko`
Contains the script to numerically compute the Gram matrix of the minimal Gross lattice of all elliptic curves given $D_1\leq \left\lceil 4\sqrt{\frac{p}{3}}\right\rceil$ and test if they belong to the isogeny graph of CSIDH for $p=419>37$.
- One can run the script by navigating to this directory and then executing `$ sage csidh419kaneko.sage`.
- The `ref/` directory contains the output `csidh419kaneko.txt` with all the Gram matrices. From this, one can observe that for the fixed prime $p=419$, we get minimal Gross lattice for $D_1 \in \{3, 4, 7, 12, 15, 16, 20, 23, 27, 28, 35, 36, 39, 43, 47\}$ which matches with the $D_1$ values obtained from `toy_Ibukiyama/`, except for $43$ and $47$.


## System Requirements

The code was written and tested on a system with the following specifications:

**CPU:** 12th Gen Intel i7-12800HX (24)

**Memory:** 15843MiB

**OS:** Ubuntu 22.04.5 LTS on Windows 11 x86_64 (WSL)

**Software:** SageMath version 10.3 and Magma V2.28-8
