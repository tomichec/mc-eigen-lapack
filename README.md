# mc-eigen-lapack

## Description

This program computes eigenvalues of a transition rates matrix
corresponding to a Markov chain model and saves them into the files. 

 * `vm_INa_LAPACK.dat` : grid of membrane voltage
 * `evals_INa_LAPACK.dat` : eigenvalues
 * `left_evecs_INa_LAPACK.dat` : left eigenvector matrix
 * `right_evecs_INa_LAPACK.dat` : right eigenvector matrix

Those
files can be then used for the exponential integration from a cellular
model of the Markov Chain matrix using Matrix Rush-Larsen method.

## Dependencies

This program requires [C compiler](http://www.gnu.org/software/gcc/)
and [LAPACK](http://www.netlib.org/lapack/) with the LAPACKe API for
C.

Optionally [Make](https://www.gnu.org/software/make/) for a simple
compilation by `make`.


## Compilation

The program is compiled by a command

```
gcc -Wall -I/usr/include/lapacke/ -llapacke -llapack -lm inaEigenLAPACK -o inaEigenLAPACK
```

## Running

Run the program by a command: ` ./inaEigenLAPACK ` which creates files
in the same directory where it is run. You can refer to those files
from the cellular model.

## License

Free software under [GNU GPLv3](http://www.gnu.org/licenses/gpl.txt).
