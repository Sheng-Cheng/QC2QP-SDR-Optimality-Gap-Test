# QC2QP-Semidefinite-Relaxation-Optimality-Gap-Test
This is a MATLAB script for testing whether a quadratic program with two quadratic constraints (QC2QP) can be solved using its semidefinite relaxation, in which case we say there is no optimality gap. The details of this test and the methematical proof behind can be found here:<br />
[An Optimality Gap Test for a Semidefinite Relaxation of a Quadratic Program with Two Quadratic Constraints](https://arxiv.org/abs/1907.02989)<br />

If you think this test is helpful to your research/application, please cite:<br />
@article{cheng2019optimalitygap,<br />
  title={An Optimality Gap Test for a Semidefinite Relaxation of a Quadratic Program with Two Quadratic Constraints},<br />
  author={Cheng, Sheng and Martins, Nuno C},<br />
  journal={arXiv preprint arXiv:1907.02989},<br />
  year={2019}<br />
}

## How to use this test

Just download the code and run [QC2QP_SDR_optimalityGap_test.m](https://github.com/Sheng-Cheng/QC2QP-SDR-Optimality-Gap-Test/blob/master/QC2QP_SDR_optimalityGap_test.m) in MATLAB.

## Prerequisites

You need [CVX](http://cvxr.com/cvx/) and [YALMIP](https://yalmip.github.io/) to run this test.



## Running the test
The QC2QP problem has the following form:

![alt text][image1]

[image1]: https://github.com/Sheng-Cheng/QC2QP-SDR-Optimality-Gap-Test/blob/master/images/optimization%20problem.png "QC2QP problem"

To run the optimality gap test, all the coefficients in the quadratic objective function and quadratic constraints are written in the homogeneous form, i.e., 

![alt text][image2]

[image2]: https://github.com/Sheng-Cheng/QC2QP-SDR-Optimality-Gap-Test/blob/master/images/homogeneous%20quadratic%20form.png "Homogeneous quadratic form"

Note that the test requires **Slater's condition to hold for the (primal) problem and its dual problem**. If the problem data does not meet Slater's condition for either the (primal) problem or the dual problem, then the test will abort and throw an alert.


### Syntax

`status = QC2QP_SDR_optimalityGap_test(M0,M1,M2)` <br />
`status = QC2QP_SDR_optimalityGap_test(M0,M1,M2,epsilon2)` <br />
`status = QC2QP_SDR_optimalityGap_test(M0,M1,M2,epsilon2,verbosity)` <br />

### Description
`status = QC2QP_SDR_optimalityGap_test(M0,M1,M2)` returns the status whether the QC2QP (specified by the problem data `M0`, `M1`, and `M2`) has no optimality gap. The output is a struct that contains fields regarding the solutions of the semi-definite relaxation and the Lagrange dual of the QC2QP (see details in the section 'Output arguments'). 

`status = QC2QP_SDR_optimalityGap_test(M0,M1,M2,epsilon2)` defines the tolerance to use in the rank computation. For example, the rank of a matrix A is the number of singular values of A that are greater than `epsilon2`. 

`status = QC2QP_SDR_optimalityGap_test(M0,M1,M2,epsilon2,verbosity)` selects the verbosity of the function. When `verbosity` is set to 1, complete output messages from the solvers CVX and fmincon will be displayed in the command window. When `verbosity` is set to 0, all the solver messages are muted. 

### Example

Consider the following problem:

![alt text][image3]

[image3]: https://github.com/Sheng-Cheng/QC2QP-SDR-Optimality-Gap-Test/blob/master/images/Example%20problem.png "Example problem"

In this case, the input arguments `M0`, `M1`, and `M2` are in the homogeneous quadratic form:

![alt text][image4]

[image4]: https://github.com/Sheng-Cheng/QC2QP-SDR-Optimality-Gap-Test/blob/master/images/Example%20homogeneous%20data.png "Example data"

Run the test with default tolerance `epsilon2` = 1e-5 and the verbosity option, <br />`status = QC2QP_SDR_optimalityGap_test(M0,M1,M2,[],1);` <br /> we get the following messages:

![alt text][image5]

[image5]: https://github.com/Sheng-Cheng/QC2QP-SDR-Optimality-Gap-Test/blob/master/images/MATLAB%20output%20message.png "MATLAB output message"

We can tell that there is no optimality gap by examing the value of `status.noOptimalityGap` which is 1 in this case. To verify this result, we can examine that the primal optimal value `status.primal_optval` and the optimal value of the semidefinite relaxation `status.SP_optval` both equal -9.603839.

Try the [QC2QP\_SDR_optimalityGap\_test\_example.m](https://github.com/Sheng-Cheng/QC2QP-SDR-Optimality-Gap-Test/blob/master/QC2QP_SDR_optimalityGap_test_example.m) to get familiar with the test. In this script, we provide 10 groups of sample problem data for users to get a taste. Just uncomment any group of data you want to try and hit 'run'.

### Input arguments
`M0` &mdash; coefficients of the objective function in the homogeneous form. The matrix `M0` must be real and symmetric.

`M1` and `M2` &mdash; coefficients of the constraint functions in the homogeneous form. The matrices `M1` and `M2` must be real and symmetric.

`epsilon2` &mdash; tolerance for the rank computation. The default value is 1e-5.

`verbosity` &mdash; amount of messages displayed during the test. Allowed value is either 0 or 1. The default value is 0.

### Output argument
The output is a struct which contains the following fields:

* `epsilon2` &mdash; tolerance for the rank computation
* `SlaterSP` &mdash; indicator for whether Slater's condition holds for the semidefinite relaxation (SP). 1 for holding and 0 for not.
* `SlaterSD` &mdash; indicator for whether Slater's condition holds for the dual problem (SD). 1 for holding and 0 for not.
* `noOptimalityGap` &mdash; indicator of whether there exists an optimality gap. 1 for no optimality gap, 0 for optimality gap, and NaN for undetermined.
* `exception.code` &mdash; indicator for whether an exception appears:
  * `0` for no exception;
  * `1` for problem ill-posed;
  * `2` for rank(X_hat) >= 3;
  * `3` for error showing up when rank(Z_hat) < n-1 while rank(X_hat) == 2
  * `4` for alpha1 undetermined
  * Please email the author (cheng@terpmail.umd.edu) about any exception (or error) if the exception code is not 0
* `exception.M0` &mdash; value of `M0` that causes an exception. It only shows up when `exception.code` is not 0.
* `exception.M1` &mdash; value of `M1` that causes an exception. It only shows up when `exception.code` is not 0.
* `exception.M2` &mdash; value of `M2` that causes an exception. It only shows up when `exception.code` is not 0.

If the problem data yield Slater's condition holding for both (SP) and (SD), then the following fields are included.

* `SP_optval` &mdash; optimal value of (SP)
* `SP_optsol` &mdash; optimal solution of (SP), X\_hat
* `SP_x1` &mdash; rank-one decomposition of X\_hat
* `SP_x2` &mdash; rank-one decomposition of X\_hat (this field only exists if rank(X\_hat) == 2)
* `SD_optval` &mdash; optimal value of (SD)
* `SD_optsol.Z` &mdash; optimal solution of (SD), Z\_hat
* `SD_optsol.y0` &mdash; optimal solution of (SD), y\_hat0
* `SD_optsol.y1` &mdash; optimal solution of (SD), y\_hat1
* `SD_optsol.y2` &mdash; optimal solution of (SD), y\_hat2
* `rank_X` &mdash; rank of the solution to (SP) with tolerance epsilon2
* `rank_Z` &mdash; rank of the solution to (SD) with tolerance epsilon2

If the problem data yield no optimality gap, then the following fields are included.

* `primal_optval` &mdash; optimal value of the primal problem
* `primal_optsol` &mdash; optimal solution of the primal problem

---

# Dataset

The dataset contains QC2QP instances of dimension n from 2 to 7. Each dimension contains 1000 randomly generated feasible nonconvex QC2QP instances.

## Explanation

The dataset is stored in the file [QC2QP_dataset.mat](https://github.com/Sheng-Cheng/QC2QP-SDR-Optimality-Gap-Test/blob/master/QC2QP_dataset.mat), which contains the variable QC2QP_dataset. Under QC2QP_dataset, there are 5 fields corresponding to the instances in each dimension.

![alt text][image6]

[image6]: https://github.com/Sheng-Cheng/QC2QP-SDR-Optimality-Gap-Test/blob/master/images/field_dim.png "QC2QP_dataset"

We take the two-dimensional instances 'dim2' for illustration. Under this field, there are 
* `instance` &mdash; a cell variable that stores the randomly generated feasible nonconvex instances accessible by `QC2QP_dataset.dim2.instance{k}` for k from 1 to 1000.
* `noOptimalityGapIndex` &mdash; a vector that provides the indicator of no-optimality-gap for each instance, e.g., the k-th element of this vector being 1 indicates that the instance specified by `QC2QP_dataset.dim2.instance{k}.M0`, `QC2QP_dataset.dim2.instance{k}.M1`, and `QC2QP_dataset.dim2.instance{k}.M2` has no optimality gap; and being 0 indicates the opposite. <br />

For the convenience of screening, each instance contains its own optimality gap indicator, e.g., `QC2QP_dataset.dim2.instance{k}.noOptimalityGap` being 1 means that there is no optimality gap in this instance and 0 otherwise.

---

# Author

**[Sheng Cheng](https://terpconnect.umd.edu/~cheng)**

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details

