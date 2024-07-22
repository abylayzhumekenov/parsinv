# parsinv

C MPI project for parallel selected inversion of sparse precision matrices and learning hyperparameters of latent Gaussian models.
For selected inversion, the code implements a hybrid approach discussed in [the paper](https://arxiv.org/abs/2309.05435).
Paper for approximating the gradient of the log-determinant is coming. In short, latent Gaussian models are defined as

$$
\begin{aligned}
    y|x,\theta  &\sim   \mathcal{N}(Ax,Q_y^{-1}) \\
    x|\theta    &\sim   \mathcal{N}(0,Q_x^{-1}) \\
    \theta      &=      \pi(\theta)
\end{aligned}
$$

The program uses gradient descent to find the mode of the Laplace approximation of posterior hyperparameters

$$
    \pi(\theta|y) \approx \frac{\pi(y|x,\theta)\pi(x|\theta)\pi(\theta)}{\pi_G(x|\theta,y)}
$$

where $\pi_G(x|\theta,y)$ is a Gaussian approximation to the full conditional. The latter has a distribution

$$
    x|\theta,y  \sim    \mathcal{N}(\mu(\theta), Q(\theta))
$$

with $Q = Q_x + A^TQ_yA$. From the optimal hyperparameter configuration $\hat{\theta}$, we can extract means and 
marginal variances of the full conditional latent field with precision $Q(\hat{\theta})$

$$
    \mu = Q^{-1}b \quad\text{and}\quad \text{diag}\Sigma = \text{diag}Q^{-1}
$$

for some right hand side $b$.


## Dependencies

* mpi
* openmp
* [pardiso](https://panua.ch/pardiso/) (requires license, not the MKL version)
* [petsc](https://petsc.org/release/) (configure with `--download-mumps --download-scalapack --download-metis --download-parmetis`)
* pcg_random (installed automatically)

After installing the dependencies, add `PETSC_DIR`, `PETSC_ARCH`, `PARDISO_DIR` to `LD_LIBRARY_PATH` and export everything.
The license path `PARDISO_LIC_PATH` and `PARDISOLICMESSAGE` have to be exported as well, the latter can be set to 1 to suppress messages.


## Download and compile

To download the source code, run

```
git clone https://github.com/abylayzhumekenov/parsinv.git
```

Compile by running `make` or `make all` in the project directory.
Clean using `make clean` and reset completely using `make distclean`.
To compile a debug version, run `make clean` and then `make DEBUG=-DDEBUG` or `make all DEBUG=-DDEBUG`.
To go back to release version, run `make clean` and `make` or `make all`.


## Run with options

Run as 
```
mpiexec [mpi_options] ./bin/main [parsinv_options]
```

MPI options:
* `-n` number of processes

Program options:
* `-ns` number of samples
* `-nn` number of neighbors
* `-ds` direct solver   (`0 = MUMPS`, `1 = PARDISO`)
* `-v` verbose mode     (`0 = FALSE`, `1 = TRUE`)
* `-p` profiling mode   (`0 = FALSE`, `1 = TRUE`)

Note that `-v 1 -p 1` works only in debug build, produced by `make DEBUG=-DDEBUG`. By default, they are ignored for efficiency.


## Input and output

Run

```
Rscript R/simulation/generate.R [rscript_options]
```

from the project directory to generate precision matrices for simulated example. Alternatively,
write your own R script for your application and save the matrices into `data` folder in the main directory.
You can check the simulation example, which uses helper functions from `R` folder to do so.
For the simulation, the arguments can be set through command line options as:

* `-nt` temporal dimension
* `-ns` spatial dimension (approximate)

The program can save the output as binary vectors, where the first 64 bits (PETSc header) can be ignored.
Note that on linux machines, R and PETSc have swapped endianness, this must be taken into account when reading and writing binary files.
