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
* [petsc](https://petsc.org/release/) (configure with `--download-mumps --download-scalapack --download-metis --download-parmetis`)

After installing the dependencies, add `PETSC_DIR`, `PETSC_ARCH`, `PARDISO_DIR` to `LD_LIBRARY_PATH` and export everything.


## Download and compile

To download the source code, run

```
git clone https://github.com/abylayzhumekenov/parsinv.git
```

Compile as a shared library by running `make`, `make all` or `make debug` in the project directory.
Clean using `make clean`. To swtich between debug and release versions, run `make clean release` or `make clean debug`.


## Compile examples

Compile examples in subfolders by changing working directory `cd examples/ex#` and running `make`.


## Generate data

Generate necessary data for examples by further changing the directory `cd R` and 
running `Rscript generate.R [rscript_options]` with the following options

* `-ns` latent spatial size
* `-nt` latent temporal size
* `-ms` data spatial size
* `-mt` data spatial size
* `-res1` spatial mesh convexity
* `-res2` spatial mesh size
* `-res3` spatial boundary size

This will generate FEM matrices as well as the data objects in the `examples/ex#/data` folder in binary format, 
where the first 64 bits (PETSc header) can be ignored. Note that on linux machines, R and PETSc have swapped endianness, 
this must be taken into account when reading and writing binary files.

You can run `fitinla.R` script with the same options to approximate hyperparameters for small to medium sized examples.


## Run with options

Run examples above as
```
mpiexec [mpi_options] ./bin/main [parsinv_options]
```

MPI options:
* `-n` number of processes

Program options:
* `-ni` number of iterations
* `-ns` number of samples
* `-no` overlap size (temporal)
* `-gd` use gradient descent (else Newton)
* `-lr` learning rate
* `-dr` decay rate
* `-ee` epsilon
* `-rt` gradient tolerance
