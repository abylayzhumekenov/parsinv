# parsinv

C MPI project for parallel selected inversion of sparse precision matrices and learning hyperparameters of latent Gaussian models.
The code implements the approach discussed in the paper (see [zhumekenov2023parallel](https://arxiv.org/abs/2309.05435)).

The models of interest are latent Gaussian models

$$
\begin{aligned}
    y|x,\theta  &\sim   \mathcal{N}(Ax,Q_y^{-1}) \\
    x|\theta    &\sim   \mathcal{N}(0,Q_x^{-1}) \\
    \theta      &\sim   \pi(\theta)
\end{aligned}
$$

with the latent field that is a solution to a non-separable space-time SPDE (see [lindgren2020diffusion](https://arxiv.org/abs/2006.04917)).
We approximate the marginal posterior of the hyperparameters using the Laplace approximation (see [rue2009approximate](https://rss.onlinelibrary.wiley.com/doi/10.1111/j.1467-9868.2008.00700.x))

$$
    \pi(\theta|y) \approx \frac{\pi(y|x,\theta)\pi(x|\theta)\pi(\theta)}{\pi_G(x|\theta,y)}
$$

where $\pi_G(x|\theta,y)$ is a Gaussian approximation to the full conditional. The latter has a distribution

$$
    x|\theta,y  \sim    \mathcal{N}(\mu(\theta), Q(\theta))
$$

with $Q = Q_x + A^TQ_yA$ and $\mu=Q^{-1}A^TQ_yy$.

To find the mode of the Laplace approximation, we use the a quasi-Newton method with an approximate diagonal Hessian. 
The gradient directions are computed from stochastic estimates of the trace $tr(Q^{-1}\partial Q)$, which involves (stochastic) selective inversion of the matrix $Q$.
Various options of the optimizer are given at the end of the document.


## Dependencies

* mpi
* cmake
* [petsc](https://petsc.org/release/) 

If the MPI is available from modules, load it as `module load mpi`. Same with Cmake.

To download PETSc, you can run the following in your preferred location
```
git clone -b release https://gitlab.com/petsc/petsc.git petsc
```

Then go the root folder `petsc` and configure as
```
./configure --download-mumps --download-scalapack --download-metis --download-parmetis
```

After configuration is complete, compile the library by substituting `/path/to/petsc` with the actual path
```
make PETSC_DIR=/path/to/petsc PETSC_ARCH=arch-linux-c-debug all
```

Check with
```
make PETSC_DIR=/path/to/petsc PETSC_ARCH=arch-linux-c-debug check
```

After installing everything, do not forget to add the `libpetsc.so` location to `LD_LIBRARY_PATH`. For example,
```
export PETSC_DIR=/path/to/petsc
export PETSC_ARCH=arch-linux-c-debug
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PETSC_DIR/$PETSC_ARCH/lib
```


## Download and compile

To download the source code of the project, run

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

Finally, run examples in `examples/ex#` folders as
```
mpiexec [mpi_options] ./bin/main [parsinv_options]
```

MPI options:
* `-n` number of processes

Program options:
* `-ni` number of optimization iterations (`default 100`)
* `-ns` number of samples (`default 10`)
* `-no` temporal overlap size, must be `> 0` (`default 1`)
* `-nr` number of iterations before restarting the learning rate (`default ni`)
* `-gd` use gradient descent if `1`, else use Newton `0` (`default 0`)
* `-lr` initial learning rate (`default 0.5`)
* `-dr` exponential decay rate (`default 1.0`)
* `-dp` polynomial decay rate, overrides `-dr` option if `> 0.5` (`default 0.0`)
* `-ee` epsilon for numerical derivative (`default 5e-3`)
* `-rt` gradient norm relative tolerance (`default 1e-3`)
* `-at` gradient norm absolute tolerance (`default 1e-3`)
* `-hh` initial hyperparameters (4 values) (`default 0.0 0.0 0.0 0.0`)

Some remarks:
* In general, `lr` should be proportional to the number of samples `ns`, but divided by a sufficiently large constant.
If `ns` is small, the variance of the gradient becomes large, and the optimizer will often overshoot. Therefore, one should consider
choosing smaller learning rate `lr` or introduce a decay with `dr` and `dp`.
* The scheduled learning rate is achieved via decay parameter `dr` (or `dp`) and the restart time `nr`.
* The optimizer automatically switches to gradient ascent/descent update if the Hessian is not negative/positive definite or too small.
* To prevent jumping on bad values randomly, update in each direction is limited to have a magnitude equal to `lr`. Therefore, restarting might be beneficial at later stages.

