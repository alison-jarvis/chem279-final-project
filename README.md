# CNDO/2 Molecular Property Classification

## Code Usage

### Container

This project used container based development, with the already established `berkeley-chem-179-279/dev:latest` image. To launch into this container, run `interactive.sh`. 

### Input Data

The data corresponding to the molecules of interest are contained in the `input` directory. To add a molecule, you must create a .xyz file in the `atoms` directory, which contains the optimized structure of the molecule in question. You must also create a corresponding .json file in the `input` directory which references the name of the .xyz file, and the number of $\alpha$ and $\beta$ electrons. Additionally, you must ensure that the basis functions required for the atoms within the molecule are contained within the `basis` folder. 

### Pipeline Run Instructions

Three main scripts allow you to run the pipeline. The script `build.sh` will compile the code. The script `run_pipeline.sh` will run the property classification pipeline for every molecule within the `input` folder. Finally, use the `clean.sh` script to remove previous outputs before re-building. All outputs from the pipeline will be written to a folder called `pipeline_outputs`. 

### Field Magnitude Run Instructions

Building will also create an executable called `test`, which corresponds to the test in which we varied the magnitude of the electric field and re-computed the polarizability results, to see if they were consistent across different field magnitude values. To complete this test, run `test.sh`. The results will be written out in csv file format to a folder `testing_outputs`. 

## Previous Dipole Moment Calculation

The permanent dipole is a static property of the optimized molecular structure. We can use an atom-centered approximation (for now) which essentially sums the electron cloud over each atom using the CNDO/2 SCF density matrices. For unrestricted CNDO/2 we get the that total density matrix is:

$$P_{\mu \nu} = P_{\mu \nu}^{\alpha} + P_{\mu \nu}^{\beta}$$

From the density matrix, for each atom $A$ we can calculate the net electron population on atom $A$ (which we call $N_A$) as:

$$N_A = \sum_{\mu \in A} P_{\mu \mu}$$

Then for every atom, we can calculate the net atomic charge, which is:

$$q_A = Z_A - N_A$$

Where $Z_A$ is the nuclear charge for that atom. We can then calculate the dipole moment, $\mu$, as:

$$\mathbf{\mu} = \sum_{A} q_A \mathbf{R_A}$$

Here, $\mathbf{R_A}$ is the vector center of each atom $A$. The dipole moment $\mathbf{\mu}$ will be a vector corresponding to the dipole in each cartesian direction. $\mu = |\mathbf{\mu}|$ represents the overall strength of the dipole. 


## Permanent Dipole AO method More Details

The CNDO/2 energy is approximately:

$$E_{CNDO/2} = \sum_{\mu \nu} P_{\mu \nu} (h_{\mu \nu} + f_{\mu \nu}) + E_{nuc}$$

The dipole moment ($\mu$) is the derivative of the energy with respect to a change in applied external electric field, which we call $F$. 

In equation form:

$$\mu_{x} = - \frac{\partial E}{\partial F_x}$$

So we need $\frac{\partial E}{\partial F_x}$ for CNDO/2. As a result of Hellmann-Feynman stationarity of the SCF solution (what that is not sure), this is:

$$\frac{\partial E}{\partial F_x} = \sum_{\mu \nu} P_{\mu \nu} \frac{\partial h_{\mu \nu}}{\partial F_x} + \frac{\partial E_{nuc-field}}{\partial F_x}$$

We can break this into two terms, since they are evaluated differently:

$$-\frac{\partial E}{\partial F_x} = \mathbf{\mu}_x = \mathbf{\mu}_{x}^{elec} + \mathbf{\mu}_{x}^{nuc}$$

Where the nuclear contribution is the sum of effective nuclear charges over all atomic positions:

$$\mathbf{\mu}_{x}^{nuc} = \sum_{A} Z_A \mathbf{R_{A,x}}$$

And the electronic contribution is:

$$\mathbf{\mu}_{x}^{elec} = \sum_{\mu \nu} P_{\mu \nu} \frac{\partial h_{\mu \nu}}{\partial F_x}$$

We know $P_{\mu \nu}$, which is the total density matrix from unrestricted hartree fock:

$$P_{\mu \nu} = P^{\alpha}_{\mu \nu} + P^{\beta}_{\mu \nu}$$

We need to define $\frac{\partial h_{\mu \nu}}{\partial F_x}$. Let's call this:

$$\frac{\partial h_{\mu \nu}}{\partial F_x} = D^x_{\mu \nu}$$

These terms in theory will be: (why though)

$$D^x_{\mu \nu} = \langle \chi_{\mu} | x | \chi_{\nu} \rangle$$

Where $\chi_{\mu}$ and $\chi_{\nu}$ are the atomic orbital basis functions and $x$ is position operator. 

For our basis functions, we need to work out what $D^x_{\mu \nu}$ terms will be. 

For two primitive cartesian gaussians:

$$g_A = (x - X_A)^{l_A} (y - Y_A)^{m_A} (z - Z_A)^{n_A} e^{\alpha |r - R_A|^2}$$

$$g_B = (x - X_B)^{l_B} (y - Y_B)^{m_B} (z - Z_B)^{n_B} e^{\beta |r - R_B|^2}$$

The $x$ dipole primitive integral is:

$$D_x^{AB} = \int g_A(r) x g_B(r) dr$$

Factor into $x/y/z$:

$$D_x^{AB} = D_x^{AB,x} S_y^{AB} S_z^{AB}$$

Where:

$$D_x^{AB,x} = \int (x-X_A)^{l_A} x (x - X_B)^{l_B} e^{-\alpha (x-X_A)^2 - \beta (x-X_B)^2} dx$$

Math, then:

$$D_x^{AB} = [S_x^{{AB}} (l_{A} + 1, l_{B}) + X_{A} S_x^{{AB}}(l_{A}, l_{B})] S_y^{AB} S_z^{AB}$$

Broken down into the full $D^x_{\mu \nu}$ matrix of terms:

$$D^x_{\mu \nu} = \sum_{kl} d_{k \mu} d_{l \nu} N_{k \mu} N_{l \nu} D_x^{kl}$$

Where $D_x^{kl}$ is calculated by the $D_x^{AB}$ equation. 


## Polarizability Calculation

In order to calculate the polarizability tensor, we need to modify the core hamiltonian to account for a uniform external electric field. For a full ab initio method, the change in the hamiltonian due to a uniform external electric field $\mathbf{F}$ is:

$$\hat{H}^{field} = \mathbf{F} \cdot \mathbf{r}$$

In matrix and full operator form, this would be:

$$\hat{H}_{\mu \nu}^{field} = \langle \psi_{\mu} | \mathbf{F} \cdot \mathbf{r} | \psi_{\nu} \rangle$$

Because we are using neglect of differential overlap theory though, we assume that the orbitals on different atoms don't interact significantly. Therefore, for our CNDO/2 approximation, we can simplify to:

$$\Delta H_{\mu \nu}^{field} = 0, \mu \neq \nu$$

$$\Delta H_{\mu \mu}^{field} = \mathbf{F} \cdot \mathbf{R_A}$$

So $\Delta H_{\mu \mu}^{field}$ only contains diagonal terms in the following form, where $\mu \in A$:

$$\Delta H_{\mu \mu}^{field} = F_x x_A + F_y y_A + F_z z_A$$

We have to apply this change in the hamiltonian due to the external electric field to both the fock matrix and the core hamiltonian, so that the SCF solution will produce updated density matrices and energies. The modified fock and hamiltonian matrices as a function of $\mathbf{F}$ are:

$$F_{\mu \mu}(\mathbf{F}) = F_{\mu \mu}(0) + \Delta H_{\mu \mu}^{field}$$

$$H_{\mu \mu}(\mathbf{F}) = H_{\mu \mu}(0) + \Delta H_{\mu \mu}^{field}$$

So this allows us to run SCF with an arbitrary uniform external electric field applied, which will result in modified density matrices $P_{\mu \nu}^{\alpha}$ and $P_{\mu \nu}^{\beta}$. We can use these to calculate the dipole with an applied electric field. Using this principle, we can estimate the polarizability, which is conceptually:

$$\alpha_{ij} = \left(\frac{\partial \mu_i}{\partial F_j}\right)_{\mathbf{F}=0}$$

We don't have a closed form solution for this derivative, but we can use central finite differences to approximate it, in the form:

$$\alpha_{ij} \approx \frac{\mu_i(+F_j) - \mu_i(-F_j)}{2F}$$

The indices $i$ and $j$ refer to the cartesian direction components of the field/dipole vectors. Therefore, the final result will be a $3x3$ polarizability tensor in the format:

$$\alpha = \begin{pmatrix}
\alpha_{xx} & \alpha_{xy} & \alpha_{xz} \\
\alpha_{yx} & \alpha_{yy} & \alpha_{yz} \\
\alpha_{zx} & \alpha_{zy} & \alpha_{zz} \\
\end{pmatrix}$$

We can use this tensor to calculate additional derived properties. 


## Derived Properties Calculation

The first property is the dipole magnitude, which quantifies permanent charge asymmetry. As discussed, this is:

$$|\mu| = \sqrt{\mu_x^2 + \mu_y^2 + \mu_z^2}$$

The second property is the isotropic polarizability, which quantifies how easily a molecule is polarized by an external electric field. This is calculated as:

$$\bar{\alpha} = \frac{\alpha_{xx} + \alpha_{yy} + \alpha_{zz}}{3}$$

The third property is the ansiotropic polarizability, which quantifies how directional the response to the external electric field is. This is calculated as:

$$\Delta \alpha = \left[ \frac{(\alpha_{xx} - \alpha_{yy})^2 + (\alpha_{yy} - \alpha_{zz})^2 + (\alpha_{zz} - \alpha_{xx})^2 + 6(\alpha_{xy}^2 + \alpha_{yz}^2 + \alpha_{xz}^2)}{2}  \right]^{\frac{1}{2}}$$

