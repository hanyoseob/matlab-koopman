# Koopman Theory

### Title

[Koopman theory](https://arxiv.org/abs/1607.07076)

### Abstract

We consider the application of [**Koopman theory to nonlinear partial differential equations**](https://arxiv.org/abs/1607.07076). We demonstrate that the observables chosen for constructing the Koopman operator are critical for enabling an accurate approximation to the nonlinear dynamics. If such observables can be found, then the dynamic mode decomposition algorithm can be enacted to compute a finite-dimensional approximation of the Koopman operator, including its eigenfunctions, eigenvalues and Koopman modes. Judiciously chosen observables lead to physically interpretable spatio-temporal features of the complex system under consideration and provide a connection to manifold learning methods. We demonstrate the impact of observable selection, including kernel methods, and construction of the Koopman operator on two canonical, nonlinear PDEs: Burgers' equation and the nonlinear Schrödinger equation. These examples serve to highlight the most pressing and critical challenge of Koopman theory: a principled way to select appropriate observables.

## Lectures

[Paper link](https://arxiv.org/abs/1607.07076):
Koopman theory to nonlinear partial differential equations

[Youtube link](https://youtu.be/9jZ2MsLnl0U): Koopman theory + Embeddings

## Theory

Dynamical process is formulated as follows:

$$\frac{d\vec{\mathbf{x}}}{dt} = f(\vec{\mathbf{x}}, t, \mu),$$

where $\vec{\mathbf{x}}$ defines a measurements, $t$ is a time, $\mu$ is a parametrical dependence, and $f$ indicates a system. Since the system $f$ is too complex and/or combined as well as nonlinear, it is not clear the system $f$ what is. In other words, we do not the system $f$. A lot of data $\vec{\mathbf{x}}$ is measured from the system $f$ although the system $f$ is not clear, the complex dynamical system $f$ can be approximated as follows:

$$\frac{d\vec{\mathbf{x}}}{dt} \approx A\vec{\mathbf{x}}$$

where $A$ defines a linear dynamical system which is a low-rank structure.

---

In [dynamic mode decomposition (DMD)](https://github.com/hanyoseob/matlab-DMD), **nonlinear dynamical system $f$** is performed to a least-square fitting into approximate linear dynamical system using **finite dimensionanl measurements $\vec{\mathbf{x}}$**. However, since a dynamics of the finite-dimensional measurements $\vec{\mathbf{x}}$ follows the nonlinear manifold defined by nonlinear dynamical system $f$, fitting errors exist in the approximate linear dynamical system.

To overcome the fitting errors, a koopman theory proposes a koopman operator $\mathcal{K}$ with observable function $g$. Now, we will construct new variable $\vec{\mathbf{y}}$ as follows:

$$\vec{\mathbf{y}} = g(\vec{\mathbf{x}}).$$

The observable function $g$ transforms a finite-dimensional vector space $\vec{\mathbf{x}} \in \mathbb{R}^n$ into a infinite-dimensional function space $g(\vec{\mathbf{x}}) \in \mathbb{R}^\infty$. In embedded space, the koopman operator $\mathcal{K}$ sets up a linear dynamical system on the observable function $g$:

$$
\mathcal{K} g(\vec{\mathbf{x}}_{k}) = g(\vec{\mathbf{x}}_{k+1}).
$$

On the new coordinate $\vec{\mathbf{y}}$, the above linear dynamical system can be refomulated as follows:

$$
\frac{d\vec{\mathbf{y}}}{dt} = \mathcal{K}\vec{\mathbf{y}}.
$$

Whereas the DMD was formulated as **approximate linear dynamical system** in the finite-dimensional vector space $\vec{\mathbf{x}} \in \mathbb{R}^{n}$:

$$
\frac{d\vec{\mathbf{x}}}{dt} \approx A\vec{\mathbf{x}},
$$

the koopman theory formulates exact linear dynamical system $\frac{d\vec{\mathbf{y}}}{dt} = \mathcal{K}\vec{\mathbf{y}}$ in the infinite-dimensional function space $\vec{\mathbf{y}} = g(\vec{\mathbf{x}}) \in \mathbb{R}^{\infty}$.

---

In this section, we will show how nonlinear system is organized into the linear system based on the koopman theory.

A measurement $\vec{\mathbf{x}}$ is defined as:

$$\vec{\mathbf{x}} = [\mathbf{x}_1, \mathbf{x}_2]^{\rm{T}}.$$

A nonlinear dynamical system $\frac{d\vec{\mathbf{x}}}{dt} = f(\vec{\mathbf{x}}, t)$ is formulated as:

$$
\frac{d\mathbf{x}_1}{dt} = \dot{\mathbf{x}}_1 = \mu \mathbf{x}_1,
$$

$$
\frac{d\mathbf{x}_2}{dt} = \dot{\mathbf{x}}_2 = \lambda(\mathbf{x}_2 - \mathbf{x}_1^2).
$$

For the nonlinear dynamical system $f$ with measurement $\vec{\mathbf{x}}$, new variable $\vec\mathbf{y}$ will be constructed by using observable function $g$ as follows:

$$
g(\vec\mathbf{x})=\vec\mathbf{y} = [\mathbf{y}_1, \mathbf{y}_2, \mathbf{y}_3]^{\rm{T}}.
$$

Each variables is defined as:

$$
\mathbf{y}_1 = \mathbf{x}_1,
$$

$$
\mathbf{y}_2= \mathbf{x}_2,
$$

$$
\mathbf{y}_3 = \mathbf{x}_1^2.
$$

In the koopman embedding coordinate $\vec\mathbf{y}$, its dynamical system is calculated as follows:

$$
\dot{\mathbf{y}_1} = \dot{\mathbf{x}_1} = \mu \mathbf{x}_1 = \mu \mathbf{y}_1,
$$

$$
\dot{\mathbf{y}_2} = \dot{\mathbf{x}_2} = \lambda(\mathbf{x}_2 - \mathbf{x}_1^2) = \lambda(\mathbf{y}_2 - \mathbf{y}_3),
$$

$$
\dot{\mathbf{y}_3} = 2\mathbf{x}_1\dot{\mathbf{x}_1} = 2\mathbf{x}_1\mu\mathbf{x}_1 = 2\mu\mathbf{x}_1^2 = 2\mu\mathbf{y}_3.
$$

Above equations are reformulated as matrix form:

$$
\frac{d\vec{\mathbf{y}}}{dt} =
\begin{bmatrix}
 \mu &    0    & 0 \\
  0  & \lambda & -\lambda \\
  0  &    0    & 2\mu
\end{bmatrix}
\vec{\mathbf{y}} = \mathcal{K} \vec\mathbf{y}.
$$

Thanks to the koopman embedding ($\vec\mathbf{x} \rarr \vec\mathbf{y}$), the original nonlinear dynamical system $f$ was converted to the linear dynamical system $\mathcal{K}$.

---

Previous example was simply linearized using single coordinate embedding. However, it is difficult to linearize a typical nonlinear systems using simple coordinate embedding, and it requires high-dimensional embedding, sometimes infinite dimensional embedding.

This example changed quadratic term from $\mathbf{x}_1$ to $\mathbf{x}_2$.

$$
\frac{d\mathbf{x}_1}{dt} = \dot{\mathbf{x}}_1 = \mu \mathbf{x}_1,
$$

$$
\frac{d\mathbf{x}_2}{dt} = \dot{\mathbf{x}}_2 = \lambda(\mathbf{x}_2^2 + \mathbf{x}_1).
$$

Similar to previous one, the nonlinear dynamical system was performed with a coordinate embedding form $\vec\mathbf{x}$ to $\vec\mathbf{y}$.

Each variables $\vec\mathbf{y}$ is defined as:

$$
\mathbf{y}_1 = \mathbf{x}_1,
$$

$$
\mathbf{y}_2 = \mathbf{x}_2,
$$

$$
\mathbf{y}_3 = \mathbf{x}_2^2.
$$

In the koopman embedding coordinate $\vec\mathbf{y}$, its dynamical system is calculated as follows:

$$
\dot{\mathbf{y}_1} = \dot{\mathbf{x}_1} = \mu \mathbf{x}_1 = \mu \mathbf{y}_1,
$$

$$
\dot{\mathbf{y}_2} = \dot{\mathbf{x}_2} = \lambda(\mathbf{x}_2^2 - \mathbf{x}_1) = \lambda(\mathbf{y}_3 - \mathbf{y}_1),
$$

$$
\dot{\mathbf{y}_3} = 2\mathbf{x}_2\dot{\mathbf{x}_2} = 2\mathbf{x}_2\lambda(\mathbf{x}_2^2 - \mathbf{x}_1) = 2\lambda(\mathbf{x}_2^3 - \mathbf{x}_1\mathbf{x}_2).
$$

For this example, the dynamical system was not linearized using single coordinate embedding $\mathbf{y}_3$, and new nonliear variables are generated such as $\mathbf{x}_2^3$ and $\mathbf{x}_1\mathbf{x}_2$. To linearize according to new variables $\mathbf{x}_2^3$ and $\mathbf{x}_1\mathbf{x}_2$, additional coordinate embeddings are performed as follows:

$$
\mathbf{y}_4 = \mathbf{x}_2^3,
$$

$$
\mathbf{y}_5 = \mathbf{x}_1\mathbf{x}_2.
$$

In addition, these cascade coordinate embeddings are performed consecutively until all nonlinear variables are linearized.

Specifically, the koopman embedding is conducted for all multi fixed points. In other hands, if the dynamical system takes a single fixed point, the system is linear and does not require coordinate embedding.

---

So far we have briefly reviewed how to linearize from a finite-dimensional nonlinear dynamical system to an infinite-dimensional linear dynamical system based on koopman theory.

However, we does not cover how to find an optimal observation function $g$. Actually, it is difficult to determine the proper observation function $g$ and the observation function $g$ is sensitively related to a performance of a koopman operator. Therefore, to determine the observation function $g$, we need expertise in the dynamical system.

---

An algorithm is closly similar to [DMD](https://github.com/hanyoseob/matlab-DMD).

First, measured vector space $\vec\mathbf{x}$ is transformed to function space $\vec\mathbf{y}$ using observation function $g$.

Next, Using the observables $\vec\mathbf{y}$, a linear dynamical system $A_Y$ is constructed as follows:

$$ 
\bar{Y} = 
\begin{bmatrix}
  &  &  &  \\
  &  &  &  \\
\rm{y}_1 & \rm{y}_2 & \cdots & \rm{y}_{m-1}\\
  &  &  &  \\
  &  &  &  
\end{bmatrix}.
$$

Another matrix shifted by 1 time step is defined as:

$$ 
\bar{Y}' = 
\begin{bmatrix}
  &  &  &  \\
  &  &  &  \\
\rm{y}_2 & \rm{y}_3 & \cdots & \rm{y}_{m}\\
  &  &  &  \\
  &  &  &  
\end{bmatrix}.
$$

Therefore, the linear dynamical system $A_Y$ is satisfied with the relationship below:

$$\bar{Y}' = A_Y \bar{Y},$$

where $\bar{Y}'$ and $\bar{Y}$ are the future state of $\bar{Y}$ and the current state, respectively.

The linear dynamical system $A_Y$ can be extracted using a pseudo inverse $\bar{Y}^{\dagger}$ of $\bar{Y}$:

$$A_Y = \bar{Y}' \bar{Y}^{\dagger}.$$

After configuring the linear dynamical system $A_Y$, the DMD is applied to extract the koopman eigenvalues $\Phi$ and modes $\Lambda$ (or $\Omega = \log{\Lambda}$).

Using the eigen vectors $\Phi$ and the eigen values $\Lambda$, the solution of the observables $\vec\mathbf{y}^*$ can be calculated as:

$$
\rm{y}_t^* = \Phi e ^{\Omega t} \rm{b} = \sum_{k=1}^{r} \phi_k e^{\omega_k t}b_k.
$$

Finally, reconstructed measurement $\vec\mathbf{x}^*$ can be transformed from the observables using inverse observation function $g^{-1}$:

$$
\vec\mathbf{y}^* = g(\vec\mathbf{x}^*) \rarr \vec\mathbf{x}^* = g^{-1}(\vec\mathbf{y})^*.
$$
