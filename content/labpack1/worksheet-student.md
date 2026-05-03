# Lab Pack #1: From Boys Functions to SCF Convergence

**Lab Pack:** 1 - From Boys to Orbitals
**Version:** 1.0
**Last Updated:** 2026-04-05
**Estimated Time:** 60-90 minutes

---

## Introduction

Welcome to Lab Pack #1 of the Interactive Quantum Chemistry Playground (IQCP). In this guided exploration, you will build intuition for three fundamental concepts that underpin all molecular electronic structure calculations:

1. **Boys functions** - special mathematical functions that appear in molecular integral evaluation
2. **Rys quadrature** - a numerical integration technique optimized for quantum chemistry
3. **Self-Consistent Field (SCF) convergence** - the iterative process that finds molecular orbitals

By the end of this lab, you will understand not just *what* these quantities are, but *why* they matter for practical quantum chemistry calculations.

---

## Learning Objectives

By completing this lab, you will be able to:

1. **LO1:** Describe the qualitative behavior of Boys functions $F_m(T)$ as $T$ varies, and explain why different computational regimes are necessary for numerical stability.

2. **LO2:** Explain the relationship between Rys quadrature order and integration accuracy, and select an appropriate order to meet a specified error target.

3. **LO3:** Interpret SCF energy convergence plots and explain how DIIS (Direct Inversion in the Iterative Subspace) accelerates convergence.

4. **LO4:** Connect numerical parameters (quadrature order, convergence thresholds, DIIS settings) to computational outcomes and cost.

5. **LO5:** Export run artifacts that document your computational explorations for reproducibility.

6. **LO6:** Apply physical and mathematical reasoning to predict computational behavior before running calculations.

---

## What You Will Need

- A modern web browser (Chrome, Firefox, Safari, or Edge)
- This worksheet (print or digital)
- Access to the IQCP web application at: https://iqcp.dev
- Approximately 60-90 minutes of uninterrupted time

---

## What to Submit

At the end of this lab, you should have:

- [ ] Written answers to all numbered questions (24 total)
- [ ] Three exported run artifacts (one per section)
- [ ] Your responses to the synthesis questions

Your instructor will provide specific submission instructions.

---

## How to Use Deep Links

Throughout this worksheet, you will see links that look like:

> **[Click here to open this configuration](https://iqcp.dev/v1/boys?run=...)**

These "deep links" will open IQCP with specific parameter settings already configured. This ensures you see exactly what the worksheet describes. If a link does not work, you can manually set the parameters using the controls in the application.

---

## Section 1: Introduction and Warm-Up (~5 min)

Before diving into the interactive exploration, take a moment to recall some foundational concepts.

The **Boys function** $F_m(T)$ is defined as:

$$F_m(T) = \int_0^1 t^{2m} e^{-Tt^2} \, dt$$

This integral appears naturally when evaluating molecular integrals over Gaussian-type orbitals (GTOs). The parameter $m$ is a non-negative integer related to angular momentum, and $T$ is a positive real number related to distances between atomic centers.

**Quick recall questions (not graded, just to warm up your thinking):**

- What is the value of the integral when T = 0?
- What happens to the integrand as T becomes very large?
- Why might direct numerical integration be challenging for this function?

When you are ready, proceed to Section 2.

---

## Section 2: Boys Function Exploration (~15-20 min)

In this section, you will explore how the Boys function $F_m(T)$ behaves as you vary its parameters, and discover why quantum chemistry codes use different computational methods for different parameter ranges.

### Step 2.1: Getting Started

Open IQCP to the Boys Function Lab (Module C) with the default settings.

**[Open: Default Boys View (m=0, T=1.0)](https://iqcp.dev/v1/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgNYAVIg1YosMAO79seAOYKQAX1aksOyTMXEYAD3ZyCFgETk5AA)**

Familiarize yourself with the interface:
- The **Controls Panel** on the left allows you to change m and T
- The **Main Display** shows the computed value of $F_m(T)$
- Use the **Mode Toggle** to switch between "Explain" and "Internals" views

### Step 2.2: Small T Behavior (Series Regime)

Now explore what happens when T is small.

**[Open: Small T (m=0, T=0.5, Internals Mode)](https://iqcp.dev/v1/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgNYAVfUICsrFFhgB3ftjwBzBSAC+rUlh2SZi4ljx06HgIchouLkA)**

Set m = 0 and T = 0.5, then switch to "Internals" mode.

Observe:
- What computational regime is being used?
- What is the computed value of $F_0(0.5)$?

Now try an even smaller T value:

**[Open: Very Small T (m=0, T=0.01)](https://iqcp.dev/v1/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgNYAVfWIasUWGAHd+2PAHMFIAL6tSWHZJmLiMAB7scghYBM7OQA)**

Set T = 0.01 and observe the result.

**Q2.1:** What value does $F_0(T)$ approach as $T$ approaches 0? Express this as a simple fraction in terms of $m$.

*Your answer:* _______________________________________________

### Step 2.3: Large T Behavior (Recurrence Method)

Now explore what happens when T becomes larger.

**[Open: Moderate T (m=0, T=15.0, Internals Mode)](https://iqcp.dev/v1/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgNYAVIgwCsrFFhgB3ftjwBzBSAC+rUlh2SZi4ljx06HgIchouLkA)**

Set T = 15.0 and observe:
- Which computational method is now being used?
- How does this differ from what you might see for small T with higher m values?

**[Open: Large T (m=0, T=35.0, Internals Mode)](https://iqcp.dev/v1/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgNYAVIgGYArKxRYYAd37Y8AcwUgAvq1JYdkmYuJY8OnQ8BDkNV1cgA)**

Now increase T to 35.0 and observe.

**Q2.2:** As $T$ increases, $F_0(T)$ approaches zero. At approximately what $T$ value does $F_0(T)$ become smaller than 0.01? (You can use the slider to find this.)

*Your answer:* _______________________________________________

**Q2.3:** In the Internals panel, what computational method is used for T = 0.5? For T = 15.0? For T = 35.0? (Note: The method may depend on BOTH T and m. For m=0, observe what method IQCP uses.)

*Your answer:*
- T = 0.5 (m=0): _______________________________________________
- T = 15.0 (m=0): _______________________________________________
- T = 35.0 (m=0): _______________________________________________

**Bonus exploration:** Try changing $m$ to 5 or 10. Does the method change for small $T$ values? What does this tell you about how method selection works?

### Step 2.4: Effect of Order m

The Boys function has an order parameter $m$. Let us see how this affects the function.

**[Open: Higher Order m (m=5, T=10.0)](https://iqcp.dev/v1/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgFZWAFSIMBrFFhgB3ftjwBzBSAC+rUlh2SZi4ljx06HgIchouLkA)**

Set m = 5 and T = 10.0.

**Q2.4:** Compare $F_0(10.0)$ and $F_5(10.0)$. Which is larger? Why do you think higher-order Boys functions decay more rapidly with increasing $T$?

*Your answer:* _______________________________________________

_______________________________________________

### Step 2.5: Regime Boundaries

**[Open: Regime Boundary (m=0, T=10.0, Internals Mode)](https://iqcp.dev/v1/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgNYAVIgwMgUWGAHd+2PAHMFIAL6tSWHZJmLiMAB7scghYBM7OQA)**

Set T = 12.0 (near the series/recurrence boundary).

**Q2.5:** Why do you think quantum chemistry codes use different computational methods (series expansion vs. recurrence relation) for different parameter values, rather than using a single method for all cases? Consider both $T$ and $m$ in your answer.

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

### Step 2.6: Theoretical Regimes vs. Implementation Methods

The lecture notes describe **three theoretical regimes** for computing Boys functions:

| Regime | Condition | Theoretical Method |
|--------|-----------|-------------------|
| Small $T$ | $T < 25$ | Series expansion |
| Moderate $T$ | $25 \le T < 30+5m$ | erf formula + upward recurrence |
| Large $T$ | $T \ge 30+5m$ | Asymptotic expansion |

However, IQCP (following libcint) uses a **two-method implementation** with $m$-dependent turnover points.

**[Open: Compare m=0 at T=35 (Internals Mode)](https://iqcp.dev/v1/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgNYAVIgGYArKxRYYAd37Y8AcwUgAvq1JYdkmYuJY8OnQ8BDkNV1cgA)**

Set $m=0$ and $T=35.0$ in Internals mode. According to the theoretical table, this would be in the "Large $T$ / Asymptotic" regime (since $35 \ge 30+5 \times 0 = 30$).

**[Open: Compare m=5 at T=45 (Internals Mode)](https://iqcp.dev/v1/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgFZWAFSIAWQyBRYYAd37Y8AcwUgAvq1JYdkmYuJY8OnQ8BDkNV1cgA)**

Now set $m=5$ and $T=45.0$. According to the theoretical table, this would be in the "Moderate $T$" regime (since $30+5 \times 5 = 55$, and $45 < 55$).

**Q2.6:** Looking at the Internals panel, what computational method does IQCP actually use for $m=0$ at $T=35$? For $m=5$ at $T=45$? According to the theoretical 3-regime model, what theoretical regime would these fall into?

*Your answer:*
- m=0, T=35: IQCP method: _______________ Theoretical regime: _______________
- m=5, T=45: IQCP method: _______________ Theoretical regime: _______________

**Q2.7:** The theoretical description in the lecture notes describes 3 regimes (series, moderate/erf-recurrence, asymptotic), but IQCP implements only 2 methods. Why might implementations combine the moderate and large $T$ methods rather than implementing all three separately?

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

---

### Checkpoint: Boys Function Artifact

Before moving on, export a run artifact to document your exploration.

**[Open: Checkpoint State (m=3, T=9.5)](https://iqcp.dev/v1/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgGZWAFSIBOIQFZWKLDADu-bHgDmCkAF9WpLDskzFxFh4dOh4CHIa7u5AA)**

1. Set m = 3 and T = 9.5
2. Click the **Export** button
3. Save the artifact file as `boys-artifact.json`

- [ ] I have exported my Boys function artifact

---

## Section 3: Rys Quadrature Exploration (~15-20 min)

Rys quadrature is a specialized numerical integration technique used to evaluate molecular integrals. In this section, you will explore how the quadrature order affects accuracy and learn to choose appropriate settings.

### Step 3.1: Introduction to Rys Quadrature

**[Open: Default Rys View (n=3, T=10.0)](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEAzKwAqRBgNY10AcxhV+DGAFoAbCAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

Open the Rys Quadrature Lab (Module D) and observe the interface:
- The **Controls Panel** lets you adjust the quadrature order $n$ and the parameter $T$
- The **Roots/Weights Table** shows the quadrature points and their weights
- The **Error Curve** (when visible) shows reconstruction error vs. order

Rys quadrature computes $n$ roots $(t_i)$ and $n$ weights $(w_i)$ such that:

$$\int_0^1 f(t^2) \, e^{-Tt^2} \, dt \approx \sum_{i=1}^n w_i \, f(t_i^2)$$

### Step 3.2: Inspecting Roots and Weights

**[Open: Roots and Weights Inspection (n=5, T=10.0, Internals Mode)](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEArKwAqRBgNY10AcxhV+DGAFoAbCAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

Set n = 5 and T = 10.0, then switch to Internals mode.

**Q3.1:** Examine the roots $(t_i)$ and weights $(w_i)$ displayed in the table. Are all roots strictly between 0 and 1? Are all weights positive?

*Your answer:*
- Roots in (0, 1)? _______________________________________________
- Weights positive? _______________________________________________

**[Open: High Quadrature Order (n=7, T=10.0, Internals Mode)](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEA7KwAqRBgNY10AcxhV+DGAFoAbCAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

Increase the order to n = 10 and compare.

**Q3.2:** How does the number of quadrature points affect the computation? Think about both accuracy and computational cost.

*Your answer:* _______________________________________________

_______________________________________________

### Step 3.3: Error vs. Quadrature Order

The reconstruction error measures how accurately the quadrature can reproduce the exact integral moments.

**[Open: Error Curve at T=10.0 (n=5, T=10.0, Explain Mode)](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEArKwAqRBgNY10AcxhV+DGAFoAbCAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

In Explain mode, observe the error curve or error information displayed.

**Q3.3:** Looking at the error information, what is the approximate maximum reconstruction error for n=3? For n=5? For n=7?

*Your answer:*
- n=3: _______________________________________________
- n=5: _______________________________________________
- n=7: _______________________________________________

**[Open: Error Curve at T=25.0 (n=3, T=25.0)](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEAzKwAqRAEwBWVjXQBzGFX4MYAWgBsIAL6tSWIqCnTFxDAAHuxyCFgEbm5AA)**

Change T to 25.0 and observe how the error behavior changes.

### Step 3.4: Choosing Optimal Quadrature Order

**[Open: Target 1e-6 Accuracy (n=5, T=10.0)](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEArKwAqRBgNY10AcxhV+DGAFoAbCAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

With T = 10.0, use the target accuracy selector to set the target to 1e-6.

**[Open: Target 1e-8 Accuracy (n=5, T=15.0)](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEArKwAqRBvpA10AcxhV+DGAFoAHCAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

Now change the target to 1e-8.

**Q3.4:** At T = 10.0, what is the minimum quadrature order needed to achieve 1e-8 accuracy? What about 1e-6 accuracy?

*Your answer:*
- For 1e-8: _______________________________________________
- For 1e-6: _______________________________________________

**Q3.5:** How does the recommended quadrature order change when T increases from 10 to 25? Why do you think this happens?

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

### Step 3.5: Shell Quartet and Root Count Rule

The number of Rys quadrature roots needed depends on the **total angular momentum** $L$ of the shell quartet. The lecture notes provide the formula:

$$n_r = \lfloor L/2 \rfloor + 1$$

where $L = l_A + l_B + l_C + l_D$ is the sum of angular momenta for the four shells in an electron repulsion integral $(ab|cd)$.

**[Open: Shell Quartet Selector (Internals Mode)](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEAzKwAqRBgNY10AcxhV+DGAFoAbCAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

Using the shell quartet selector in the Internals panel, try different shell combinations.

**Q3.6:** Using the shell quartet selector, set the quartet to (pp|pp). What is the total angular momentum $L$? What quadrature order $n$ is automatically selected by IQCP? Verify this matches the formula $n_r = \lfloor L/2 \rfloor + 1$.

*Your answer:*
- Total angular momentum L = _______________
- IQCP selected order n = _______________
- Formula check: floor(___/2) + 1 = _______________
- Do they match? _______________

**[Open: (dd|pp) Shell Quartet](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEALKwAqRBgNY10AcxhV+DGAFoAbCAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

Now try (dd|pp).

**Q3.7:** For the (dd|pp) shell quartet: What is $L$? What quadrature order does IQCP select? According to the root count rule, is this the minimum required order?

*Your answer:*
- L = _______________
- IQCP selected order = _______________
- Minimum required by formula = _______________

### Step 3.6: Algorithm 5.1 Pipeline

The lecture notes describe **Algorithm 5.1** for computing Rys nodes and weights from moments:

1. Compute moments: $\mu_k(T) = 2F_k(T)$ for $k = 0, 1, \ldots, 2n-1$
2. Form Hankel matrix $H_{ij} = \mu_{i+j}(T)$
3. Form shifted Hankel matrix $H^{(1)}_{ij} = \mu_{i+j+1}(T)$
4. Cholesky factorize: $H = LL^T$, then $C = L^{-1}$
5. Build Jacobi matrix: $J = C H^{(1)} C^T$
6. Eigendecomposition: nodes = eigenvalues, weights = $\mu_0 \cdot (V_{0i})^2$

**[Open: Algorithm 5.1 Internals (n=3, T=10)](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEAzKwAqRBgNY10AcxhV+DGAFoAbCAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

Set n=3 and T=10, then examine the Algorithm 5.1 pipeline in the Internals panel.

**Q3.8:** Looking at the Algorithm 5.1 internals for $T=10$, $n=3$:
- What are the first 3 moments $\mu_0, \mu_1, \mu_2$? (These are $2F_k(T)$ values)
- What is the dimension of the Hankel matrix $H$?

*Your answer:*
- $\mu_0 = 2F_0(10) =$ _______________
- $\mu_1 = 2F_1(10) =$ _______________
- $\mu_2 = 2F_2(10) =$ _______________
- Hankel matrix dimension: ___ x ___

---

### Checkpoint: Rys Quadrature Artifact

Export a run artifact for the Rys Quadrature Lab (Module D).

**[Open: Checkpoint State (n=5, T=15.0, target=1e-8)](https://iqcp.dev/v1/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEArKwAqRBvpA10AcxhV+DGAFoAHCAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

1. Set T = 15.0 with target accuracy 1e-8
2. Note the recommended order
3. Click **Export** and save as `rys-artifact.json`

- [ ] I have exported my Rys quadrature artifact

---

## Section 4: SCF Convergence Exploration (~20-30 min)

The Self-Consistent Field (SCF) method is the workhorse of computational quantum chemistry. In this section, you will observe how SCF iteratively finds molecular orbitals and how DIIS acceleration dramatically improves convergence.

### Step 4.1: Understanding the SCF Process

**[Open: H2 Default Run (medium convergence, DIIS enabled)](https://iqcp.dev/v1/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn8oAJkpUcAZgDmSdKIAsExAA99ddEQCsA1hHwp+OBhpLFI4CRCsMCIqMhgAX1ZSLA1JGUViGGd2OQQsAnj4oA)**

Open the SCF Sandbox (Module E) with the default H2 molecule.

The SCF procedure works as follows:
1. Start with an initial guess for the density matrix
2. Build the Fock matrix from the density
3. Diagonalize to get new orbitals
4. Form a new density from occupied orbitals
5. Check for convergence; if not converged, return to step 2

Observe the iteration table and convergence plot.

**Q4.1:** For H2 with medium convergence and DIIS enabled, how many iterations does it take to converge? What is the final energy?

*Your answer:*
- Iterations: _______________________________________________
- Final energy: _______________________________________________ Hartree

### Step 4.2: Effect of DIIS on Convergence

Now let us compare convergence with and without DIIS.

**[Open: H2 Tight Convergence WITHOUT DIIS](https://iqcp.dev/v1/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn8oAJkpUcAZgDmSdKIAsExAA99ddEQCsA1hHwp+KiwrKCoJaSwsMCJlBDkwGABfVlIsDUkZRWIYZ3Y5BCwCRMSgA)**

Set convergence to "tight" and DIIS to OFF. Run the calculation and observe:
- How many iterations?
- What is the convergence pattern?

**[Open: H2 Tight Convergence WITH DIIS](https://iqcp.dev/v1/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn8oAJkpUcAZgDmSdKIAsExAA99ddEQCsA1hHwp+KiwrKCoJaSwsMCIqMhgAX1ZSLA1JGUViGGd2OQQsAnj4oA)**

Now enable DIIS and run again.

**Q4.2:** For H2 with tight convergence, how many iterations does SCF take without DIIS? With DIIS?

*Your answer:*
- Without DIIS: _______________________________________________
- With DIIS: _______________________________________________

**Q4.3:** Looking at the energy vs. iteration plot, describe the difference in convergence patterns between the two cases. How does DIIS change the shape of the convergence curve?

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

### Step 4.3: Larger Molecule - Water

**[Open: H2O WITHOUT DIIS](https://iqcp.dev/v1/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn8oAJhyUqOAMwBzCYgAe+uuiIBWAawj4U-ODGlYpHASAVhgRMoIcmAwAL6spFgakjKKxDD27HIIWASxsUA)**

Now switch to H2O (water molecule) without DIIS.

**[Open: H2O WITH DIIS](https://iqcp.dev/v1/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn8oAJhyUqOAMwBzCYgAe+uuiIBWAawj4U-ODGlYpHASAVhgRFRkMAC+rKRYGpIyisQw9uxyCFgE0dFAA)**

Now enable DIIS for H2O.

**Q4.4:** For H2O, what is the final RHF energy? Does the calculation converge in both cases (with and without DIIS)?

*Your answer:*
- Final energy: _______________________________________________ Hartree
- Converges without DIIS? _______________________________________________
- Converges with DIIS? _______________________________________________

### Step 4.4: Inspecting SCF Internals

**[Open: H2 Matrix Inspection (Internals Mode)](https://iqcp.dev/v1/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn8oAJkpUcAZgDmSdKIAsExAA99ddEQCsA1hHwp+OBhpLFI4CRCsMCIqMhgAX1ZSLA1JGUViLDx3PAQ5aPj4oA)**

Switch to Internals mode and examine the matrices.

**Q4.5:** In the Internals mode, examine the Fock matrix $F$. Is it symmetric ($F_{ij} = F_{ji}$)? Why is symmetry of the Fock matrix physically important?

*Your answer:* _______________________________________________

_______________________________________________

### Step 4.5: Orbital Energies

Still in the H2 calculation, look at the orbital energies (eigenvalues of the Fock matrix).

**Q4.6:** What is the HOMO energy for H2? What is the LUMO energy? What does the HOMO-LUMO gap tell you about the molecule?

*Your answer:*
- HOMO energy: _______________________________________________
- LUMO energy: _______________________________________________
- HOMO-LUMO gap significance: _______________________________________________

### Step 4.6: A More Complex System

**[Open: LiH System (Internals Mode)](https://iqcp.dev/v1/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn85WKJSo4AzAHMJiAB7666IgFYBrCPhT84MaVlI4CX8sMCIqMhgAX1ZSLA1JGUViGDt2OQQsAiiooA)**

Try the LiH (lithium hydride) system.

**Q4.7:** Based on your observations in this section, when would you recommend using DIIS? Are there any situations where DIIS might not help or could cause problems?

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

> **Going Deeper:** The SCF Sandbox also supports DFT methods (LDA, B3LYP), multiple basis sets (STO-3G through cc-pVDZ), geometry optimization, potential energy surface scanning, and 3D orbital visualization. These features are explored in **Lab Pack #2: 3D Exploration, PES, and Orbitals**.

---

### Checkpoint: SCF Artifact

Export your final SCF run artifact.

**[Open: Checkpoint State (H2O with DIIS)](https://iqcp.dev/v1/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn8oAJhyUqOAMwBzCYgAe+uuiIBWAawj4U-ODGlYpHASAVhgRFRkMAC+rKRYGpIyisQw9uxyCFgE0dFAA)**

1. Run H2O with medium convergence and DIIS enabled
2. Switch to Internals mode to ensure matrix data is captured
3. Click **Export** and save as `scf-artifact.json`

- [ ] I have exported my SCF artifact

---

## Section 5: Synthesis and Reflection (~5-10 min)

Now that you have explored all three modules, let us connect the concepts.

### Connecting the Concepts

The Boys function, Rys quadrature, and SCF are deeply connected in quantum chemistry:

1. **Boys functions** appear when evaluating molecular integrals over Gaussian basis functions
2. **Rys quadrature** uses Boys function values as moments to compute efficient quadrature rules
3. **SCF calculations** require evaluating many molecular integrals at each iteration

Understanding each piece helps you appreciate the computational challenges of electronic structure theory.

---

**Q5.1:** The Boys function $F_m(T)$ appears in nuclear attraction integrals where $T$ depends on the distance between nuclei and basis function centers. Based on what you learned about Boys function behavior:

a) What happens to the integrand when two Gaussian centers are very close (small $T$)?
b) What happens when they are far apart (large $T$)?
c) Why might special numerical care be needed at the regime boundaries?

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

---

**Q5.2:** Rys quadrature is used to evaluate two-electron integrals $(ij|kl)$, where each integral can require many quadrature points for high accuracy. If a calculation requires $10^{-10}$ accuracy and you need to evaluate 10,000 integrals:

a) How would you estimate the total number of quadrature point evaluations needed?
b) How does your understanding of the order-accuracy relationship help predict computational cost?
c) Why might adaptive quadrature order selection (varying order based on the integral's $T$ value) be valuable?

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

---

**Q5.3:** DIIS (Direct Inversion in the Iterative Subspace) dramatically accelerates SCF convergence. Based on your observations:

a) Summarize in 2-3 sentences what DIIS does to improve convergence.
b) The DIIS method works by extrapolating from previous Fock matrices. Why might this be more effective than simple iteration?
c) Under what circumstances might standard SCF iteration (without DIIS) still be useful?

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

---

## Final Deliverables Checklist

Before submitting, verify that you have completed:

**Section 2: Boys Function (7 questions)**
- [ ] **Q2.1:** Boys function limiting value as T approaches 0
- [ ] **Q2.2:** T value where F_0(T) < 0.01
- [ ] **Q2.3:** Computational regimes at different T values
- [ ] **Q2.4:** Comparison of F_0 and F_5
- [ ] **Q2.5:** Why different computational methods are used
- [ ] **Q2.6:** IQCP methods vs. theoretical 3-regime model
- [ ] **Q2.7:** Why implementations combine moderate/large T methods

**Section 3: Rys Quadrature (8 questions)**
- [ ] **Q3.1:** Roots and weights properties
- [ ] **Q3.2:** Effect of quadrature points on accuracy/cost
- [ ] **Q3.3:** Reconstruction errors at different orders
- [ ] **Q3.4:** Minimum orders for target accuracies
- [ ] **Q3.5:** Effect of T on recommended order
- [ ] **Q3.6:** Shell quartet (pp|pp) root count verification
- [ ] **Q3.7:** Shell quartet (dd|pp) root count verification
- [ ] **Q3.8:** Algorithm 5.1 moments and Hankel matrix

**Section 4: SCF (7 questions)**
- [ ] **Q4.1:** H2 default iterations and energy
- [ ] **Q4.2:** H2 iterations with/without DIIS
- [ ] **Q4.3:** Convergence pattern description
- [ ] **Q4.4:** H2O final energy and convergence
- [ ] **Q4.5:** Fock matrix symmetry
- [ ] **Q4.6:** HOMO-LUMO energies and gap
- [ ] **Q4.7:** DIIS recommendations

**Section 5: Synthesis (3 questions)**
- [ ] **Q5.1:** Boys function in integral evaluation
- [ ] **Q5.2:** Quadrature and computational cost
- [ ] **Q5.3:** DIIS summary and explanation

**Artifacts:**
- [ ] `boys-artifact.json` (m=3, T=9.5)
- [ ] `rys-artifact.json` (n=5, T=15.0, target=1e-8)
- [ ] `scf-artifact.json` (H2O with DIIS)

---

## Appendix A: Troubleshooting

**Problem:** Deep links do not load the expected state
**Solution:** Try clearing your browser cache, or manually enter the parameters shown in the link description.

**Problem:** The calculation seems stuck or unresponsive
**Solution:** Wait a few seconds; larger systems may take time. If it persists, refresh the page and try again.

**Problem:** Export button does not appear or does not work
**Solution:** Ensure the calculation has completed before exporting. Check that your browser allows file downloads from the site.

**Problem:** Numbers look different from what my classmate sees
**Solution:** Small numerical differences (beyond the 10th decimal place) are normal due to floating-point arithmetic. Focus on agreement within reasonable tolerances.

---

## Appendix B: Glossary

**Boys function** ($F_m(T)$): A special function defined as an integral over the unit interval, appearing in molecular integral evaluation over Gaussian basis functions.

**DIIS (Direct Inversion in the Iterative Subspace):** An acceleration technique that extrapolates from previous iterations to improve SCF convergence.

**Fock matrix:** The effective one-electron Hamiltonian in Hartree-Fock theory, which depends on the electron density and must be solved self-consistently.

**HOMO (Highest Occupied Molecular Orbital):** The highest-energy orbital that contains electrons in the ground state.

**LUMO (Lowest Unoccupied Molecular Orbital):** The lowest-energy orbital that does not contain electrons in the ground state.

**Quadrature:** A numerical method for approximating definite integrals using weighted sums of function values at specific points.

**RHF (Restricted Hartree-Fock):** A variant of Hartree-Fock theory where alpha and beta electrons share the same spatial orbitals, used for closed-shell molecules.

**Rys quadrature:** A Gaussian quadrature scheme optimized for integrals appearing in quantum chemistry, where the weight function includes an exponential factor.

**SCF (Self-Consistent Field):** An iterative method where the electron density is updated until it becomes consistent with the orbitals derived from it.

---

## Appendix C: References

1. Shavitt, I. (1963). "The Gaussian Function in Calculations of Statistical Mechanics and Quantum Mechanics." *Methods in Computational Physics*, Vol. 2, pp. 1-45.

2. Dupuis, M., Rys, J., and King, H.F. (1976). "Evaluation of molecular integrals over Gaussian basis functions." *Journal of Chemical Physics*, 65(1), 111-116.

3. Pulay, P. (1980). "Convergence acceleration of iterative sequences. The case of SCF iteration." *Chemical Physics Letters*, 73(2), 393-398.

4. Pulay, P. (1982). "Improved SCF convergence acceleration." *Journal of Computational Chemistry*, 3(4), 556-560.

5. Szabo, A. and Ostlund, N.S. (1996). *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*. Dover Publications.

---

**Lab Pack #1 Complete**

Thank you for completing this guided exploration of quantum chemistry fundamentals. The concepts you explored today - Boys functions, Rys quadrature, and SCF convergence - form the computational foundation for essentially all molecular electronic structure calculations.

---

*IQCP Lab Pack #1 v1.0 | Interactive Quantum Chemistry Playground | https://iqcp.dev*
