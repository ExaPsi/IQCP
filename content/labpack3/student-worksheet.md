# Lab Pack #3: Basis Functions, Integrals, and Electron Density

**Lab Pack:** 3 - Computational Layers of Quantum Chemistry
**Version:** 1.0
**Last Updated:** 2026-04-05
**Estimated Time:** 60 minutes (sections total ~60 min; synthesis adds ~5 min buffer)
**Prerequisite:** Lab Packs #1 and #2 completed; familiarity with Modules C-E

---

## Introduction

Welcome to Lab Pack #3 of the Interactive Quantum Chemistry Playground (IQCP). In Lab Packs #1 and #2, you built a foundation for quantum chemistry computation: you explored the numerical engines (Boys functions, Rys quadrature, SCF convergence), saw molecules in 3D, scanned potential energy surfaces, and visualized molecular orbitals as isosurfaces.

But there are critical computational layers you have not yet examined directly. When the SCF engine builds a Fock matrix, where do the numbers come from? What do the basis functions actually look like? How does the computer assemble the electron density from a density matrix and a set of basis functions?

In this lab, you will peel back three layers of the computation:

1. **Basis functions** -- the mathematical building blocks that approximate atomic orbitals. You will see what a contracted Gaussian looks like, compare minimal and split-valence basis sets, and confront the question: is a basis function the same as an atomic orbital?

2. **One- and two-electron integrals** -- the matrix elements that encode physics (overlap, kinetic energy, nuclear attraction, electron repulsion). You will inspect matrix heatmaps, trace how the Fock matrix is assembled step by step, and see how the density matrix weights electron repulsion integrals.

3. **Electron density** -- the observable quantity that connects the abstract wavefunction to physical reality. You will visualize total density isosurfaces, examine 2D cross-sections, and use difference density maps to see where electrons accumulate when atoms form bonds.

Throughout this lab, you will follow a **Predict-Observe-Explain** (POE) approach. Before each observation, you will write a prediction. Predictions that turn out to be wrong are just as valuable as correct ones -- they reveal where your mental model needs updating.

---

## Learning Objectives

By completing this lab, you will be able to:

1. **LO13 (Understand):** Describe the components of a contracted Gaussian basis function and explain how it differs from a single Gaussian primitive.

2. **LO14 (Analyze):** Compare minimal, split-valence, and polarized basis sets by examining radial profiles; explain how basis set flexibility affects SCF energy and computational cost.

3. **LO15 (Analyze):** Identify the physical meaning of overlap (S), kinetic energy (T), and nuclear attraction (V) integrals, and predict how interatomic distance affects their magnitudes.

4. **LO16 (Analyze):** Trace Fock matrix construction from the core Hamiltonian H^core and the two-electron matrix G(P); explain the density matrix's role in weighting electron repulsion integrals.

5. **LO17 (Analyze):** Interpret electron density isosurfaces and 2D cross-sections; predict where density accumulates in bonding regions.

6. **LO18 (Analyze):** Analyze difference density maps to identify charge accumulation and depletion regions, and connect these patterns to bonding concepts.

---

## What You Will Need

- A modern web browser (Chrome, Firefox, Safari, or Edge)
- This worksheet (print or digital)
- Access to the IQCP web application at: **https://iqcp.dev**
- Approximately 60 minutes of uninterrupted time

---

## What to Submit

At the end of this lab, you should have:

- [ ] Written answers to all numbered questions (26 total: Q1.1-Q1.8, Q2.1-Q2.8, Q3.1-Q3.8, Q4.1-Q4.2)
- [ ] Screenshot deliverables:
  - Section 1: STO-3G H 1s radial profile with primitives; STO-3G vs. 6-31G comparison for O 2s
  - Section 2: H2 overlap matrix heatmap with cell selected; H2O overlap matrix; Fock build trace
  - Section 3: H2O total density isosurface; H2O density cross-section; H2 difference density
- [ ] One exported density run artifact (JSON)
- [ ] Your responses to the synthesis questions (Section 4)

Your instructor will provide specific submission instructions.

---

## How to Navigate IQCP

Throughout this worksheet, you will find step-by-step navigation instructions in **indented blocks** telling you exactly how to configure each IQCP module. Follow the numbered steps in order.

**Lab Pack #3 uses three different IQCP modules:**
- **Module A** (Basis Explorer) at https://iqcp.dev/v1/basis -- for Section 1
- **Module B** (Integral Inspector) at https://iqcp.dev/v1/integrals -- for Section 2
- **Module E** (SCF Sandbox) at https://iqcp.dev/v1/scf -- for Section 3

**Tip:** You can keep multiple browser tabs open -- one for each module -- so you can switch between them quickly.

---

## Section 1: Basis Function Anatomy and Comparison (~15 min)

**Target learning objectives:** LO13 (contracted Gaussian anatomy), LO14 (basis set comparison)

In Lab Packs #1 and #2, you used basis sets like STO-3G without examining the basis functions themselves. What do these functions actually look like? In this section, you will visualize radial profiles, explore the anatomy of a contracted Gaussian, and compare different basis sets.

---

### Activity 1.1: Exploring the H 1s Radial Profile in STO-3G

**Predict:** Before viewing the radial profile, write your prediction:

> **Q1.1:** A single Gaussian function has the form g(r) = N exp(-alpha r^2). It is a smooth bell-shaped curve. The STO-3G basis set for hydrogen uses three of these Gaussians combined (contracted) into one basis function. Predict: will the contracted function look like (a) a wider bell curve, (b) a sharper, more peaked curve, or (c) three separate bumps? Sketch what you think the radial profile will look like.

*Your prediction and sketch:*



**Observe:**

> 1. Navigate to **https://iqcp.dev/v1/basis** (Module A: Basis Explorer).
> 2. In the periodic table grid on the left, click **H** (Hydrogen).
> 3. From the **Basis Set** dropdown, select **STO-3G**.
> 4. In the **Shell List** below, click the **1s** shell row.
> 5. Examine the radial profile plot on the right. Notice the three individual primitive Gaussians (dashed lines) and the total contracted function (solid line).

**Explain:**

> **Q1.2:** Describe the radial profile you observe. How does the shape of the contracted function compare to the individual primitives? Which primitive contributes most near the nucleus? Which contributes most at larger distances?

*Your answer:*



---

### Activity 1.2: Understanding Contraction

**Explain (continued):**

> **Q1.3:** STO-3G uses 3 primitives contracted together to approximate a Slater-type orbital (which has a cusp at r = 0 and decays exponentially). Why not just use a single Gaussian? What does adding more primitives with different exponents achieve? What is the tradeoff?

*Your answer:*



> **Q1.4: Misconception Check** -- A student says: "A basis function is the same thing as an atomic orbital." Do you agree with this statement? Consider: does the oxygen 2p basis function in STO-3G have exactly the same shape as a hydrogen-like 2p orbital? What is the relationship between basis functions and atomic orbitals?

*Your answer:*



**Observe:**

> 1. With **H / STO-3G / 1s** still selected in Module A, look below the radial profile plot for the **Exponent Sliders** panel. (These sliders appear automatically when a shell is selected and comparison mode is off.)
> 2. Drag one of the exponent sliders to modify a primitive exponent. Observe how the radial profile changes in real time. The modified profile appears overlaid on the original.

> Record: When you increase the exponent of the tightest (largest-exponent) primitive, does the contracted function become more peaked or more diffuse near the nucleus?

*Your observation:*



---

### Activity 1.3: Comparing STO-3G and 6-31G

**Predict:** Before comparing basis sets, write your prediction:

> **Q1.5:** The 6-31G basis set is a "split-valence" basis -- it uses two sets of contractions for valence shells instead of one. Predict: will the oxygen 2s radial profile in 6-31G extend further from the nucleus, closer to the nucleus, or be about the same as in STO-3G?

*Your prediction:*



**Observe:**

> 1. In Module A (https://iqcp.dev/v1/basis), click **O** (Oxygen) in the periodic table grid.
> 2. From the **Basis Set** dropdown, select **STO-3G**.
> 3. In the Shell List, click the **2s** shell row.
> 4. Below the radial profile plot, find the **Advanced** section and click the **Compare Basis Sets** button to activate comparison mode.
> 5. In the comparison controls that appear, add **6-31G** to the comparison list. If a shell type dropdown appears, ensure **2s** is selected.
> 6. Examine the overlaid radial profiles for STO-3G and 6-31G side by side.

**Explain:**

> **Q1.6:** Describe the differences between the O 2s radial profiles in STO-3G and 6-31G. How many contractions does each basis set use for the valence 2s? What additional flexibility does the split-valence approach provide?

*Your answer:*



> **Q1.7:** What does "split-valence" mean in 6-31G? The name encodes the basis set structure: "6" primitives for the core, "31" means the valence is split into a 3-primitive contraction and a 1-primitive (uncontracted) function. Why split the valence shells but not the core?

*Your answer:*



> **Q1.8: Misconception Check** -- A student claims: "Adding more basis functions always makes the answer more accurate, so we should always use the biggest basis set available." Evaluate this claim. What practical consideration limits the benefit of adding more basis functions? (Hint: think about computational cost scaling and the concept of diminishing returns.)

*Your answer:*



---

### Checkpoint: Section 1 Deliverables

Before moving on, verify that you have completed:

- [ ] Screenshot 1.A: STO-3G H 1s radial profile showing primitive decomposition (dashed lines) and contracted function (solid line)
- [ ] Screenshot 1.B: Basis set comparison -- STO-3G vs. 6-31G for O 2s shell
- [ ] Written answers to Q1.1 through Q1.8

---

## Section 2: Integral Inspection and Fock Tracing (~20 min)

**Target learning objectives:** LO15 (integral physical meaning), LO16 (Fock matrix construction)

In Lab Pack #1, you saw the SCF procedure use matrices like H^core and the Fock matrix F, but you did not examine where these matrix elements come from. In this section, you will open the hood: inspect overlap, kinetic, and nuclear attraction integrals; trace how the Fock matrix is built step by step; and see how electron repulsion integrals are weighted by the density matrix.

---

### Activity 2.1: Exploring Integral Matrices for H2

**Predict:** Before viewing any integrals, write your prediction:

> **Q2.1:** Three fundamental one-electron integral matrices appear in quantum chemistry:
> - **S** (overlap): measures spatial overlap between basis functions
> - **T** (kinetic energy): matrix elements of the kinetic energy operator
> - **V** (nuclear attraction): matrix elements of the nuclear attraction potential
>
> For the H2 molecule (2 basis functions, one 1s per atom), which matrix do you predict will have the largest off-diagonal element (in absolute value): S, T, or V? Why?

*Your prediction:*



**Observe:**

> 1. Navigate to **https://iqcp.dev/v1/integrals** (Module B: Integral Inspector).
> 2. Ensure **Preset** mode is selected (the left button in the input mode toggle).
> 3. From the **Molecule** dropdown, select **H2**.
> 4. From the **Basis Set** dropdown, select **STO-3G**.
> 5. The integral matrices will compute automatically. In the matrix tab bar, click **S** (overlap) if it is not already selected.
> 6. View the 2x2 overlap matrix heatmap. Click on the off-diagonal cell S(1,2) to see the primitive breakdown panel below, which shows how each pair of Gaussian primitives contributes to the total overlap integral.

> Record the value of S(1,2): _______________

**Observe:**

> 1. With H2 / STO-3G still loaded in Module B, click the **T** tab in the matrix tab bar to view the kinetic energy matrix.
> 2. Then click the **V** tab to view the nuclear attraction matrix. Compare the magnitudes of the off-diagonal elements across S, T, and V.

> Record the off-diagonal values:
> - T(1,2): _______________
> - V(1,2): _______________

**Explain:**

> **Q2.2:** Which matrix (S, T, or V) had the largest off-diagonal element for H2? Why does this matrix dominate? Think about what each integral physically measures.

*Your answer:*



> **Q2.3:** If the H-H distance were doubled (from 1.4 to 2.8 bohr), predict how the overlap integral S(1,2) would change. Would it increase, decrease, or stay the same? Justify your answer using the definition of the overlap integral: S_ij = integral phi_i(r) phi_j(r) dr, where the basis functions are centered on different atoms.

*Your answer:*



---

### Activity 2.2: Exploring Integral Matrices for H2O

**Predict:** Before viewing the H2O integrals, write your prediction:

> The H2O molecule has 7 basis functions in STO-3G: 5 on oxygen (1s, 2s, 2px, 2py, 2pz) and 1 on each hydrogen (1s). Before viewing: which pairs of basis functions do you expect to have the largest overlaps -- bonded pairs (O-H) or non-bonded pairs (H-H)?

*Your prediction:*



**Observe:**

> 1. In Module B (https://iqcp.dev/v1/integrals), change the **Molecule** dropdown to **H2O** (keep **STO-3G** as the basis set).
> 2. Click the **S** tab in the matrix tab bar to view the overlap matrix.
> 3. Notice the block structure of the 7x7 matrix: the upper-left 5x5 block corresponds to oxygen-oxygen interactions, and the off-diagonal blocks show oxygen-hydrogen interactions.

**Explain:**

> **Q2.4:** Examine the H2O overlap matrix. Identify the matrix elements that correspond to:
> - O-H bonded pairs (which row/column combinations?)
> - H-H non-bonded pair (which element?)
>
> Record approximate values for one O-H overlap and the H-H overlap. Are these consistent with what you know about the relative distances between these atoms?

*Your answer:*

| Atom Pair | Matrix Element | Approximate Value |
|-----------|----------------|-------------------|
| O-H (bonded) | S( , ) | |
| H-H (non-bonded) | S( , ) | |

---

### Activity 2.3: Fock Build Tracing

The Fock matrix F is the central quantity in the SCF procedure. It is built from two pieces:

```
F = H^core + G(P)
```

where H^core = T + V (the one-electron integrals you just examined) and G(P) encodes electron-electron repulsion, weighted by the density matrix P.

**Observe:**

> 1. In Module B (https://iqcp.dev/v1/integrals), change the **Molecule** dropdown back to **H2** (with **STO-3G** selected).
> 2. Wait for the integral matrices to load. Then scroll down past the heatmap to the **Advanced** section.
> 3. Click the **Fock Build Tracer** button to expand the Fock build panel.
> 4. Click the **Run SCF + Decompose Fock Matrix** button. This runs a tight-convergence RHF calculation to obtain the density matrix, then decomposes the Fock matrix.
> 5. Once the decomposition completes, step through the Fock matrix construction using the step controls:
>    - **Step 1:** H^core = T + V (one-electron terms)
>    - **Step 2:** G(P) (two-electron terms weighted by density matrix)
>    - **Step 3:** F = H^core + G(P) (complete Fock matrix)
> 6. At each step, examine the matrix displayed and note how the values change.

**Explain:**

> **Q2.5:** Write out the formula F = H^core + G(P) in your own words. What physical interactions does H^core capture? What does G(P) add?

*Your answer:*



> **Q2.6:** The G(P) matrix depends on the density matrix P, which changes at each SCF iteration. Why does this mean the Fock matrix must be rebuilt at every iteration? What would happen if you froze the density matrix and stopped updating it?

*Your answer:*



> **Q2.7:** In the Fock build trace for H2, examine one electron repulsion integral (ERI) quartet that contributes to F(1,1). The contribution is weighted by the density matrix element P(lambda, sigma). How does the density matrix determine which ERI contributions matter most? (Hint: if P(lambda, sigma) is near zero, what happens to that ERI's contribution?)

*Your answer:*



> **Q2.8:** A student says: "The Fock matrix is just H^core plus the electron repulsion integrals." What is missing from this statement? (Hint: there are two things missing -- one involves the density matrix, and one involves exchange.)

*Your answer:*



---

### Checkpoint: Section 2 Deliverables

Before moving on, verify that you have completed:

- [ ] Screenshot 2.A: H2 overlap matrix heatmap with off-diagonal cell selected, showing primitive breakdown panel
- [ ] Screenshot 2.B: H2O overlap matrix heatmap showing block structure
- [ ] Screenshot 2.C: Fock build trace at step 2 (G(P) contribution) for H2
- [ ] Written answers to Q2.1 through Q2.8

---

## Section 3: Electron Density and Difference Density (~25 min)

**Target learning objectives:** LO17 (density interpretation), LO18 (difference density analysis)

You have now seen the basis functions (Section 1) and the integrals that encode physics (Section 2). The SCF procedure combines all of these to produce a density matrix P and a set of orbital energies. But what does the electron density actually look like in space? In this section, you will visualize the total electron density as an isosurface and cross-section, then use difference density maps to see where electrons move when atoms form bonds.

---

### Activity 3.1: Total Density Isosurface for H2O

**Predict:** Before viewing the density, make a prediction:

> **Q3.1:** Sketch where you expect the electron density to be highest for H2O. Consider: oxygen has 8 electrons and each hydrogen has 1 -- where will the density concentrate? Will it be evenly distributed, or will one region dominate?

*Your prediction and sketch:*



**Observe:**

> 1. Navigate to **https://iqcp.dev/v1/scf** (Module E: SCF Sandbox).
> 2. In the **Single Point** tab (the default), select **H2O** from the **Molecule** dropdown and **STO-3G** from the **Basis Set** dropdown.
> 3. Click **Run SCF** and wait for the calculation to converge (look for a green convergence indicator).
> 4. After convergence, scroll down in the left sidebar to find the **Density Visualization** panel. Check the **Show** checkbox to enable it.
> 5. Ensure the **Density Mode** is set to **Total** (the left button).
> 6. The total density isosurface appears in the 3D viewer at the default isovalue (0.05 e/bohr^3).
> 7. Rotate the molecule in the 3D viewer. Notice the shape of the isosurface: where is it largest? Where does it extend the furthest?

**Explain:**

> Does the isosurface shape match your prediction from Q3.1? Where is the density highest -- near the oxygen nucleus, in the O-H bonding region, or near the hydrogens?

*Your observation:*



---

### Activity 3.2: Density Cross-Section

**Observe:**

> 1. With H2O density still displayed in Module E, scroll down in the Density Visualization panel to find the **Cross-Section Plane** selector.
> 2. Click the **XZ** button to select the XZ plane (the molecular plane for H2O).
> 3. A 2D density cross-section plot appears below the 3D viewer. Examine the color-coded density map.
> 4. Adjust the **Plane Position** slider to move the cross-section slice through the molecule.
> 5. Notice the sharp peaks near the nuclear positions and the density between the bonded atoms.

**Explain:**

> **Q3.2:** In the 2D cross-section, there are sharp density peaks near the nuclear positions. Why does the electron density spike near the nuclei? (Hint: think about the Coulomb attraction between electrons and nuclei.)

*Your answer:*



> **Q3.3:** Compare the density in the O-H bonding region to the density in the H-H non-bonding region (across from the oxygen). Is there more density in the bonding region? What does this tell you about how electrons participate in bonding?

*Your answer:*



---

### Activity 3.3: Isovalue Exploration

**Observe:** Remain in the total density view for H2O in Module E. In the **Density Visualization** panel, locate the **Density Isovalue** slider.

> 1. Set the isovalue to 0.05 (default) and note the surface size in the 3D viewer
> 2. Decrease the isovalue to 0.01 -- observe how the surface changes
> 3. Increase the isovalue to 0.10 -- observe again
> 4. Return to 0.05

**Explain:**

> **Q3.4: Misconception Check** -- When you decreased the isovalue, the density isosurface expanded to enclose a larger volume. Does this mean there are more electrons? Explain what the isovalue threshold actually represents. (The isosurface shows all points where rho(r) = isovalue. Decreasing the isovalue means you are showing regions of lower density that extend further from the nuclei.)

*Your answer:*



---

### Activity 3.4: Difference Density for H2

Now you will use the difference density to visualize what happens to the electron density when two atoms come together to form a molecule. The difference density is defined as:

```
Delta-rho(r) = rho_molecule(r) - rho_promolecule(r)
```

where the promolecule is the sum of spherically-averaged atomic densities placed at the molecular positions -- what the density would look like if the atoms were nearby but not interacting.

**Predict:** Before viewing the difference density, make a prediction:

> **Q3.5:** When two hydrogen atoms form an H2 molecule, electrons rearrange. Predict: where will electrons accumulate (compared to the promolecule)? Where will they deplete? Sketch your prediction on a simple diagram of two H atoms.

*Your prediction and sketch:*



**Observe:**

> 1. In Module E (https://iqcp.dev/v1/scf), change the **Molecule** dropdown to **H2** (keep **STO-3G**).
> 2. Click **Run SCF** and wait for convergence.
> 3. In the **Density Visualization** panel, check **Show** and switch the **Density Mode** to **Difference** (the right button).
> 4. You will see two types of isosurface in the 3D viewer:
>    - **Green solid (accumulation):** regions where Delta-rho > 0 -- electrons have moved INTO this region upon bond formation
>    - **Red translucent (depletion):** regions where Delta-rho < 0 -- electrons have moved AWAY from this region
> 5. If the surfaces are not visible, try lowering the **Difference Isovalue** slider (e.g., to 0.005 or 0.002).
> 6. Rotate the molecule and identify where each region appears.

**Explain:**

> Based on your observation: where do electrons accumulate when two H atoms form H2? Where are they depleted? Is this consistent with the idea that a covalent bond involves electron sharing between atoms?

*Your observation:*



---

### Activity 3.5: Difference Density for H2O

**Observe:**

> 1. In Module E, change the **Molecule** dropdown to **H2O** (keep **STO-3G**).
> 2. Click **Run SCF** and wait for convergence.
> 3. In the **Density Visualization** panel, ensure **Show** is checked and **Density Mode** is set to **Difference**.
> 4. View the H2O difference density and compare the pattern to what you saw for H2.

**Explain:**

> **Q3.6:** The difference density uses a "promolecule" reference -- the density of non-interacting atoms placed at the molecular geometry. What is a promolecule, and why is it the right reference for computing the difference density? (Why not use the density of isolated atoms at infinite separation?)

*Your answer:*



> **Q3.7:** Compare the H2 and H2O difference density maps. In H2O, where does charge accumulate relative to the promolecule? Do you observe any asymmetry in the accumulation pattern along the O-H bonds? Is the pattern consistent with the O-H bond being polar (oxygen is more electronegative than hydrogen)?

*Your answer:*



> **Q3.8:** A student says: "The solid regions in the difference density map contain all the bonding electrons; the translucent regions contain none." Correct this statement. What do the solid and translucent regions actually represent? (Hint: the difference density shows changes relative to the promolecule, not absolute electron counts.)

*Your answer:*



---

### Checkpoint: Section 3 Deliverables

Before moving on:

1. Export a run artifact: click the **Export** button in the Module E header toolbar (top-right, next to "Share") and save as `density-artifact.json`
2. Take the required screenshots

- [ ] Screenshot 3.A: H2O total density isosurface
- [ ] Screenshot 3.B: H2O density 2D cross-section (XZ plane)
- [ ] Screenshot 3.C: H2 difference density showing accumulation (solid) and depletion (translucent) regions
- [ ] Exported `density-artifact.json`
- [ ] Written answers to Q3.1 through Q3.8

---

## Section 4: Synthesis and Reflection (~5 min)

Now connect what you have learned across all three sections.

### Connecting Basis Functions, Integrals, and Density

The three topics in this lab form a computational pipeline:

1. **Basis functions** (Section 1) are the mathematical building blocks -- contracted Gaussians with specific exponents and coefficients.
2. **Integrals** (Section 2) are computed from pairs (or quartets) of basis functions. They encode the physics: overlap, kinetic energy, nuclear attraction, and electron repulsion.
3. **Electron density** (Section 3) is the final physical observable, constructed from the density matrix P and the basis functions: rho(r) = sum P_ij phi_i(r) phi_j(r).

Each layer builds on the previous one. The choice of basis set determines the integrals, the integrals determine the Fock matrix and SCF solution, and the SCF solution determines the electron density.

---

**Q4.1:** Trace the computational path from basis functions to electron density. Suppose you changed the basis set from STO-3G to 6-31G for H2O. At each layer of the pipeline, describe one specific thing that would change:

- Basis functions: _______________________________________________
- Integrals: _______________________________________________
- Fock matrix / SCF result: _______________________________________________
- Electron density: _______________________________________________

**Q4.2:** You explored three levels of representation today: mathematical (exponent tables, matrix elements, density values), graphical (radial profiles, heatmaps, cross-section color maps), and spatial (3D isosurfaces). Which representation was most useful for understanding bonding, and why? Was there a concept that only became clear when you saw it in a particular representation?

*Your answer:*



---

## Final Deliverables Checklist

Before submitting, verify that you have completed:

**Section 1: Basis Function Exploration (8 questions)**
- [ ] **Q1.1:** Predicted contracted Gaussian shape
- [ ] **Q1.2:** Described STO-3G H 1s radial profile and primitive contributions
- [ ] **Q1.3:** Explained why STO-3G uses 3 primitives
- [ ] **Q1.4:** Evaluated "basis function = atomic orbital" claim (misconception check)
- [ ] **Q1.5:** Predicted 6-31G vs. STO-3G difference for O 2s
- [ ] **Q1.6:** Described STO-3G vs. 6-31G radial profile differences
- [ ] **Q1.7:** Explained split-valence meaning in 6-31G
- [ ] **Q1.8:** Evaluated "more basis functions = always more accurate" claim (misconception check)

**Section 2: Integral Inspection and Fock Tracing (8 questions)**
- [ ] **Q2.1:** Predicted which matrix has largest off-diagonal element
- [ ] **Q2.2:** Identified dominant matrix and explained why
- [ ] **Q2.3:** Predicted effect of doubling H-H distance on S(1,2)
- [ ] **Q2.4:** Identified O-H and H-H overlaps in H2O matrix with values
- [ ] **Q2.5:** Explained F = H^core + G(P) in own words
- [ ] **Q2.6:** Explained why Fock matrix must be rebuilt each iteration
- [ ] **Q2.7:** Explained density matrix weighting of ERI contributions
- [ ] **Q2.8:** Identified what is missing from "F = H^core + ERIs"

**Section 3: Electron Density and Difference Density (8 questions)**
- [ ] **Q3.1:** Predicted density distribution for H2O
- [ ] **Q3.2:** Explained density peaks near nuclei
- [ ] **Q3.3:** Compared bonding vs. non-bonding density regions
- [ ] **Q3.4:** Explained isovalue meaning (misconception check)
- [ ] **Q3.5:** Predicted and observed H2 difference density pattern
- [ ] **Q3.6:** Explained promolecule reference and its purpose
- [ ] **Q3.7:** Compared H2 and H2O difference density patterns
- [ ] **Q3.8:** Corrected misconception about difference density regions

**Section 4: Synthesis (2 questions)**
- [ ] **Q4.1:** Traced computational path from basis functions to density
- [ ] **Q4.2:** Reflected on which representation best supported understanding

**Screenshots and Artifacts:**
- [ ] Screenshot 1.A: STO-3G H 1s radial profile with primitive decomposition
- [ ] Screenshot 1.B: STO-3G vs. 6-31G comparison for O 2s
- [ ] Screenshot 2.A: H2 overlap matrix with cell selected and primitive breakdown
- [ ] Screenshot 2.B: H2O overlap matrix heatmap
- [ ] Screenshot 2.C: Fock build trace at step 2 for H2
- [ ] Screenshot 3.A: H2O total density isosurface
- [ ] Screenshot 3.B: H2O density 2D cross-section (XZ plane)
- [ ] Screenshot 3.C: H2 difference density (accumulation and depletion)
- [ ] Exported `density-artifact.json`

---

## Appendix A: Troubleshooting

**Problem:** The radial profile plot in Module A is empty or does not show primitives
**Solution:** Verify that the correct element and basis set are selected. The primitive decomposition (dashed lines) appears when the shell is selected in the shell table. Click on a shell row to select it.

**Problem:** The comparison mode in Module A does not show two profiles
**Solution:** Ensure comparison mode is toggled ON and that a second basis set is selected from the comparison dropdown. Both basis sets must support the selected element.

**Problem:** The integral matrix heatmap in Module B shows all zeros or does not appear
**Solution:** Ensure a molecule and basis set are selected from the dropdowns. Integrals compute automatically once a valid molecule/basis combination is selected. If using a non-standard combination (e.g., CH4 with 6-31G*), the computation may take a moment as integrals are computed on-the-fly.

**Problem:** Clicking a cell in the integral matrix does not show the primitive breakdown
**Solution:** Click directly on the colored cell in the heatmap. The primitive breakdown panel appears below or to the right of the matrix view. If the cell is on the diagonal, you will see self-overlap primitives.

**Problem:** The Fock build trace does not show step-by-step progression
**Solution:** First click the **Fock Build Tracer** button in the Advanced section to expand the panel. Then click **Run SCF + Decompose Fock Matrix** to compute the decomposition. Once complete, use the step controls (numbered buttons or slider) to step through: Step 1 shows H^core, Step 2 shows G(P), and Step 3 shows the final F = H^core + G(P).

**Problem:** The density isosurface does not appear or is very faint
**Solution:** Ensure the SCF calculation has converged (check for a green convergence indicator). If the isosurface is not visible, try decreasing the isovalue using the slider. The default isovalue of 0.05 should produce a visible surface for most molecules.

**Problem:** The difference density shows no accumulation or depletion regions
**Solution:** Difference density values are much smaller than total density values. Try decreasing the difference density isovalue (e.g., to 0.005 or 0.002). The default isovalue for difference density is typically lower than for total density.

**Problem:** The IQCP page does not load or shows an error
**Solution:** Clear your browser cache and refresh the page. Ensure you are using a modern browser (Chrome, Firefox, Safari, or Edge) with JavaScript enabled. If a module appears stuck on "Initializing WASM module...", wait a few seconds for the compute engine to load.

---

## Appendix B: Glossary

**Basis function:** A mathematical function used to approximate molecular orbitals in quantum chemistry. In practice, basis functions are contracted Gaussians -- linear combinations of primitive Gaussians with fixed exponents and coefficients.

**Contracted Gaussian:** A linear combination of primitive Gaussian functions: phi(r) = sum d_i g(alpha_i, r). The exponents alpha_i and coefficients d_i are pre-determined and fixed during the SCF calculation.

**Core Hamiltonian (H^core):** The one-electron part of the Fock matrix, containing kinetic energy and nuclear attraction contributions: H^core = T + V. It does not depend on the electron density.

**Difference density (Delta-rho):** The difference between the molecular electron density and the promolecule density: Delta-rho(r) = rho_molecule(r) - rho_promolecule(r). Positive values indicate electron accumulation; negative values indicate depletion.

**Electron density (rho):** The probability distribution for finding an electron at position r, computed from the density matrix and basis functions: rho(r) = sum P_ij phi_i(r) phi_j(r). Units: electrons per cubic bohr.

**Electron repulsion integral (ERI):** A four-center integral (mu nu | lambda sigma) describing the Coulomb interaction between two charge distributions, each formed from a pair of basis functions.

**Fock matrix (F):** The effective one-electron Hamiltonian including electron-electron interactions: F = H^core + G(P). It depends on the density matrix P and must be rebuilt at each SCF iteration.

**Isovalue:** A threshold value used to define an isosurface. For density visualization, the isosurface shows all points where rho(r) equals the isovalue. The electron density extends beyond this surface -- the boundary is a chosen threshold, not a physical edge.

**Kinetic energy integral (T):** A matrix element of the kinetic energy operator between two basis functions: T_ij = <phi_i | -1/2 nabla^2 | phi_j>. Measures the kinetic energy associated with electron motion.

**Nuclear attraction integral (V):** A matrix element of the nuclear attraction potential between two basis functions: V_ij = sum_A -Z_A <phi_i | 1/r_A | phi_j>. Measures the Coulomb attraction between electrons and nuclei.

**Overlap integral (S):** Measures the spatial overlap between two basis functions: S_ij = integral phi_i(r) phi_j(r) dr. Ranges from 0 (no overlap) to 1 (identical functions). Off-diagonal values decrease with increasing interatomic distance.

**Polarization function:** A higher-angular-momentum function added to a basis set to describe how electron distributions distort in a molecular environment. For example, the asterisk in 6-31G* indicates d-type polarization functions on non-hydrogen atoms.

**Primitive Gaussian:** A single Gaussian function g(alpha, r) = N exp(-alpha r^2), characterized by its exponent alpha and its center position. The exponent controls the width: large alpha gives a tight (narrow) function, small alpha gives a diffuse (wide) function.

**Promolecule:** A reference electron density formed by placing spherically-averaged free-atom densities at the molecular atomic positions, without allowing them to interact. Used as the baseline for computing difference density.

**Split-valence basis:** A basis set that uses multiple contractions for valence shells to provide additional flexibility for describing bonding. For example, 6-31G splits each valence shell into two functions (a 3-primitive contraction and a 1-primitive function), while using a single 6-primitive contraction for the core.

---

## Appendix C: References

1. Szabo, A. and Ostlund, N.S. (1996). *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*. Dover Publications.

2. Hehre, W.J., Stewart, R.F., and Pople, J.A. (1969). "Self-consistent molecular-orbital methods. I. Use of Gaussian expansions of Slater-type atomic orbitals." *Journal of Chemical Physics*, 51(6), 2657-2664.

3. Kozma, R. and Russell, J. (2005). "Students becoming chemists: Developing representational competence." In J.K. Gilbert (Ed.), *Visualization in Science Education*, 121-145. Springer.

4. Tsaparlis, G. and Papaphotis, G. (2009). "High-school students' conceptual difficulties and attempts at conceptual change: The case of basic quantum chemical concepts." *International Journal of Science Education*, 31(7), 895-930.

5. White, R.T. and Gunstone, R.F. (1992). *Probing Understanding*. Falmer Press.

6. Wiggins, G. and McTighe, J. (2005). *Understanding by Design* (2nd ed.). ASCD.

---

**Lab Pack #3 Complete**

Thank you for completing this guided exploration of basis functions, integrals, and electron density. The concepts you explored today -- how basis functions become integrals, how integrals build the Fock matrix, and how the SCF solution produces the electron density -- form the computational backbone of modern quantum chemistry. Every ab initio calculation you will ever encounter follows this same pipeline.

---

*IQCP Lab Pack #3 v1.0 | Interactive Quantum Chemistry Playground | https://iqcp.dev*
