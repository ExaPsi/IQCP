# Lab Pack #2: 3D Visualization, PES Scanning, and Orbital Exploration

**Lab Pack:** 2 - 3D Exploration, PES, and Orbitals
**Version:** 1.0
**Last Updated:** 2026-04-05
**Estimated Time:** 60 minutes (sections total ~65 min; most students finish in 50-60 min)
**Prerequisite:** Lab Pack #1 completed; familiarity with Module E (SCF Sandbox)

---

## Introduction

Welcome to Lab Pack #2 of the Interactive Quantum Chemistry Playground (IQCP). In Lab Pack #1, you explored the numerical foundations of quantum chemistry: Boys functions, Rys quadrature, and SCF convergence. You worked with matrices, energies, and convergence plots -- but the molecules themselves remained abstract.

In this lab, you will connect those numerical results to the physical and spatial reality of molecules. You will:

1. **See molecules in 3D** and connect their geometry to the overlap matrix elements you computed in Lab Pack #1.
2. **Scan potential energy surfaces** to discover why chemical bonds form and where the SCF method breaks down.
3. **Visualize molecular orbitals** as 3D isosurfaces and confront the question: does an orbital have a definite edge?

Throughout this lab, you will follow a **Predict-Observe-Explain** (POE) approach. Before each observation, you will make a prediction. This is not a test -- predictions that turn out to be wrong are just as valuable as correct ones, because they reveal where your mental model needs updating.

---

## Learning Objectives

By completing this lab, you will be able to:

1. **LO7 (Analyze):** Predict how total electronic energy changes as a function of internuclear distance for a diatomic molecule, and explain the physical origin of the potential energy minimum.

2. **LO8 (Analyze):** Distinguish bonding, nonbonding, and antibonding molecular orbitals from their isosurface shapes, and relate orbital character to molecular properties.

3. **LO9 (Apply):** Given an SCF overlap matrix and a 3D molecular view, identify which matrix element corresponds to a specific atom pair, and predict whether increasing the interatomic distance will increase or decrease that element.

4. **LO10 (Evaluate):** Evaluate the tradeoff between computational cost and accuracy by comparing SCF results across different basis sets, and explain why larger basis sets yield lower energies (variational principle).

5. **LO11a (Analyze):** Identify that RHF produces physically incorrect behavior at large bond distances and recognize this as a limitation of the method.

6. **LO11b (Evaluate, graduate extension):** Explain the RHF dissociation failure in terms of single-determinant limitations.

7. **LO12 (Understand):** Explain what an isovalue represents and predict how the isosurface shape changes as the isovalue is adjusted.

---

## What You Will Need

- A modern web browser with WebGL support (Chrome, Firefox, Safari, or Edge)
- This worksheet (print or digital)
- Access to the IQCP web application at: **https://iqcp.dev**
- Approximately 60 minutes of uninterrupted time

---

## What to Submit

At the end of this lab, you should have:

- [ ] Written answers to all numbered questions (26 total: Q1.1-Q1.8, Q2.1-Q2.8, Q3.1-Q3.8, Q4.1-Q4.2)
- [ ] Screenshot deliverables: 3D molecular view (Section 1), PES curve (Section 2), orbital isosurfaces (Section 3)
- [ ] Two exported run artifacts (PES scan artifact, orbital visualization artifact)
- [ ] Your responses to the synthesis questions (Section 4)

Your instructor will provide specific submission instructions.

---

## How to Navigate IQCP

All activities in this lab take place in **Module E (SCF Sandbox)**. To get started:

1. Open [https://iqcp.dev/v1/scf](https://iqcp.dev/v1/scf) in your browser.
2. Module E has four **workflow tabs** across the top: **Single Point**, **Optimize**, **PES Scan**, and **Compare**. Each section of this worksheet will tell you which tab to use.
3. The left side panel contains molecule and basis set selectors, SCF settings, and action buttons. The right side displays the **3D viewer**, **results panels**, and **plots**.
4. After running an SCF calculation, the Results panel offers two view modes: **Explain** (educational summary) and **Internals** (overlap, Fock, and density matrices).

Each activity step includes explicit instructions for which molecule, basis set, tab, and settings to use. Follow them in order.

---

## Section 1: 3D Molecular Exploration (~15 min)

**Target learning objectives:** LO9 (spatial-symbolic bridging)

In Lab Pack #1, you ran SCF calculations and examined overlap matrices, Fock matrices, and orbital energies. But which matrix element corresponds to which pair of atoms? In this section, you will connect the numbers to the molecule.

### Step 1.1: Predicting Molecular Geometry

Before looking at any 3D visualization, make a prediction.

**Predict:**

**Q1.1:** Sketch the geometry of the water molecule (H2O) in the space below. Label each atom. What bond angle do you predict between the two O-H bonds? (Hint: think about VSEPR or what you know from general chemistry.)

*Your prediction:*

Bond angle: _______________ degrees

Sketch:




### Step 1.2: Observing the 3D Structure

Now open the 3D viewer to see the actual geometry.

**Open Module E (SCF Sandbox):** Navigate to [https://iqcp.dev/v1/scf](https://iqcp.dev/v1/scf). In the **Single Point** tab, select **H2O** from the molecule dropdown and **STO-3G** as the basis set. Click **Run SCF**. Once the calculation completes, the 3D viewer panel ("3D Structure & Orbitals") will appear on the right. Check the **Labels** checkbox in the viewer header to display atom labels.

**Observe:**

- Rotate the molecule by clicking and dragging. Zoom with the scroll wheel.
- Notice the CPK coloring: oxygen is red, hydrogen is white.
- Toggle atom labels ON if they are not already visible.
- Examine the bond angle visually.

**Explain:**

**Q1.2:** Does the 3D geometry match your prediction from Q1.1? What is the approximate bond angle shown? If your prediction was different, what assumption led to the discrepancy?

*Your answer:* _______________________________________________

_______________________________________________

### Step 1.3: Running SCF and Examining the Overlap Matrix

Now run an SCF calculation and look at the overlap matrix alongside the 3D view.

**Using your H2O calculation from Step 1.2:** After the SCF calculation has converged, locate the **Results** panel below the convergence plot. Click the **Internals** toggle (next to "Explain") to switch to the internals view, which displays the overlap matrix **S** alongside other matrices.

**Predict:**

**Q1.3:** The overlap matrix S measures how much each pair of basis functions overlaps in space. Before examining the matrix: which pairs of atoms do you expect to have the largest overlap? Which pairs should have the smallest? (Consider which atoms are bonded and which are far apart.)

*Your prediction:*
- Largest overlap between: _______________________________________________
- Smallest overlap between: _______________________________________________

**Observe:**

Examine the overlap matrix S in the Internals view of the Results panel. The matrix elements S_ij describe the spatial overlap between basis functions i and j.

For H2O with STO-3G:
- Basis functions 1-5 are centered on oxygen (1s, 2s, 2px, 2py, 2pz)
- Basis function 6 is the hydrogen-1 1s orbital
- Basis function 7 is the hydrogen-2 1s orbital

**Q1.4:** Look at the overlap matrix elements. What is the approximate value of the overlap between:
- An oxygen basis function and a bonded hydrogen (e.g., S_2,6)? _______________
- The two hydrogen atoms (S_6,7)? _______________
- The diagonal elements (e.g., S_1,1)? _______________

**Explain:**

**Q1.5:** Why are the diagonal elements of the overlap matrix all equal to 1.0? What would it mean physically if a diagonal element were not 1.0?

*Your answer:* _______________________________________________

_______________________________________________

### Step 1.4: Comparing Molecules

Now switch to the H2 molecule for comparison.

**Switch to H2:** In the **Single Point** tab, change the molecule dropdown to **H2** (keep STO-3G as the basis set). Click **Run SCF**. Once converged, switch the Results panel to **Internals** mode to see the overlap matrix.

**Predict:**

**Q1.6:** H2 has only 2 basis functions (one 1s per atom). The overlap matrix is therefore 2x2. Before looking: predict the approximate value of the off-diagonal element S_1,2. Will it be closer to 0.0, 0.5, or 1.0?

*Your prediction:* _______________________________________________

**Observe:** Run SCF and examine the overlap matrix for H2.

**Explain:**

**Q1.7:** What is the actual value of S_1,2 for H2? How does this compare to the hydrogen-hydrogen overlap S_6,7 in water? Explain the difference in terms of the interatomic distances visible in the 3D views.

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

**Q1.8:** Based on what you have observed, predict: if you were to increase the H-H bond length in H2, would the off-diagonal overlap element S_1,2 increase or decrease? Explain your reasoning using both the 3D view and the mathematical meaning of the overlap integral.

*Your answer:* _______________________________________________

_______________________________________________

---

### Checkpoint: Section 1 Deliverables

Before moving on, verify that you have completed:

- [ ] Answers to Q1.1 through Q1.8
- [ ] Screenshot of the H2O 3D view showing atom labels and the bent geometry

---

## Section 2: Potential Energy Surface Scanning (~20 min)

**Target learning objectives:** LO7 (geometry-energy connection), LO11a (dissociation limit awareness), LO11b (dissociation explanation, graduate extension)

In Lab Pack #1, you computed the SCF energy at a single molecular geometry. But what happens to the energy as you change the geometry? In this section, you will scan the bond length of H2 and discover the potential energy surface (PES) -- the energy landscape that determines molecular structure.

### Step 2.1: Predicting the PES Shape

Before running any calculation, make a prediction.

**Predict:**

**Q2.1:** Sketch what you think the energy vs. bond length curve looks like for H2. On your sketch:
- Label the x-axis "Bond length R (bohr)" with a range of 0.5 to 5.0
- Label the y-axis "Energy (Hartree)"
- Mark where you think the energy minimum is
- Show what happens at very short R (atoms very close together)
- Show what happens at very large R (atoms far apart)

Sketch:




### Step 2.2: Running the PES Scan

Now run the actual scan and compare.

**Open the PES Scan tab:** In Module E, click the **PES Scan** tab in the workflow tab bar at the top of the page. Select **H2** as the molecule and **STO-3G** as the basis set in the left-side controls panel.

**Observe:**

1. For H2 (a diatomic molecule), the coordinate selector automatically selects a **Bond** scan between the two hydrogen atoms
2. Set the scan range: R_min = **0.5** bohr, R_max = **5.0** bohr, N points = **20**
3. Leave the scan mode as **Rigid** (the default)
4. Click **Start Scan** and watch the energy curve appear point by point in the PES plot on the right
5. Note the equilibrium marker that appears at the minimum
6. After the scan completes, try dragging the **Scan Geometry Viewer** slider below the plots to animate the molecule through the scan geometries

Record the following from the completed scan:

- Equilibrium bond length R_eq: _______________ bohr
- Equilibrium energy E_eq: _______________ Hartree
- Energy at R = 0.5 bohr: _______________ Hartree
- Energy at R = 5.0 bohr: _______________ Hartree

**Explain:**

**Q2.2:** Compare the computed PES curve to the sketch you made in Q2.1. What features did you predict correctly? What surprised you?

*Your answer:* _______________________________________________

_______________________________________________

### Step 2.3: Understanding the Energy Minimum

**Q2.3:** The potential energy curve has a minimum near R = 1.4 bohr. Why does this minimum exist? Consider two competing effects:
- At very short R, what happens to the nuclear-nuclear repulsion?
- At intermediate R, what stabilizes the system?

Explain in your own words why there is an optimal bond length.

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

**Q2.4:** The repulsive wall at short R rises steeply. Is this primarily due to (a) electron-electron repulsion, (b) nuclear-nuclear repulsion, or (c) both? Justify your answer.

*Your answer:* _______________________________________________

_______________________________________________

### Step 2.4: The Dissociation Limit Problem

Now focus on the right side of the PES curve, where R is large.

**Predict:**

**Q2.5:** When two hydrogen atoms are infinitely far apart, they should behave as two independent atoms. The energy of two isolated H atoms is 2 x E(H) = 2 x (-0.4666 Ha) = -0.9332 Ha. Look at the energy your PES scan gives at R = 5.0 bohr. Is it close to -0.9332 Ha?

*Your observation:*
- Energy at R = 5.0 bohr: _______________________________________________
- Expected for two isolated H atoms: -0.9332 Ha
- Difference: _______________________________________________

**Observe:**

Expand the **Computational Notes** section below the PES curve plot (click on it if collapsed). You should see a warning note about **RHF dissociation limits** when the scan includes large bond distances.

**Explain:**

**Q2.6:** The energy at large R does not approach the correct limit of two isolated hydrogen atoms. Instead, it is too high. This is a known limitation of the Restricted Hartree-Fock (RHF) method. Based on the computational note and what you know about RHF:
- What physical process does RHF fail to describe correctly at large R?
- Why might requiring alpha and beta electrons to share the same spatial orbitals cause problems when the bond is stretched?

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

> **Going Deeper (Graduate Extension, LO11b):** RHF uses a single Slater determinant, which means the wavefunction has the form |phi_1 phi_2 ... phi_N|. At large R, the correct H2 wavefunction is a superposition of sigma_g^2 and sigma_u^2 configurations. Explain why a single determinant cannot represent this superposition, and what type of method (e.g., multi-reference, CASSCF) would be needed.
>
> *Your answer (optional):* _______________________________________________
>
> _______________________________________________

### Step 2.5: Basis Set Comparison on the PES

**Compare basis sets:** Switch to the **Compare** tab. Select **H2** as the molecule, then check both **STO-3G** and **6-31G** in the basis set checklist. Click **Compare** to run SCF calculations with each basis set. Compare the converged energies in the results table.

Alternatively, you can run two separate PES scans: one with STO-3G and one with 6-31G, noting the equilibrium energy from each.

**Q2.7:** Compare the equilibrium energy of H2 with STO-3G vs. 6-31G (or a larger basis set if available). Which gives a lower (more negative) energy? Why does a larger basis set always give a lower energy? (Hint: think about the variational principle from Lab Pack #1.)

*Your answer:* _______________________________________________

_______________________________________________

**Q2.8:** A student claims: "A bigger basis set always gives a better answer, so we should always use the biggest possible basis." Do you agree? What practical considerations might limit this approach?

*Your answer:* _______________________________________________

_______________________________________________

> **Going Deeper (optional):** The PES Scan tab supports more than just bond-length scans on diatomics. For polyatomic molecules (e.g., H2O), you can select **Bond**, **Angle**, or **Dihedral** as the coordinate type, choose which atoms define the coordinate, and toggle between **Rigid** and **Relaxed** scan modes. Try scanning the O-H bond length or H-O-H angle in water if you have time.

---

### Checkpoint: Section 2 Deliverables

Before moving on:

1. Export the PES scan artifact: click the **Export** button (in the Module E header bar, near "Share") and save as `pes-artifact.json`
2. Take a screenshot of the completed PES curve

- [ ] Answers to Q2.1 through Q2.8
- [ ] Screenshot of the PES curve with equilibrium marker
- [ ] Exported `pes-artifact.json`

---

## Section 3: Orbital Visualization (~25 min)

**Target learning objectives:** LO8 (orbital interpretation), LO12 (isovalue interpretation)

So far in this lab, you have seen molecular geometry in 3D and explored how energy depends on geometry. But what do the computed molecular orbitals actually look like? In this section, you will visualize MO isosurfaces and learn what these shapes mean.

### Step 3.1: Predicting Orbital Shapes

Before viewing any orbitals, make predictions.

**Predict:**

**Q3.1:** H2O has 7 basis functions (STO-3G) and 5 occupied molecular orbitals. The lowest-energy MO (MO 1) is built primarily from the oxygen 1s atomic orbital.

Predict: will MO 1 look like (a) a small sphere centered on oxygen, (b) a shape spread across the entire molecule, or (c) two lobes pointing along the O-H bonds?

*Your prediction:* _______________________________________________

*Your reasoning:* _______________________________________________

### Step 3.2: Exploring Core and Valence Orbitals

**View H2O orbitals:** Return to the **Single Point** tab. Select **H2O** with **STO-3G** and click **Run SCF** (or use your earlier converged result). Once the calculation finishes, the **Orbitals** sidebar appears to the right of the 3D viewer. In the orbital selector dropdown, choose **MO 1** (the lowest-energy orbital).

**Observe:**

- Examine the isosurface shape of MO 1
- Note: the isosurface shows where the orbital wavefunction magnitude |psi| exceeds a threshold value (the isovalue)
- Rotate the molecule to see the orbital from different angles

**Q3.2:** Describe the shape of MO 1. Is it centered on one atom or spread across the molecule? Does your observation match your prediction from Q3.1? Why is this orbital almost unaffected by bonding?

*Your answer:* _______________________________________________

_______________________________________________

Now select the **HOMO** (highest occupied molecular orbital).

**Select the HOMO:** In the orbital selector dropdown, choose the orbital labeled **HOMO** (MO 5 for H2O with STO-3G). The isosurface will update in the 3D viewer.

**Q3.3:** Describe the shape of the H2O HOMO. How is it different from MO 1? Based on its shape, would you classify the HOMO as bonding, nonbonding, or antibonding? Explain your classification.

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

### Step 3.3: Bonding vs. Antibonding Orbitals in H2

To build clearer intuition about bonding character, examine the simplest case: H2.

**Switch to H2 orbitals:** In the **Single Point** tab, change the molecule to **H2** (STO-3G) and click **Run SCF**. Once converged, select **MO 1** in the orbital selector (the sigma_g bonding orbital).

**Predict:**

**Q3.4:** H2 has 2 MOs. MO 1 is the bonding orbital (sigma_g) and MO 2 is the antibonding orbital (sigma_u*). Before viewing MO 2: predict how the shape of the antibonding orbital will differ from the bonding orbital. Where do you expect to see a node (a surface where psi = 0)?

*Your prediction:* _______________________________________________

_______________________________________________

**Select MO 2:** In the orbital selector dropdown, switch to **MO 2** (the sigma_u* antibonding orbital). This is the LUMO for H2.

**Observe:**

- Notice that the isosurface has two distinct lobes
- Notice the rendering distinction: positive lobes are solid, negative lobes are translucent
- The node between the two lobes is where the orbital wavefunction is zero

**Explain:**

**Q3.5:** Compare the shapes of MO 1 (sigma_g) and MO 2 (sigma_u*) for H2. Fill in the comparison table:

| Feature | MO 1 (sigma_g) | MO 2 (sigma_u*) |
|---------|----------------|------------------|
| Number of lobes | | |
| Electron density between nuclei | High / Low | High / Low |
| Node between atoms? | Yes / No | Yes / No |
| Character | Bonding / Antibonding | Bonding / Antibonding |

### Step 3.4: The Isovalue Slider -- Does an Orbital Have an Edge?

This activity targets a common misconception: that orbitals have sharp boundaries.

**Return to H2O HOMO:** Switch back to **H2O** (STO-3G) in the molecule dropdown, run SCF, and select the **HOMO** in the orbital selector. Locate the **isovalue slider** below the orbital selector in the sidebar (default value: 0.03 a.u.).

**Predict:**

**Q3.6:** As you decrease the isovalue (make the threshold smaller), do you predict the isosurface will (a) shrink, (b) expand, or (c) stay the same size?

*Your prediction:* _______________________________________________

**Observe:**

1. Set the isovalue to 0.03 (default) and observe the shape
2. Decrease the isovalue to 0.01 -- observe how the surface changes
3. Increase the isovalue to 0.08 -- observe again
4. Return to 0.03

**Explain:**

**Q3.7:** The isosurface shows all points in space where |psi(x,y,z)| = isovalue. As you adjusted the slider:
- What happened when you decreased the isovalue to 0.01?
- What happened when you increased it to 0.08?
- Does the orbital have a definite, physical edge? Or is the boundary you see simply a chosen threshold? Explain.

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

> **Misconception check:** A student says, "The orbital is the colored surface I see on screen. Electrons live inside that surface and cannot exist outside it." What is wrong with this statement?

### Step 3.5: Orbital Classification Summary

Return to H2O and examine several orbitals systematically.

**Examine all H2O orbitals:** Using the same H2O (STO-3G) calculation, step through each orbital in the orbital selector dropdown: MO 1, MO 2, MO 3, MO 4, and MO 5 (HOMO). Examine MOs 1 through 5 (all occupied orbitals) and classify each one.

**Q3.8:** Complete the orbital classification table for H2O (STO-3G):

| MO Index | Approximate Description | Bonding / Nonbonding / Antibonding | Key Visual Feature |
|----------|------------------------|-------------------------------------|-------------------|
| 1 | | | |
| 2 | | | |
| 3 | | | |
| 4 | | | |
| 5 (HOMO) | | | |

---

### Checkpoint: Section 3 Deliverables

Before moving on:

1. Export an orbital visualization artifact: click the **Export** button (in the Module E header bar) and save as `orbital-artifact.json`
2. Take screenshots of: (a) H2 sigma_g orbital, (b) H2 sigma_u* orbital, (c) H2O HOMO at isovalue 0.03

- [ ] Answers to Q3.1 through Q3.8
- [ ] Screenshots of H2 bonding and antibonding orbitals
- [ ] Screenshot of H2O HOMO
- [ ] Exported `orbital-artifact.json`

---

## Section 4: Synthesis and Reflection (~5 min)

Now connect what you have learned across all three sections.

### Connecting 3D Structure, Energy, and Orbitals

The three topics in this lab are deeply connected:

1. **3D molecular geometry** determines the distances between atoms
2. **The PES** shows how the total energy depends on geometry -- the minimum tells us the equilibrium structure
3. **Molecular orbitals** describe how electrons distribute themselves across the atoms for a given geometry

At the equilibrium geometry (the PES minimum), the electrons have found the arrangement that minimizes the total energy. The shapes of the occupied orbitals at that geometry are the result.

---

**Q4.1:** In Section 2, you found that the H2 PES has a minimum near R = 1.4 bohr. In Section 3, you saw that the bonding orbital (sigma_g) has high electron density between the nuclei.

Connect these observations: how does the shape of the bonding orbital explain why the energy minimum exists? What would happen to the orbital shape -- and to the stabilization energy -- if R were much larger than the equilibrium value?

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

**Q4.2:** Reflect on the three types of representation you used today:

| Representation | What It Shows | Example |
|----------------|---------------|---------|
| 3D molecular viewer | Atomic positions and bonds | H2O bent geometry |
| PES curve | Energy vs. geometry parameter | H2 bond-length scan |
| Orbital isosurface | Electron distribution for one MO | H2O HOMO |

A student who has only seen the PES curve might think of a molecule as two balls on a spring. A student who has only seen the orbital might think of a molecule as a colored cloud. How does combining all three representations give you a more complete picture of a molecule than any single representation alone?

*Your answer:* _______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

---

## Final Deliverables Checklist

Before submitting, verify that you have completed:

**Section 1: 3D Molecular Exploration (8 questions)**
- [ ] **Q1.1:** Predicted H2O geometry and bond angle
- [ ] **Q1.2:** Compared 3D view to prediction
- [ ] **Q1.3:** Predicted largest and smallest overlap elements
- [ ] **Q1.4:** Recorded overlap matrix values
- [ ] **Q1.5:** Explained diagonal elements of overlap matrix
- [ ] **Q1.6:** Predicted H2 off-diagonal overlap
- [ ] **Q1.7:** Compared H2 and H2O hydrogen-hydrogen overlaps
- [ ] **Q1.8:** Predicted effect of bond stretching on overlap

**Section 2: Potential Energy Surface Scanning (8 questions)**
- [ ] **Q2.1:** Sketched predicted PES curve for H2
- [ ] **Q2.2:** Compared computed PES to prediction
- [ ] **Q2.3:** Explained origin of energy minimum
- [ ] **Q2.4:** Identified cause of repulsive wall
- [ ] **Q2.5:** Compared large-R energy to isolated atom limit
- [ ] **Q2.6:** Explained RHF dissociation failure
- [ ] **Q2.7:** Compared STO-3G and 6-31G energies (variational principle)
- [ ] **Q2.8:** Evaluated "bigger basis is always better" claim

**Section 3: Orbital Visualization (8 questions)**
- [ ] **Q3.1:** Predicted shape of MO 1 for H2O
- [ ] **Q3.2:** Described core orbital and compared to prediction
- [ ] **Q3.3:** Classified H2O HOMO character
- [ ] **Q3.4:** Predicted antibonding orbital shape
- [ ] **Q3.5:** Completed sigma_g vs. sigma_u* comparison table
- [ ] **Q3.6:** Predicted isovalue slider effect
- [ ] **Q3.7:** Explained isovalue meaning and orbital boundaries
- [ ] **Q3.8:** Completed H2O orbital classification table

**Section 4: Synthesis (2 questions)**
- [ ] **Q4.1:** Connected bonding orbital shape to PES minimum
- [ ] **Q4.2:** Reflected on combining multiple representations

**Artifacts and Screenshots:**
- [ ] Screenshot: H2O 3D view with atom labels
- [ ] Screenshot: H2 PES curve with equilibrium marker
- [ ] Screenshot: H2 sigma_g and sigma_u* orbitals
- [ ] Screenshot: H2O HOMO isosurface
- [ ] Exported `pes-artifact.json`
- [ ] Exported `orbital-artifact.json`

---

## Appendix A: Troubleshooting

**Problem:** The 3D viewer shows "WebGL required" instead of a molecule
**Solution:** Try a different browser (Chrome or Firefox recommended). Ensure hardware acceleration is enabled in your browser settings.

**Problem:** I cannot find the controls described in a step
**Solution:** Make sure you are on the correct workflow tab (Single Point, Optimize, PES Scan, or Compare). The tab bar is at the top of the Module E page.

**Problem:** The PES scan appears to hang or takes very long
**Solution:** Try reducing the number of scan points (e.g., from 20 to 10). Larger molecules and bigger basis sets take longer per point. For this lab, H2 with STO-3G should complete a 20-point scan quickly.

**Problem:** Orbital isosurfaces look blocky or pixelated
**Solution:** This is normal at the default grid resolution. The orbital extends smoothly in reality -- the grid is an approximation.

**Problem:** The orbital isosurface disappears when I increase the isovalue
**Solution:** At high isovalues, the isosurface encloses only a very small region near the nucleus. Try decreasing the isovalue back toward the default (0.03).

**Problem:** I see two different colors or opacities on the orbital lobes
**Solution:** This is intentional. Positive lobes (where psi > 0) are rendered as solid surfaces, and negative lobes (where psi < 0) are rendered as translucent. This distinguishes the sign of the wavefunction without relying on color alone.

---

## Appendix B: Glossary

**Antibonding orbital:** A molecular orbital with a node between bonded atoms, resulting in decreased electron density in the bonding region. Occupying an antibonding orbital destabilizes the molecule.

**Bonding orbital:** A molecular orbital with increased electron density between bonded atoms, which stabilizes the molecule and contributes to bond formation.

**CPK coloring:** A standard color scheme for atoms in molecular visualization (Corey-Pauling-Koltun). Hydrogen is white, carbon is gray, nitrogen is blue, oxygen is red.

**Dissociation limit:** The energy of a molecule as the bond length approaches infinity. For H2, the correct limit is the energy of two isolated hydrogen atoms.

**HOMO (Highest Occupied Molecular Orbital):** The highest-energy orbital that contains electrons in the ground state. Often determines chemical reactivity.

**Isovalue:** A threshold value used to define an isosurface. For orbital visualization, the isosurface shows all points where |psi(r)| equals the isovalue. The orbital extends beyond this surface.

**Isosurface:** A 3D surface connecting all points in space where a function has a constant value. Analogous to a contour line on a topographic map, but in three dimensions.

**LUMO (Lowest Unoccupied Molecular Orbital):** The lowest-energy orbital that does not contain electrons in the ground state.

**Nonbonding orbital:** A molecular orbital that is localized on a single atom or group and does not significantly contribute to bonding between atoms. Lone pairs are typically nonbonding.

**Overlap matrix (S):** A matrix whose elements S_ij measure the spatial overlap between basis functions i and j. Diagonal elements are 1.0 (normalized functions); off-diagonal elements range from 0 (no overlap) to approaching 1 (nearly identical functions).

**PES (Potential Energy Surface):** The energy of a molecule as a function of its geometry. For a diatomic, this is a curve of energy vs. bond length. For polyatomic molecules, PES scans can vary bond lengths, bond angles, or dihedral angles.

**RHF (Restricted Hartree-Fock):** A variant of Hartree-Fock theory where alpha and beta electrons share the same spatial orbitals. Works well near equilibrium but fails at the dissociation limit for many molecules.

**Variational principle:** The energy computed from any trial wavefunction is always greater than or equal to the true ground-state energy. Larger basis sets provide more variational freedom, yielding lower (more accurate) energies.

---

## Appendix C: References

1. Szabo, A. and Ostlund, N.S. (1996). *Modern Quantum Chemistry: Introduction to Advanced Electronic Structure Theory*. Dover Publications.

2. Johnstone, A.H. (1991). "Why is science difficult to learn? Things are seldom what they seem." *Journal of Computer Assisted Learning*, 7(2), 75-83.

3. Tasker, R. and Dalton, R. (2006). "Research into practice: Visualisation of the molecular world using animations." *Chemistry Education Research and Practice*, 7(2), 141-159.

4. Tsaparlis, G. and Papaphotis, G. (2009). "High-school students' conceptual difficulties and attempts at conceptual change: The case of basic quantum chemical concepts." *International Journal of Science Education*, 31(7), 895-930.

5. Wiggins, G. and McTighe, J. (2005). *Understanding by Design* (2nd ed.). ASCD.

---

**Lab Pack #2 Complete**

Thank you for completing this guided exploration of 3D molecular visualization, potential energy surfaces, and orbital visualization. The concepts you explored today -- connecting molecular geometry to energy and understanding what orbitals look like in space -- are central to how chemists think about molecules and chemical bonding.

---

*IQCP Lab Pack #2 v1.0 | Interactive Quantum Chemistry Playground | https://iqcp.dev*
