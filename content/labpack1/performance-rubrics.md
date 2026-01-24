# Lab Pack #1: Performance Task Rubrics

**Lab Pack:** 1 - From Boys to Orbitals
**Version:** 1.0
**Last Updated:** 2026-01-18
**Document Type:** Assessment Instruments (CONFIDENTIAL - Instructor Use Only)

---

## Overview

This document provides detailed analytic rubrics for three performance tasks that complement the concept check and worksheet questions. Performance tasks require students to apply their understanding to novel scenarios, demonstrating higher-order thinking and practical competency with the IQCP tool.

### Performance Task Summary

| Task | Module | Time | Points | Learning Outcomes |
|------|--------|------|--------|-------------------|
| **PT-Boys** | Boys Function | 10-15 min | 14 | LO1, LO4, LO6 |
| **PT-Rys** | Rys Quadrature | 10-15 min | 14 | LO2, LO4, LO5 |
| **PT-SCF** | SCF Sandbox | 10-15 min | 14 | LO3, LO4, LO5, LO6 |
| **Total** | - | 30-45 min | **42** | - |

### Rubric Scale (5-Level Analytic)

All rubrics use the following consistent scale:

| Score | Level | General Description |
|-------|-------|---------------------|
| **4** | Exemplary | Complete, accurate, insightful; exceeds expectations |
| **3** | Proficient | Mostly complete with minor errors; meets expectations |
| **2** | Developing | Partial understanding; significant gaps or errors |
| **1** | Beginning | Minimal evidence of understanding; major errors |
| **0** | No Response | No attempt or completely incorrect |

---

## PT-Boys: Boys Function Regime Prediction

### Task Description

**Given:** A specific (m, T) parameter combination not explored in the worksheet.

**Task:**
1. Predict which computational method (series or recurrence) will be used, based on the m-dependent turnover point
2. Verify your prediction using IQCP
3. Explain any discrepancy between prediction and observation

**Evidence Required:**
- Written prediction with reasoning (before using IQCP)
- Screenshot showing the method used (from Internals mode)
- Written explanation comparing prediction to actual result

**Time Allocation:** 10-15 minutes

### Example Deep Links

**Scenario A (T = 9.5, m = 3):**
> **[Open: Scenario A Parameters](http://localhost:5173/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgGZWAFSIBOIQFZWKLDADu-bHgDmCkAF9WpLDskzFxFh4dOh4CHIa7u5AA)**

**Scenario B (T = 28.0, m = 0):**
> **[Open: Scenario B Parameters](http://localhost:5173/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgNYAVIgCYAHKxRYYAd37Y8AcwUgAvq1JYdkmYuJY8OnQ8BDkNV1cgA)**

**Scenario C (T = 14.0, m = 8):**
> **[Open: Scenario C Parameters](http://localhost:5173/boys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj8ARjgCeYCSvVFQcIgA5WAFSIMALKxRYYAd37Y8AcwUgAvq1JYdkmYuJY8OnQ8BDkNV1cgA)**

### Task Prompt (Student Version)

---

**PT-Boys: Method Prediction Task**

Your instructor will provide a (T, m) combination. Complete the following:

**Part A: Prediction (Before opening IQCP)**

Based on what you learned about Boys function computational methods, predict which method will be used:

- [ ] Series expansion
- [ ] Recurrence relation (erf + upward recurrence)

**Explain your reasoning** (2-3 sentences): Why did you choose this method? What is the approximate turnover point for your assigned m value?

_______________________________________________

_______________________________________________

_______________________________________________

**Part B: Theoretical Regime Identification**

The lecture notes describe **three theoretical regimes**:
- Small T (T < 25): Series expansion
- Moderate T (25 <= T < 30+5m): erf + upward recurrence
- Large T (T >= 30+5m): Asymptotic expansion

For your assigned (T, m), which theoretical regime does this fall into?

- [ ] Small T (Series)
- [ ] Moderate T (erf + recurrence)
- [ ] Large T (Asymptotic)

Calculate the threshold: 30 + 5m = 30 + 5*___ = ___

Explain your classification: _______________________________________________

**Part C: Verification (Using IQCP)**

1. Open IQCP to the Boys module with your assigned parameters
2. Switch to Internals mode
3. Capture a screenshot showing the method used
4. Record the computed F_m(T) value: _______________

**Part D: Explanation**

Did your IQCP method prediction match the actual method?

- [ ] Yes, my prediction was correct
- [ ] No, the actual method was different

Does the IQCP method match the theoretical regime you identified? (Note: IQCP may use "Recurrence" for both Moderate and Large T theoretical regimes)

_______________________________________________

**If your prediction was correct:** Explain what specific knowledge from the worksheet helped you make this prediction. How did knowing the m-dependent turnover help?

**If your prediction was incorrect:** Explain why you think the actual method was used instead. What would you need to change about your understanding of the turnover points?

_______________________________________________

_______________________________________________

_______________________________________________

---

### Analytic Rubric: PT-Boys (14 points total)

#### Dimension 1: Prediction Accuracy (4 points)

| Score | Criteria |
|-------|----------|
| **4** | Correct method prediction with explicit reference to m-dependent turnover points |
| **3** | Correct method prediction with general reasoning about T and m relationship |
| **2** | Incorrect prediction but demonstrates knowledge of both methods (series/recurrence) |
| **1** | Incorrect prediction with minimal reasoning |
| **0** | No prediction or completely incorrect reasoning |

**Scoring Notes:**
- IQCP uses TWO methods: Series (T < turnover(m)) and Recurrence (T >= turnover(m))
- Turnover points: m=0-1: 0, m=2: 0.87, m=5: 2.11, m=10: 4.05, m=20: 7.84, m=30: 11.58
- Accept predictions that are correct for the given (m, T) scenario
- Award partial credit for recognizing that method selection depends on BOTH T and m
- **Note:** There is NO asymptotic regime in IQCP

**Bonus Challenge (Theory vs. Implementation):**
For advanced students, ask them to also identify which theoretical regime the parameters fall into:
- Small T (T < 25): Series expansion
- Moderate T (25 <= T < 30+5m): erf + upward recurrence
- Large T (T >= 30+5m): Asymptotic expansion

Then discuss why IQCP combines the moderate and large T theoretical regimes into a single "Recurrence" method.

#### Dimension 2: Verification Quality (4 points)

| Score | Criteria |
|-------|----------|
| **4** | Clear screenshot showing Internals mode with method visible; correct F_m(T) value recorded |
| **3** | Screenshot present with method identifiable; F_m(T) value recorded |
| **2** | Screenshot present but method partially visible or value missing |
| **1** | Screenshot provided but wrong mode or parameters |
| **0** | No screenshot or completely unusable evidence |

**Scoring Notes:**
- Screenshot should clearly show the "Method" indicator in Internals mode (Series or Recurrence)
- F_m(T) value should be recorded to at least 4 significant figures
- Accept equivalent numerical notation (e.g., 2.34e-5 = 0.0000234)

#### Dimension 3: Reasoning/Explanation (6 points)

| Score | Criteria |
|-------|----------|
| **6** | Insightful explanation connecting m-dependent turnover to method selection; demonstrates understanding of why turnover varies with m |
| **5** | Clear explanation with correct connections to turnover points; minor detail missing |
| **4** | Adequate explanation demonstrating understanding of the two-method system |
| **3** | Partial explanation with some correct elements but gaps in reasoning |
| **2** | Limited explanation; shows awareness of concepts but weak connections |
| **1** | Minimal explanation; mostly incorrect or irrelevant |
| **0** | No explanation or completely incorrect |

**For Correct Predictions:**
- Full credit requires explaining WHY the method is appropriate based on the m-dependent turnover (e.g., "For m=5, turnover is about 2.1, and T=3.0 is above this, so recurrence is used")
- Partial credit for stating the prediction was correct without substantial reasoning

**For Incorrect Predictions:**
- Full credit still available for thoughtful analysis of why the actual method was chosen
- Must demonstrate learning from the discrepancy (e.g., "I assumed fixed boundaries, but IQCP uses m-dependent turnover")
- Partial credit for acknowledging error without substantive analysis

---

## PT-Rys: Quadrature Order Selection

### Task Description

**Given:** A target accuracy requirement (e.g., 1e-8) and either a T value OR a shell quartet type.

**Task:**
1. Determine the minimum quadrature order needed to achieve the target accuracy
2. Use IQCP's error curve or shell quartet selector to justify your selection
3. For shell quartet scenarios: verify using the root count rule n_r = floor(L/2) + 1
4. Export a run artifact documenting your selection

**Evidence Required:**
- Screenshot of the error curve or shell quartet selector showing your order selection
- Written justification for the chosen order (including root count verification if applicable)
- Exported run artifact file

**Time Allocation:** 10-15 minutes

### Example Deep Links

**Scenario A (T = 8.0, target = 1e-6):**
> **[Open: Scenario A Parameters](http://localhost:5173/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEAzKwAqRABysa6AOYwq-BjAC0ANhABfVqSxFQU6YuIwAD3Y5BCwCFxcgA)**

**Scenario B (T = 18.0, target = 1e-8):**
> **[Open: Scenario B Parameters](http://localhost:5173/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEArKwAqRBgA5WNdAHMYVfgxgBaEyAC+rUliKgp0xcRgAHuxyCFgELi5AA)**

**Scenario C (T = 5.0, target = 1e-8):**
> **[Open: Scenario C Parameters](http://localhost:5173/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEArKwAqRfSBroA5jCr8GMALQAOEAF9WpLEVBTpi4jAAe7HIIWATOzkA)**

**Scenario D (Shell Quartet: dd|dd, verify root count):**
> **[Open: Scenario D - Shell Quartet](http://localhost:5173/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEGrACoMAzABUiAJhD6Q+oQHMYVfgxgBaAGwgAvq1JYioKdMXEseOnQ8BDk1Z2cgA)**

**Scenario E (Shell Quartet: ff|pp, verify root count):**
> **[Open: Scenario E - Shell Quartet](http://localhost:5173/rys?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj90ATzASVawqAKEGrACoMAzABUiAJhC6j+kQHMYVfgxgBaAGwgAvq1JYioKdMXEseOnQ8BDk1F2cgA)**

### Task Prompt (Student Version)

---

**PT-Rys: Quadrature Order Selection Task**

Your instructor will provide a T value and target accuracy. Complete the following:

**Given:**
- T = _______ (assigned by instructor)
- Target accuracy = _______ (assigned by instructor)

**Part A: Order Determination**

Using IQCP's Rys module:

1. Set the T value as specified
2. Examine the error curve and/or recommended order indicator
3. Identify the **minimum** quadrature order that achieves the target accuracy

**Selected order: n = _______**

**Part B: Evidence**

Capture a screenshot of the error curve that shows:
- Your selected order marked
- The error at your selected order
- Evidence that lower orders do not meet the target

*(Attach screenshot)*

**Part C: Justification**

Explain your order selection in 3-4 sentences:
- Why is this order sufficient?
- Why would a lower order be insufficient?
- Is there a reason not to use a much higher order?

_______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

**Part D: Artifact Export**

1. Export your run artifact
2. Save as `rys-performance-task.json`
3. Submit with your responses

- [ ] I have exported and will submit my artifact

---

### Alternative Task Prompt: Shell Quartet Selection (Student Version)

---

**PT-Rys: Shell Quartet Order Selection Task**

Your instructor will assign a shell quartet type. Complete the following:

**Given:**
- Shell quartet = _______ (e.g., (dd|pp), (ff|ss), etc.)

**Part A: Root Count Calculation**

Using the formula n_r = floor(L/2) + 1 where L = l_A + l_B + l_C + l_D:

1. Identify the angular momentum of each shell:
   - l_A = ___ (s=0, p=1, d=2, f=3)
   - l_B = ___
   - l_C = ___
   - l_D = ___

2. Calculate total angular momentum: L = ___ + ___ + ___ + ___ = ___

3. Apply the root count formula: n_r = floor(___/2) + 1 = ___

**Part B: IQCP Verification**

1. Open the Rys module and use the shell quartet selector
2. Set your assigned shell quartet
3. Record the order IQCP selects: n = ___
4. Does IQCP's selection match your calculation? ___

*(Attach screenshot)*

**Part C: Justification**

Explain in 2-3 sentences why this number of quadrature roots is needed for this shell quartet. What would happen if fewer roots were used?

_______________________________________________

_______________________________________________

_______________________________________________

**Part D: Artifact Export**

1. Export your run artifact
2. Save as `rys-shell-quartet-task.json`
3. Submit with your responses

- [ ] I have exported and will submit my artifact

---

### Analytic Rubric: PT-Rys (14 points total)

#### Dimension 1: Order Selection Correctness (4 points)

| Score | Criteria |
|-------|----------|
| **4** | Selected order is exactly the minimum that meets target; demonstrates efficiency consideration |
| **3** | Selected order meets target accuracy (may be 1 higher than minimum) |
| **2** | Selected order is within 2 of correct; shows understanding of relationship |
| **1** | Selected order is too low (fails accuracy) or too high (>2 above minimum) |
| **0** | No selection or completely inappropriate order |

**Scoring Notes:**
- "Minimum order" means the smallest n where max reconstruction error < target
- Accept n or n+1 for full credit if student explains conservative choice
- Common acceptable ranges (T-dependent, approximate):
  - T=8, target 1e-6: n = 5-6
  - T=18, target 1e-8: n = 7-8
  - T=5, target 1e-10: n = 9-10

#### Dimension 2: Evidence Quality (Screenshot) (4 points)

| Score | Criteria |
|-------|----------|
| **4** | Screenshot clearly shows error curve with selected order marked; error values visible; T value confirmed |
| **3** | Screenshot shows error information with selected order identifiable |
| **2** | Screenshot present but error values or order marking unclear |
| **1** | Screenshot provided but minimal relevant information visible |
| **0** | No screenshot or completely unusable |

**Scoring Notes:**
- Error curve or error table should be visible
- Selected order should be highlighted or annotated
- T value should be verifiable from screenshot

#### Dimension 3: Justification Quality (4 points)

| Score | Criteria |
|-------|----------|
| **4** | Clear explanation addressing all three points: sufficiency, lower order insufficiency, and efficiency consideration |
| **3** | Addresses two of three points clearly |
| **2** | Addresses one point clearly or two points partially |
| **1** | Minimal justification with limited reasoning |
| **0** | No justification or completely incorrect |

**Key Points for Full Credit:**
- **Sufficiency:** "Order n achieves error of [value], which is below the target of [target]"
- **Lower orders:** "Order n-1 has error [value], which exceeds the target"
- **Efficiency:** "Higher orders would work but cost more computation without benefit"

#### Dimension 4: Artifact Quality (2 points)

| Score | Criteria |
|-------|----------|
| **2** | Valid JSON artifact with correct module, T value, and order; results section populated |
| **1** | Valid artifact with minor discrepancy (e.g., slightly different T or order from reported values) |
| **0** | No artifact, invalid JSON, or major discrepancies |

**Scoring Notes:**
- Open artifact file and verify:
  - `module: "rys"` or equivalent identifier
  - T value matches task parameters
  - Order matches reported selection
  - Results/computed values present

---

## PT-SCF: Convergence Analysis

### Task Description

**Given:** A molecular system preset (e.g., H2O or LiH).

**Task:**
1. Run SCF with DIIS disabled; record iterations to convergence
2. Run SCF with DIIS enabled; record iterations to convergence
3. Compare convergence behavior and export both artifacts

**Evidence Required:**
- Two run artifacts (one without DIIS, one with DIIS)
- Comparison table of iteration counts and final energies
- Written explanation of observed differences

**Time Allocation:** 10-15 minutes

### Example Deep Links

**Scenario A: LiH System**

Without DIIS:
> **[Open: LiH without DIIS](http://localhost:5173/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn85WKJSo4AzAHMJiAB7666IgFYBrCPhT84MaVlI4CX8sMCJlBDkwGABfVlIsDUkZRWIYO3Y5BCwCGJigA)**

With DIIS:
> **[Open: LiH with DIIS](http://localhost:5173/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn85WKJSo4AzAHMJiAB7666IgFYBrCPhT84MaVlI4CX8sMCIqMhgAX1ZSLA1JGUViGDt2OQQsAiiooA)**

**Scenario B: NH3 System**

Without DIIS:
> **[Open: NH3 without DIIS](http://localhost:5173/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn88UAMyUqOEwHMJiAB7666IgFYBrCPhT84MaVlI4CX8sMCJlBDkwGABfVlIsDUkZRWIYO3Y5BCwCGJigA)**

With DIIS:
> **[Open: NH3 with DIIS](http://localhost:5173/scf?run=N4IgzgxgFgpgtgQwPoDcYCcwEsD2A7EALhHQFc8kwAXBKmVARhABoQEAHd1DbfIkAAwA6BkIEsQcHABNSAGxj9IAMwkqioMAE9q8JFmn88UAMyUqOEwHMJiAB7666IgFYBrCPhT84MaVlI4CX8sMCIqMhgAX1ZSLA1JGUViGDt2OQQsAiiooA)**

### Task Prompt (Student Version)

---

**PT-SCF: Convergence Analysis Task**

Your instructor will assign a molecular system. Complete the following:

**Assigned System:** _____________ (e.g., LiH, NH3, H2O)

**Part A: SCF Without DIIS**

1. Open IQCP to the SCF module with your assigned system
2. Set convergence to "medium" (or as specified)
3. Set DIIS to **OFF**
4. Run the calculation
5. Record results:

| Metric | Value |
|--------|-------|
| Iterations to convergence | |
| Final energy (Hartree) | |
| Converged? (Yes/No) | |

6. Export artifact as `scf-no-diis.json`

**Part B: SCF With DIIS**

1. Keep the same system and convergence settings
2. Set DIIS to **ON**
3. Run the calculation
4. Record results:

| Metric | Value |
|--------|-------|
| Iterations to convergence | |
| Final energy (Hartree) | |
| Converged? (Yes/No) | |

5. Export artifact as `scf-with-diis.json`

**Part C: Comparison Table**

| Metric | Without DIIS | With DIIS | Difference |
|--------|--------------|-----------|------------|
| Iterations | | | |
| Final Energy (Ha) | | | |
| Energy Agreement? | - | - | Yes/No |

**Part D: Explanation**

In 4-5 sentences, explain:
1. How did DIIS affect the number of iterations?
2. Did both runs converge to the same energy (within tolerance)?
3. Based on what you learned about DIIS, why does it achieve faster convergence?
4. Under what circumstances might you choose NOT to use DIIS?

_______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

_______________________________________________

**Part E: Artifact Submission**

- [ ] I have exported `scf-no-diis.json`
- [ ] I have exported `scf-with-diis.json`
- [ ] Both artifacts will be submitted with this response

---

### Analytic Rubric: PT-SCF (14 points total)

#### Dimension 1: Data Collection Accuracy (4 points)

| Score | Criteria |
|-------|----------|
| **4** | All values recorded accurately for both runs; iteration counts and energies match artifacts |
| **3** | Minor discrepancy in one value (e.g., rounding) but overall accurate |
| **2** | One significant error or missing value |
| **1** | Multiple errors or incomplete data |
| **0** | No data collected or completely incorrect |

**Scoring Notes:**
- Energy values should agree to at least 10^-6 Ha between the two runs (same converged result)
- DIIS run should have fewer iterations (typically 30-50% reduction)
- If calculation did not converge, student should note this and explain

**Expected Ranges (approximate, convergence-dependent):**

| System | Without DIIS | With DIIS | Energy (Ha) |
|--------|--------------|-----------|-------------|
| H2 | 8-15 | 4-8 | -1.1167 |
| LiH | 12-25 | 6-12 | -7.8634 |
| H2O | 15-30 | 8-15 | -74.9630 |
| NH3 | 15-30 | 8-15 | -55.4546 |

#### Dimension 2: Comparison Quality (4 points)

| Score | Criteria |
|-------|----------|
| **4** | Table complete with correct differences; explicitly notes energy agreement |
| **3** | Table complete; differences calculated correctly |
| **2** | Table mostly complete; minor calculation error in differences |
| **1** | Table incomplete but shows attempt at comparison |
| **0** | No comparison table or completely incorrect |

**Scoring Notes:**
- "Difference" in iterations should be (without - with) to show DIIS reduction
- Energy agreement should be "Yes" if values agree to ~10^-8 Ha or better
- Accept "N/A" for difference if one run did not converge

#### Dimension 3: Explanation Quality (4 points)

| Score | Criteria |
|-------|----------|
| **4** | Addresses all four points clearly; demonstrates understanding of DIIS mechanism and appropriate applications |
| **3** | Addresses three of four points clearly |
| **2** | Addresses two points clearly or three points partially |
| **1** | Addresses one point or shows minimal understanding |
| **0** | No explanation or completely incorrect |

**Key Points for Full Credit:**

1. **DIIS effect on iterations:** Should note significant reduction (often 30-50%)
2. **Energy agreement:** Should confirm same final energy (DIIS is an acceleration, not approximation)
3. **Why DIIS works:** Should mention extrapolation from history, error minimization, or avoiding oscillations
4. **When NOT to use DIIS:** Acceptable answers include:
   - Very far from solution (initial iterations)
   - Seeking multiple SCF solutions
   - Near-degenerate systems where DIIS may bias
   - Simple systems where overhead not worthwhile

#### Dimension 4: Artifact Quality (2 points)

| Score | Criteria |
|-------|----------|
| **2** | Both artifacts valid; DIIS settings match reported (on/off); system and energy match |
| **1** | One artifact valid and correct; other has minor issues |
| **0** | Missing artifacts, invalid JSON, or major discrepancies |

**Verification Checklist:**
- [ ] `scf-no-diis.json` has DIIS=false and reported iteration count
- [ ] `scf-with-diis.json` has DIIS=true and reported iteration count
- [ ] Both artifacts show same system preset
- [ ] Energies in artifacts match reported values

---

## Scoring Procedures

### Administering Performance Tasks

1. **Timing:** Allocate 10-15 minutes per task; total 30-45 minutes for all three
2. **Order:** Tasks may be completed in any order
3. **Resources:** Students should have access to IQCP and this prompt document
4. **Submission:** Collect written responses and artifact files

### Scoring Workflow

1. **Artifact verification first:** Open each artifact file and verify basic validity
2. **Evidence alignment:** Check that screenshots/artifacts match reported values
3. **Score each dimension independently** using the rubric criteria
4. **Sum dimension scores** for task total
5. **Record any notes** for borderline cases or unusual responses

### Point Allocation Summary

| Task | D1: Core Skill | D2: Evidence | D3: Reasoning | D4: Artifact | Total |
|------|---------------|--------------|---------------|--------------|-------|
| PT-Boys | 4 (Prediction) | 4 (Screenshot) | 6 (Explanation) | - | 14 |
| PT-Rys | 4 (Order) | 4 (Screenshot) | 4 (Justification) | 2 (Export) | 14 |
| PT-SCF | 4 (Data) | 4 (Comparison) | 4 (Explanation) | 2 (Artifacts) | 14 |
| **Total** | **12** | **12** | **14** | **4** | **42** |

### Converting to Course Grade

| Raw Score | Percentage | Suggested Grade |
|-----------|------------|-----------------|
| 38-42 | 90-100% | A |
| 34-37 | 80-89% | B |
| 29-33 | 70-79% | C |
| 25-28 | 60-69% | D |
| 0-24 | <60% | F |

---

## Inter-Rater Reliability Guidelines

### Calibration Procedure

Before grading student submissions:

1. **Training set:** Have all graders independently score 3-5 sample responses
2. **Comparison meeting:** Compare scores and discuss discrepancies
3. **Norming:** Establish consensus on borderline cases
4. **Reference responses:** Create anchor responses at each score level

### Consistency Guidelines

| Situation | Recommended Action |
|-----------|-------------------|
| Score differs by 1 point | Accept higher score |
| Score differs by 2+ points | Discuss and reach consensus |
| Borderline between levels | Read criteria again; award higher if student demonstrates understanding |
| Missing artifact but correct reasoning | Maximum 2-point deduction (artifact dimensions only) |
| Incorrect artifact but correct written response | Investigate; may indicate copying |
| Screenshot unclear but values correct | Award full evidence points |

### Documentation Requirements

For each graded submission, record:
- Total score by task
- Any dimension scores that were borderline (note reasoning)
- Any unusual responses for potential FAQ updates

### Reliability Statistics (Post-Grading)

After grading is complete, calculate:
- Cohen's kappa for each dimension (target: kappa > 0.7)
- Percentage exact agreement (target: > 80%)
- Range of scores per dimension

---

## Sample Student Responses

### PT-Boys: Exemplary Response (14/14)

**Given:** T = 14.0, m = 8

**Part A - Prediction:**
> "I predict the **recurrence method** will be used. Based on the worksheet exploration, IQCP uses m-dependent turnover points. For m=8, looking at the turnover pattern (m=7 has turnover ~2.89, m=8 has turnover ~3.28), my T = 14.0 is well above this threshold. Since T > turnover(m=8), the recurrence method using erf(sqrt(T)) and upward recurrence should be employed."

**Part B - Verification:**
*(Screenshot clearly showing Internals mode with "Method: Recurrence" displayed, T=14.0, m=8)*

F_8(14.0) = 7.234e-08

**Part C - Explanation:**
> "My prediction was correct. IQCP uses a two-method system where the switch from series to recurrence depends on the order m. For m=8, the turnover point is approximately 3.28. Since my T = 14.0 is well above this turnover, the recurrence method was used. This makes sense because the series expansion would require many terms to converge for T values significantly above the turnover, whereas the recurrence method efficiently computes F_0 using erf(sqrt(T)) and then applies upward recurrence."

**Scoring:**
- D1 (Prediction): 4 - Correct method with explicit reference to m-dependent turnover
- D2 (Evidence): 4 - Clear screenshot, value recorded
- D3 (Reasoning): 6 - Insightful explanation connecting m-dependent turnover to method efficiency

---

### PT-Rys: Proficient Response (11/14)

**Given:** T = 18.0, target = 1e-8

**Part A - Order Selection:**
> n = 8

**Part B - Screenshot:**
*(Screenshot shows error curve, but order not explicitly marked; error values visible)*

**Part C - Justification:**
> "Order 8 achieves an error of about 2.5e-9, which is below the 1e-8 target. Order 7 has an error around 5e-8, which is too high. I chose 8 instead of a higher order to minimize computation."

**Part D - Artifact:**
*(Valid JSON with T=18.0, n=8)*

**Scoring:**
- D1 (Order): 4 - Correct minimum order
- D2 (Evidence): 3 - Screenshot present but order not marked
- D3 (Justification): 3 - Addresses sufficiency and lower order; efficiency mentioned briefly
- D4 (Artifact): 1 - Valid but T slightly different (17.9 vs 18.0)

---

### PT-SCF: Developing Response (8/14)

**Given:** H2O system

**Part A - Without DIIS:**
| Metric | Value |
|--------|-------|
| Iterations | 18 |
| Final energy | -74.96 Ha |
| Converged? | Yes |

**Part B - With DIIS:**
| Metric | Value |
|--------|-------|
| Iterations | 9 |
| Final energy | -74.96 Ha |
| Converged? | Yes |

**Part C - Comparison Table:**
*(Table partially filled; difference column empty)*

**Part D - Explanation:**
> "DIIS reduced the iterations from 18 to 9. Both runs got the same energy. DIIS works by using information from previous steps."

**Part E - Artifacts:**
*(Only one artifact submitted - with DIIS)*

**Scoring:**
- D1 (Data): 3 - Energy reported with insufficient precision (-74.96 vs -74.9630)
- D2 (Comparison): 2 - Table incomplete (no differences calculated)
- D3 (Explanation): 2 - Two points addressed but superficially
- D4 (Artifacts): 1 - Only one artifact submitted

---

## Appendix: Learning Outcome Alignment

### Performance Task to Learning Outcome Mapping

| Task | Primary LOs | Assessment Focus |
|------|-------------|------------------|
| PT-Boys | LO1, LO6 | Method prediction requires understanding Boys function behavior and m-dependent turnover (LO1) and applying reasoning to predict outcomes (LO6) |
| PT-Boys | LO4 | Connecting T and m parameters to computational method selection (two-method system)
| PT-Rys | LO2, LO5 | Order selection requires understanding order-accuracy relationship (LO2) and documenting with artifacts (LO5) |
| PT-Rys | LO4 | Connecting quadrature order and T to computational outcomes |
| PT-SCF | LO3, LO6 | Convergence analysis requires interpreting plots (LO3) and predicting DIIS effects (LO6) |
| PT-SCF | LO4, LO5 | Connecting DIIS settings to outcomes (LO4) and documenting comparisons (LO5) |

### Bloom's Taxonomy Level by Dimension

| Dimension | Cognitive Level | Justification |
|-----------|-----------------|---------------|
| Prediction/Selection | Apply | Using learned principles in new situations |
| Evidence/Screenshots | Remember/Understand | Demonstrating observation skills |
| Reasoning/Explanation | Analyze/Evaluate | Connecting concepts and justifying decisions |
| Artifact Export | Apply | Procedural application of tool features |

---

*Lab Pack #1 Performance Rubrics v1.0 | CONFIDENTIAL - Instructor Use Only*
*Interactive Quantum Chemistry Playground | https://iqcp.dev*
