# Lab Pack #1: Grading Rubric

**Lab Pack:** 1 - From Boys to Orbitals
**Version:** 1.0
**Last Updated:** 2026-01-18
**Document Type:** Instructor Materials (CONFIDENTIAL)

---

## Point Summary

| Section | Questions | Points | Per Question |
|---------|-----------|--------|--------------|
| Section 2: Boys Function | Q2.1-Q2.7 | 35 | 5 each |
| Section 3: Rys Quadrature | Q3.1-Q3.8 | 40 | 5 each |
| Section 4: SCF Convergence | Q4.1-Q4.7 | 35 | Variable |
| Section 5: Synthesis | Q5.1-Q5.3 | 15 | 5 each |
| **Total** | **25 questions** | **125** | |

**Note:** Point values may be scaled to 100 points by instructors using the formula: (raw score / 125) x 100.

### Section 4 Detailed Breakdown

| Question | Points | Content |
|----------|--------|---------|
| Q4.1 | 4 | Iterations (2) + Energy (2) |
| Q4.2 | 4 | Without DIIS (2) + With DIIS (2) |
| Q4.3 | 6 | Convergence pattern description |
| Q4.4 | 4 | Energy (2) + Convergence assessment (2) |
| Q4.5 | 5 | Symmetry observation (2) + Explanation (3) |
| Q4.6 | 6 | HOMO (2) + LUMO (2) + Gap significance (2) |
| Q4.7 | 6 | DIIS recommendation with rationale |
| **Total** | **35** | |

---

## General Grading Guidelines

### Numerical Answers

| Credit Level | Criteria | Percentage |
|--------------|----------|------------|
| Full | Within 5% or specified tolerance of expected value | 100% |
| Partial | Correct method, minor calculation error | 50% |
| Minimal | Correct units or order of magnitude only | 25% |
| None | Wrong answer, wrong units, or no attempt | 0% |

### Conceptual Answers

| Credit Level | Criteria | Percentage |
|--------------|----------|------------|
| Full | All key points addressed with clarity | 100% |
| Substantial | Most key points, minor omissions | 75% |
| Partial | Some understanding demonstrated, incomplete | 50% |
| Minimal | Limited understanding shown | 25% |
| None | Incorrect or no attempt | 0% |

### Common Considerations

1. **Units:** Accept Hartree, Ha, E_h, or atomic units interchangeably. Only deduct points if an incorrect unit conversion is attempted.

2. **Significant Figures:** Do not penalize for excess precision. Penalize only for grossly insufficient precision (e.g., reporting -1 Ha when -1.117 Ha expected).

3. **Terminology:** Accept reasonable synonyms (e.g., "Taylor series" for "series expansion," "eigenvalue decomposition" for "diagonalization").

4. **Partial Work:** Award partial credit for correct reasoning even if the final numerical answer is incorrect.

5. **Sign Errors:** Energy values must include the correct negative sign. Deduct 50% if sign is missing or wrong.

---

## Section 2: Boys Function (25 points)

### Q2.1: Limiting Value as T Approaches 0 (5 points)

**Question:** What value does F_0(T) approach as T approaches 0? Express this as a simple fraction in terms of m.

| Points | Criteria |
|--------|----------|
| **5** | States F_m(0) = 1/(2m+1) AND correctly identifies F_0(0) = 1 |
| **4** | Correct general formula OR correct F_0(0) = 1 with reasonable justification |
| **3** | Correct for m=0 only (F_0(0) = 1) but incorrect or missing general formula |
| **2** | Recognizes the limit exists and is a finite positive number |
| **1** | Demonstrates some understanding of the integral definition |
| **0** | Incorrect answer or no attempt |

**Key Concepts Required:**
- Recognition that when T=0, exp(-T*t^2) = 1
- Evaluation of simple power integral: integral of t^(2m) from 0 to 1 = 1/(2m+1)

**Common Errors to Watch:**
- Writing 1/m instead of 1/(2m+1)
- Confusing with exponential decay behavior at large T
- Stating "infinity" or "zero" as the limit

---

### Q2.2: T Value Where F_0(T) < 0.01 (5 points)

**Question:** At approximately what T value does F_0(T) become smaller than 0.01?

| Points | Criteria |
|--------|----------|
| **5** | Value in range 7000-9000 OR explicit recognition of very slow sqrt(1/T) decay |
| **4** | Value in range 5000-10000 with reasoning |
| **3** | Correct order of magnitude (thousands) |
| **2** | Recognizes decay is much slower than exponential; value > 100 |
| **1** | Any value > 50 with attempt at reasoning |
| **0** | Value < 50, completely wrong, or no attempt |

**Expected Answer:** T ~ 7854 (from solving sqrt(pi/4T) = 0.01)

**Acceptable Range:** T in range 7000-9000 for full credit

**Key Concepts Required:**
- Understanding that F_0(T) decays as sqrt(pi/4T) for large T
- Recognition that this is much slower than exponential decay

**Instructor Note:** This is intentionally a challenging question. Many students will not explore large enough T values. Award partial credit for demonstrating the exploration process and recognizing the unexpectedly slow decay.

---

### Q2.3: Computational Methods (5 points)

**Question:** What computational method is used for T = 0.5? For T = 15.0? For T = 35.0? (with m=0)

| Points | Criteria |
|--------|----------|
| **5** | All three methods correctly identified as Recurrence (for m=0), OR demonstrates understanding that method selection depends on both T and m |
| **4** | Two of three correct, with awareness that m affects method selection |
| **3** | Two correct, shows understanding of two-method system |
| **2** | One correct, or identifies both methods (series/recurrence) exist |
| **1** | Demonstrates awareness that different methods exist |
| **0** | None correct or no attempt |

**Expected Answers (for m=0):**
- T = 0.5: **Recurrence** (turnover=0 for m=0, so always recurrence)
- T = 15.0: **Recurrence**
- T = 35.0: **Recurrence**

**For higher m values (bonus exploration):**
- With m=5 and T=1.0: **Series** (since 1.0 < turnover(5)=2.11)
- With m=10 and T=3.0: **Series** (since 3.0 < turnover(10)=4.05)

**Key Concepts Required:**
- Ability to read IQCP "Internals" mode output
- Understanding that IQCP uses TWO methods (series and recurrence), NOT three
- Recognition that method selection depends on BOTH T and m via m-dependent turnover points

**Instructor Note:** Students who answer based on fixed T=12/T=30 boundaries from textbooks should receive partial credit with feedback that IQCP uses a more sophisticated m-dependent approach. There is NO asymptotic regime in IQCP.

---

### Q2.4: Comparison of F_0(10) vs F_5(10) (5 points)

**Question:** Compare F_0(10.0) and F_5(10.0). Which is larger? Why do higher-order Boys functions decay more rapidly?

| Points | Criteria |
|--------|----------|
| **5** | Correct comparison (F_0 > F_5) AND substantive explanation involving t^(2m) factor |
| **4** | Correct comparison AND general explanation about m-dependence |
| **3** | Correct comparison, weak or incomplete explanation |
| **2** | Correct comparison only, no explanation |
| **1** | Incorrect comparison but some relevant discussion |
| **0** | Incorrect comparison or no attempt |

**Expected Values:**
- F_0(10.0) = 0.280247 (larger)
- F_5(10.0) = 7.9 x 10^-5 (much smaller)

**Key Concepts Required:**
- Understanding that the t^(2m) factor suppresses contributions at larger t
- Recognizing that the integrand is increasingly concentrated near t=0 for higher m

---

### Q2.5: Why Different Computational Methods? (5 points)

**Question:** Why do quantum chemistry codes use different methods for different parameter values?

| Points | Criteria |
|--------|----------|
| **5** | Mentions BOTH accuracy/stability AND efficiency, with understanding that method selection depends on both T and m |
| **4** | Mentions two key concepts with some specifics |
| **3** | General idea about accuracy OR efficiency trade-offs |
| **2** | Mentions that "different methods work better for different ranges" |
| **1** | Acknowledges the existence of trade-offs |
| **0** | No substantive answer |

**Key Points for Full Credit:**
1. Series converges rapidly for small T, slowly for large T
2. Recurrence (erf + upward recurrence) is numerically stable for moderate-to-large T
3. The optimal switch point (turnover) depends on BOTH T and m
4. Using wrong method leads to loss of precision or slow computation

**Note:** IQCP uses only TWO methods (series and recurrence). Accept answers mentioning "asymptotic" if student demonstrates general understanding, but note that IQCP does not implement a separate asymptotic regime.

---

### Q2.6: IQCP Methods vs. Theoretical 3-Regime Model (5 points)

**Question:** What computational method does IQCP use for m=0 at T=35? For m=5 at T=45? What theoretical regime would these fall into?

| Points | Criteria |
|--------|----------|
| **5** | Both IQCP methods correctly identified as "Recurrence" AND both theoretical regimes correctly mapped (Large T for m=0,T=35; Moderate T for m=5,T=45) |
| **4** | Methods correct, one regime mapping slightly off |
| **3** | Methods correct but regime mapping incomplete or missing |
| **2** | One method/regime pair fully correct |
| **1** | Shows understanding that IQCP differs from theoretical description |
| **0** | No correct identification |

**Expected Answers:**
- m=0, T=35: IQCP method = **Recurrence**, Theoretical regime = **Large T (Asymptotic)**
- m=5, T=45: IQCP method = **Recurrence**, Theoretical regime = **Moderate T**

**Key Concepts Required:**
- Understanding that IQCP uses 2 methods while theory describes 3 regimes
- Recognition that theoretical regime thresholds are m-dependent (30+5m)

---

### Q2.7: Why Implementations Combine Moderate/Large T Methods (5 points)

**Question:** Why might implementations combine the moderate and large T methods rather than implementing all three separately?

| Points | Criteria |
|--------|----------|
| **5** | Mentions at least TWO valid reasons: efficiency, simplicity, numerical stability of recurrence, or practical accuracy |
| **4** | One good reason with substantive explanation |
| **3** | General statement about efficiency or simplicity with some reasoning |
| **2** | Acknowledges that implementations differ from theory |
| **1** | Minimal relevant content |
| **0** | No substantive answer |

**Key Points for Full Credit:**
1. Recurrence method works well for both moderate and large T
2. Implementation simplicity (fewer code paths)
3. Recurrence naturally transitions to asymptotic behavior as exp(-T) becomes negligible
4. Careful turnover tuning eliminates need for separate asymptotic method

---

## Section 3: Rys Quadrature (40 points)

### Q3.1: Roots and Weights Properties (5 points)

**Question:** Are all roots strictly between 0 and 1? Are all weights positive?

| Points | Criteria |
|--------|----------|
| **5** | Both correct ("Yes" to both) with evidence from IQCP |
| **4** | Both correct, no explicit verification mentioned |
| **3** | One correct, one missing or partially correct |
| **2** | One correct only |
| **1** | Demonstrates awareness of property requirements |
| **0** | Both incorrect or no attempt |

**Expected Answer:**
- Roots: Yes, all t_i are in (0, 1)
- Weights: Yes, all w_i > 0

**Key Concepts Required:**
- Understanding these are mathematical guarantees of proper Gaussian quadrature
- Ability to verify by inspecting IQCP output

---

### Q3.2: Effect of Quadrature Points (5 points)

**Question:** How does the number of quadrature points affect the computation?

| Points | Criteria |
|--------|----------|
| **5** | Discusses BOTH accuracy improvement AND cost increase with reasoning |
| **4** | Discusses both aspects, one more thoroughly than the other |
| **3** | Discusses only one aspect (accuracy OR cost) thoroughly |
| **2** | Mentions both aspects superficially |
| **1** | General statement about trade-offs |
| **0** | Does not address the question |

**Key Points for Full Credit:**
- More points = higher accuracy (exact integration of higher-degree polynomials)
- More points = higher computational cost (linear scaling)
- Trade-off: choose minimum n that meets accuracy requirements

---

### Q3.3: Reconstruction Errors at Different Orders (5 points)

**Question:** What is the approximate maximum reconstruction error for n=3? For n=5? For n=7?

| Points | Criteria |
|--------|----------|
| **5** | All three within one order of magnitude of expected |
| **4** | Two correct, one within factor of 10 |
| **3** | Correct trend (error decreases with n) with reasonable magnitudes |
| **2** | One value correct or correct trend identified |
| **1** | Demonstrates attempt to read IQCP error output |
| **0** | All incorrect or no attempt |

**Expected Answers (at T=10.0):**
- n=3: ~10^-4 to 10^-5
- n=5: ~10^-5 to 10^-6
- n=7: ~10^-6 to 10^-7

**Note:** Accept order-of-magnitude estimates. Exact values depend on IQCP implementation.

---

### Q3.4: Minimum Orders for Target Accuracies (5 points)

**Question:** At T = 10.0, what minimum quadrature order is needed for 1e-8 accuracy? For 1e-6 accuracy?

| Points | Criteria |
|--------|----------|
| **5** | Both answers in acceptable ranges |
| **4** | One in acceptable range, other close |
| **3** | Correct relative ordering (1e-8 needs more points than 1e-6) with reasonable values |
| **2** | One answer correct |
| **1** | Demonstrates understanding that tighter tolerance needs more points |
| **0** | Both incorrect or no attempt |

**Acceptable Ranges:**
- For 1e-8: n = 8, 9, 10, or 11
- For 1e-6: n = 5, 6, or 7

---

### Q3.5: Effect of T on Recommended Order (5 points)

**Question:** How does the recommended quadrature order change when T increases from 10 to 25? Why?

| Points | Criteria |
|--------|----------|
| **5** | Correct observation (order decreases) AND explanation involving moment decay |
| **4** | Correct observation AND general explanation about T-dependence |
| **3** | Correct observation, weak or missing explanation |
| **2** | Correct direction but incorrect reasoning |
| **1** | Incorrect direction with some relevant discussion |
| **0** | No substantive answer |

**Expected Answer:**
- Order **decreases** as T increases
- At higher T, Boys function moments decay faster
- Fewer quadrature points suffice for the same accuracy

---

### Q3.6: Shell Quartet (pp|pp) Root Count Verification (5 points)

**Question:** For (pp|pp), what is L? What order does IQCP select? Verify the formula n_r = floor(L/2) + 1.

| Points | Criteria |
|--------|----------|
| **5** | L=4 correct, n=3 correct, formula verification shown (floor(4/2)+1=3), match confirmed |
| **4** | Three of four parts correct |
| **3** | L and n correct, formula verification weak |
| **2** | L correct only, or n correct only |
| **1** | Shows attempt at calculation |
| **0** | Incorrect values or no attempt |

**Expected Answers:**
- L = 1+1+1+1 = **4** (each p shell has l=1)
- IQCP selected order = **3**
- Formula: floor(4/2) + 1 = 2 + 1 = **3**
- Match: **Yes**

---

### Q3.7: Shell Quartet (dd|pp) Root Count Verification (5 points)

**Question:** For (dd|pp), what is L? What order does IQCP select? Is this the minimum required?

| Points | Criteria |
|--------|----------|
| **5** | L=6 correct, n=4 correct, confirms this is the minimum required |
| **4** | L and n correct, minimum confirmation weak |
| **3** | Two of three parts correct |
| **2** | L correct only |
| **1** | Shows understanding of shell angular momentum |
| **0** | Incorrect L value or no attempt |

**Expected Answers:**
- L = 2+2+1+1 = **6** (d has l=2, p has l=1)
- IQCP selected order = **4**
- Minimum required = floor(6/2) + 1 = **4** (yes, this is the minimum)

---

### Q3.8: Algorithm 5.1 Moments and Hankel Matrix (5 points)

**Question:** For T=10, n=3: What are moments m_0, m_1, m_2? What is the Hankel matrix dimension?

| Points | Criteria |
|--------|----------|
| **5** | At least 2 moments within 10% AND correct dimension (3x3) |
| **4** | All moments within 20% OR correct dimension with 1 moment correct |
| **3** | Understands m_k = 2*F_k relationship, dimension correct |
| **2** | One moment approximately correct OR dimension correct |
| **1** | Shows attempt to read Algorithm 5.1 internals |
| **0** | No correct values |

**Expected Values:**
- m_0 = 2*F_0(10) = **0.5605** (or 5.6e-01)
- m_1 = 2*F_1(10) = **0.0280** (or 2.8e-02)
- m_2 = 2*F_2(10) = **0.0042** (or 4.2e-03)
- Hankel matrix dimension: **3 x 3** (n x n)

---

## Section 4: SCF Convergence (35 points)

### Q4.1: H2 Default Iterations and Energy (4 points)

**Question:** For H2 with medium convergence and DIIS enabled, how many iterations? What is the final energy?

| Points | Criteria |
|--------|----------|
| **4** | Iterations in range 3-10 AND energy within 10^-4 of -1.1167 Ha |
| **3** | One component correct |
| **2** | Energy within 10^-3 Ha OR iterations within 2 of expected |
| **1** | Correct order of magnitude for energy (around -1 Ha) |
| **0** | Both incorrect or no attempt |

**Expected Values:**
- Iterations: 4-8 (typically 5-6)
- Energy: -1.116716909173 Hartree

---

### Q4.2: H2 Iterations With/Without DIIS (4 points)

**Question:** For H2 with tight convergence, how many iterations without DIIS? With DIIS?

| Points | Criteria |
|--------|----------|
| **4** | Both in acceptable ranges AND DIIS case has fewer iterations |
| **3** | Correct relative comparison with reasonable values |
| **2** | One value correct |
| **1** | DIIS correctly identified as faster |
| **0** | No comparison or DIIS incorrectly identified as slower |

**Acceptable Ranges:**
- Without DIIS: 8-20 iterations
- With DIIS: 4-10 iterations

---

### Q4.3: Convergence Pattern Description (6 points)

**Question:** Describe the difference in convergence patterns. How does DIIS change the curve shape?

| Points | Criteria |
|--------|----------|
| **6** | Clear description of BOTH patterns AND explanation of DIIS effect |
| **5** | Good description of both patterns, brief DIIS explanation |
| **4** | Good description of patterns, weak DIIS explanation |
| **3** | Partial description of each pattern |
| **2** | Some correct observations |
| **1** | Minimal relevant observations |
| **0** | No substantive comparison |

**Key Points for Full Credit:**

*Without DIIS:*
- Energy may oscillate
- Residual decreases slowly
- May exhibit step-like or irregular pattern

*With DIIS:*
- Energy converges monotonically (or nearly so)
- Residual drops rapidly once DIIS activates
- Smooth, steep convergence curve
- "Hockey stick" shape

*DIIS Effect:*
- Smooths oscillations by extrapolating from history
- Effectively "predicts" solution

---

### Q4.4: H2O Final Energy and Convergence (4 points)

**Question:** For H2O, what is the final RHF energy? Does it converge with and without DIIS?

| Points | Criteria |
|--------|----------|
| **4** | Energy within 10^-3 of -74.963 Ha AND correct convergence assessment for both |
| **3** | Energy correct OR convergence assessment correct for both |
| **2** | Energy within 10^-2 Ha with one convergence correct |
| **1** | Order of magnitude correct (~-75 Ha) |
| **0** | Energy significantly wrong or no attempt |

**Expected Values:**
- Energy: -74.963023138435 Hartree
- Converges without DIIS: Yes (slowly)
- Converges with DIIS: Yes (faster)

---

### Q4.5: Fock Matrix Symmetry (5 points)

**Question:** Is the Fock matrix symmetric? Why is symmetry physically important?

| Points | Criteria |
|--------|----------|
| **5** | Correct observation (symmetric) AND substantive explanation (real eigenvalues, orthogonal eigenvectors) |
| **4** | Correct observation AND mentions eigenvalue properties |
| **3** | Correct observation, superficial explanation |
| **2** | Correct observation only |
| **1** | Incorrect observation with relevant discussion |
| **0** | Incorrect observation or no attempt |

**Key Points for Full Credit:**
- Yes, F_ij = F_ji (within numerical precision)
- Hermitian operators have real eigenvalues (orbital energies must be real)
- Eigenvectors are orthogonal (MOs form orthonormal basis)

---

### Q4.6: HOMO-LUMO Energies and Gap (6 points)

**Question:** What is the HOMO energy for H2? LUMO energy? What does the gap tell you?

| Points | Criteria |
|--------|----------|
| **6** | Both energies within 0.01 Ha AND meaningful gap significance |
| **5** | Both energies correct, brief significance discussion |
| **4** | Both energies correct, weak significance |
| **3** | One energy correct AND gap significance mentioned |
| **2** | One energy correct OR gap significance well explained |
| **1** | Correct identification of HOMO/LUMO concept |
| **0** | Both energies wrong or no attempt |

**Expected Values:**
- HOMO: -0.578 Ha (accept -0.57 to -0.59)
- LUMO: +0.670 Ha (accept +0.66 to +0.68)
- Gap: ~1.25 Ha (~34 eV)

**Gap Significance:**
- Large gap indicates chemical stability
- Relates to reactivity, optical properties
- H2 is an insulator (large gap)

---

### Q4.7: DIIS Recommendations (6 points)

**Question:** When would you recommend using DIIS? When might it not help or cause problems?

| Points | Criteria |
|--------|----------|
| **6** | Multiple specific recommendations AND realistic cautions with reasoning |
| **5** | Good recommendations AND at least one specific caution |
| **4** | Good recommendations OR good cautions |
| **3** | General understanding of DIIS benefits/limitations |
| **2** | Basic statement about DIIS helping convergence |
| **1** | Minimal relevant content |
| **0** | No substantive answer |

**Key Points for Full Credit:**

*Use DIIS when:*
- Most calculations (default recommendation)
- Larger molecules (more benefit)
- Oscillating or slow convergence

*Cautions:*
- Multiple SCF solutions (may bias toward one)
- Very far from minimum (simple iteration may be more stable initially)
- Near-degenerate orbitals

---

## Section 5: Synthesis (15 points)

### Q5.1: Boys Function in Integral Evaluation (5 points)

**Question:** Three parts about T dependence in nuclear attraction integrals.

| Points | Criteria |
|--------|----------|
| **5** | All three parts addressed with physical and numerical insight |
| **4** | Two parts well addressed, one partial |
| **3** | Two parts well addressed |
| **2** | One part well addressed |
| **1** | General understanding shown |
| **0** | No substantive answer |

**Expected Content:**

*(a) Small T (close centers):*
- F_m(T) near maximum value 1/(2m+1)
- Series expansion used
- Numerically straightforward

*(b) Large T (far apart):*
- F_m(T) approaches zero
- Asymptotic expansion used
- Physical: nuclear attraction decays with distance

*(c) Regime boundaries:*
- Algorithm switching must be smooth
- Both methods should agree at boundaries
- Small errors can propagate

---

### Q5.2: Quadrature and Computational Cost (5 points)

**Question:** Three parts about estimating computational cost.

| Points | Criteria |
|--------|----------|
| **5** | All three parts with quantitative reasoning and insight |
| **4** | Two parts well addressed with reasoning |
| **3** | Two parts addressed |
| **2** | One part well addressed |
| **1** | Basic cost understanding |
| **0** | No substantive answer |

**Expected Content:**

*(a) Total evaluations:*
- N_int x n (10,000 x ~8-10 = ~80,000-100,000)

*(b) Order-accuracy relationship:*
- Cost scales linearly with n
- Doubling accuracy adds ~1-2 points
- Helps budget resources

*(c) Adaptive quadrature:*
- Different integrals have different T values
- Large-T integrals need fewer points
- Adaptive selection minimizes total cost

---

### Q5.3: DIIS Summary and Explanation (5 points)

**Question:** Three parts about DIIS mechanism and applicability.

| Points | Criteria |
|--------|----------|
| **5** | Accurate summary AND insightful explanation of both (b) and (c) |
| **4** | Good summary and one well-addressed part |
| **3** | Good summary OR two parts addressed |
| **2** | Basic understanding of DIIS mechanism |
| **1** | Minimal relevant content |
| **0** | Fundamentally incorrect description |

**Expected Content:**

*(a) Summary (2-3 sentences):*
- DIIS stores history of Fock matrices and error vectors
- Finds optimal linear combination minimizing error
- Typically reduces iterations by 40-60%

*(b) Why extrapolation is effective:*
- Simple iteration can oscillate
- DIIS uses multiple iterations to predict minimum
- "Jumps over" oscillations

*(c) When standard SCF is useful:*
- Initial iterations (no history yet)
- Very close to convergence
- Exploring multiple solutions
- Simple systems

---

## Artifact Grading (Bonus: 10 points)

Students are asked to export three run artifacts. Award bonus points for correctly exported artifacts.

| Artifact | Bonus Points | Verification |
|----------|--------------|--------------|
| Boys artifact (m=3, T=8.0) | +3 | Correct parameters in JSON |
| Rys artifact (T=15.0, target=1e-8) | +3 | Correct parameters and order shown |
| SCF artifact (H2O with DIIS) | +4 | Correct system and DIIS enabled |
| **Maximum Bonus** | **+10** | |

**Verification Criteria:**
- File is valid JSON
- Contains expected module identifier
- Parameters match worksheet instructions
- Results section populated

---

## Quick Reference: Point Values

| Question | Points | Type |
|----------|--------|------|
| Q2.1 | 5 | Formula + numerical |
| Q2.2 | 5 | Numerical exploration |
| Q2.3 | 5 | Factual (from interface) |
| Q2.4 | 5 | Comparison + explanation |
| Q2.5 | 5 | Conceptual |
| Q2.6 | 5 | Interface + theory comparison |
| Q2.7 | 5 | Conceptual (implementation rationale) |
| Q3.1 | 5 | Factual (from interface) |
| Q3.2 | 5 | Conceptual |
| Q3.3 | 5 | Numerical (from interface) |
| Q3.4 | 5 | Numerical exploration |
| Q3.5 | 5 | Observation + explanation |
| Q3.6 | 5 | Shell quartet + formula verification |
| Q3.7 | 5 | Shell quartet + formula verification |
| Q3.8 | 5 | Algorithm 5.1 moments + structure |
| Q4.1 | 4 | Factual (2+2) |
| Q4.2 | 4 | Comparison (2+2) |
| Q4.3 | 6 | Descriptive |
| Q4.4 | 4 | Factual (2+2) |
| Q4.5 | 5 | Observation + explanation |
| Q4.6 | 6 | Numerical + conceptual |
| Q4.7 | 6 | Open-ended recommendation |
| Q5.1 | 5 | Multi-part synthesis |
| Q5.2 | 5 | Multi-part synthesis |
| Q5.3 | 5 | Multi-part synthesis |
| **Total** | **125** | |

**Scaled Total:** To convert to 100 points: (raw score / 125) x 100

---

## Grade Scale Recommendation

| Percentage | Letter Grade | Description |
|------------|--------------|-------------|
| 90-100 | A | Excellent understanding |
| 80-89 | B | Good understanding |
| 70-79 | C | Adequate understanding |
| 60-69 | D | Minimal understanding |
| <60 | F | Insufficient understanding |

**Note:** Instructors may adjust this scale based on class performance and institutional policies.

---

## Inter-Rater Reliability Notes

To ensure consistent grading across multiple graders:

1. **Calibration:** Have all graders score the same 3-5 sample submissions and compare
2. **Borderline cases:** When in doubt, award the higher partial credit
3. **Conceptual answers:** Focus on presence of key concepts, not specific wording
4. **Numerical tolerance:** Use the tolerances specified in this rubric, not stricter
5. **Documentation:** Note any ambiguous cases for discussion

---

*Lab Pack #1 Grading Rubric v1.0 | CONFIDENTIAL - Instructor Use Only*
*Interactive Quantum Chemistry Playground | https://iqcp.dev*
