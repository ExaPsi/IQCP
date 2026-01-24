# Lab Pack #1: Concept Check Scoring Key

**Lab Pack:** 1 - From Boys to Orbitals
**Version:** 1.0
**Last Updated:** 2026-01-18

---

## Overview

This scoring key provides comprehensive guidance for grading the 12-item pre/post concept check assessment. It includes:

- Correct answers with detailed justifications
- Distractor analysis explaining why incorrect options are wrong
- Detailed rubrics for short-answer items
- Learning outcome alignment for each item
- Common student misconceptions and how to address them

---

## Quick Reference Answer Key

### Pre-Check (8 points total)

| Item | Type | Correct | Points |
|------|------|---------|--------|
| P1 | MC | **A** | 1 |
| P2 | MC | **B** | 1 |
| P3 | MC | **B** | 1 |
| P4 | MC | **B** | 1 |
| P5 | SA | See rubric | 2 |
| P6 | SA | See rubric | 2 |

### Post-Check (9 points total)

| Item | Type | Correct | Points |
|------|------|---------|--------|
| Q1 | MC | **A** | 1 |
| Q2 | MC | **B** | 1 |
| Q3 | MC | **D** | 1 |
| Q4 | MC | **A** | 1 |
| Q5 | MC | **B** | 1 |
| Q6 | SA | See rubric | 2 |
| Q7 | SA | See rubric | 2 |

**Total Assessment: 17 points**

---

## Pre-Check Detailed Answer Key

### Item P1: Quadrature Concepts (1 point)

**Question:** When evaluating a definite integral numerically using quadrature, increasing the number of quadrature points generally...

**Correct Answer: A** - Increases accuracy but increases computational cost

**Detailed Justification:**

This is the fundamental trade-off in numerical integration. Each additional quadrature point:
- Requires one more function evaluation (cost increases linearly)
- For Gaussian quadrature, allows exact integration of polynomials up to degree 2n-1
- Reduces truncation error for smooth functions
- The relationship is well-established: O(n) cost for O(1/n^k) error improvement

**Learning Outcome Alignment:** LO2 (Rys quadrature order-accuracy relationship)

---

#### P1 Distractor Analysis

| Option | Why Students Choose It | Why It Is Wrong |
|--------|------------------------|-----------------|
| **B** - "Decreases accuracy but increases cost" | Confusion between numerical precision (finite arithmetic) and mathematical accuracy (truncation error). Some students may recall warnings about floating-point accumulation. | While round-off error can accumulate, this effect is negligible for typical quadrature orders (n < 20). The dominant error is truncation, which DECREASES with more points. |
| **C** - "Increases accuracy and decreases cost" | Wishful thinking. Students may believe efficiency improves universally with better algorithms. | Violates the "no free lunch" principle. Every additional evaluation requires computation time. There is no magical efficiency gain. |
| **D** - "No effect on accuracy" | Complete misunderstanding of quadrature. May result from unfamiliarity with numerical methods. | Quadrature accuracy is fundamentally determined by the number of points. This is the core theorem of Gaussian quadrature. |

---

### Item P2: Iterative Methods (1 point)

**Question:** In an iterative numerical method, "convergence" means...

**Correct Answer: B** - Successive approximations approach a fixed answer

**Detailed Justification:**

Convergence is a mathematical concept describing the behavior of a sequence:
- Formally: lim(n->infinity) x_n = x* for some fixed point x*
- In practice: |x_{n+1} - x_n| < tolerance for some iteration n
- The SCF method iteratively refines the density matrix and Fock matrix until self-consistency
- Convergence does NOT imply speed, memory usage, or initial guess quality

**Learning Outcome Alignment:** LO3 (SCF convergence and DIIS acceleration)

---

#### P2 Distractor Analysis

| Option | Why Students Choose It | Why It Is Wrong |
|--------|------------------------|-----------------|
| **A** - "Runs faster with each iteration" | Confusing "convergence" (mathematical) with "speedup" (computational). May conflate with concepts like adaptive algorithms. | Convergence refers to the sequence of VALUES approaching a limit, not to computational speed. Individual iterations may take the same or varying time. |
| **C** - "Uses less memory over time" | Similar confusion with computational resources. May think iterative methods discard old data. | Memory usage is an implementation detail unrelated to convergence. DIIS actually stores MORE data (previous Fock matrices) to accelerate convergence. |
| **D** - "Initial guess was correct" | Misunderstanding causality. If the initial guess is exact, convergence is trivial (0 iterations). | Convergence describes the PROCESS of reaching the answer, not whether the starting point was already correct. The definition applies regardless of initial guess quality. |

---

### Item P3: Gaussian Integrals (1 point)

**Question:** The integral of exp(-x^2) from 0 to infinity is...

**Correct Answer: B** - A finite positive value

**Detailed Justification:**

This is the half Gaussian integral:
- Full integral: integral(-inf to +inf) of exp(-x^2) dx = sqrt(pi)
- Half integral: integral(0 to +inf) of exp(-x^2) dx = sqrt(pi)/2 approximately equals 0.886
- The integrand is always positive (exp of any real number > 0)
- The exponential decay exp(-x^2) dominates any polynomial growth, ensuring convergence
- This integral is fundamental to quantum chemistry: Gaussian basis functions have this form

**Learning Outcome Alignment:** LO1 (Boys function behavior and numerical regimes)

---

#### P3 Distractor Analysis

| Option | Why Students Choose It | Why It Is Wrong |
|--------|------------------------|-----------------|
| **A** - "Undefined (diverges)" | May confuse with integrals of polynomials or exp(+x^2). Students may not recognize how fast exp(-x^2) decays. | The Gaussian exp(-x^2) decays faster than ANY inverse polynomial. It goes to zero faster than 1/x^n for any n, guaranteeing convergence. This is proven rigorously using the squeeze theorem. |
| **C** - "Zero" | May confuse with odd functions integrated over symmetric domains, or with misremembered formulas. | The integrand exp(-x^2) is strictly positive for all x. The integral of a positive function over a positive-measure domain must be positive. |
| **D** - "Negative infinity" | Complete misunderstanding. May think integrals to infinity must diverge to plus or minus infinity. | Impossible. exp(-x^2) > 0 for all x, so the integral must be positive (or zero if the domain has measure zero, which [0,inf) does not). |

---

### Item P4: Hartree-Fock Theory (1 point)

**Question:** In quantum chemistry, the Hartree-Fock method calculates...

**Correct Answer: B** - An approximate molecular energy (variational upper bound)

**Detailed Justification:**

Key concepts in Hartree-Fock theory:
- HF assumes electrons move in an average field of other electrons (mean-field approximation)
- This neglects instantaneous electron correlation
- The variational principle guarantees E_HF >= E_exact
- Typical HF recovers ~99% of total energy but misses correlation energy
- For H2 (STO-3G), E_HF approximately equals -1.117 Ha vs. E_exact approximately equals -1.174 Ha

**Learning Outcome Alignment:** LO3 (SCF convergence and DIIS acceleration)

---

#### P4 Distractor Analysis

| Option | Why Students Choose It | Why It Is Wrong |
|--------|------------------------|-----------------|
| **A** - "Exact molecular energy" | May overestimate HF capabilities or confuse with full configuration interaction (FCI). | HF uses a single Slater determinant, which cannot capture electron correlation. Exact energy requires full many-body treatment. The "correlation energy" E_corr = E_exact - E_HF is always negative. |
| **C** - "Only kinetic energy" | Confusing HF with the kinetic energy operator T. May not understand how the Fock operator is constructed. | The Fock operator includes: kinetic energy (T), nuclear attraction (V_ne), Coulomb repulsion (J), and exchange (K). HF computes the full electronic energy, not just one component. |
| **D** - "Only nuclear repulsion" | Confusing HF with classical electrostatics. May not realize that nuclear repulsion (E_nuc) is computed separately and exactly. | Nuclear repulsion is computed EXACTLY as E_nuc = sum(Z_A * Z_B / R_AB). HF focuses on the ELECTRONIC energy. The total energy is E_total = E_elec + E_nuc. |

---

### Item P5: Numerical Regimes (2 points)

**Question:** In one or two sentences, explain why numerical algorithms might use different computational methods (e.g., series expansion vs. asymptotic formula) for different parameter ranges.

**Learning Outcome Alignment:** LO1, LO4 (Boys function regimes, parameter-outcome connections)

---

#### P5 Detailed Rubric

| Points | Criteria | Example Responses |
|--------|----------|-------------------|
| **2** | Correctly explains that different methods have different convergence, stability, or efficiency properties in different ranges. Must mention at least TWO of: (a) convergence rate, (b) numerical stability, (c) computational efficiency. | "Series expansions converge quickly for small parameters but require many terms for large values, while asymptotic formulas work well for large parameters but diverge for small ones. Using the right method for each range gives both accuracy and efficiency." |
| **1** | Partially correct: mentions that different methods work in different ranges but does not clearly explain WHY (convergence, stability, or efficiency). Only addresses one of the three factors. | "Different formulas are more accurate in different ranges." (Correct but lacks detail about WHY they are more accurate.) |
| **0** | Incorrect, irrelevant, or no meaningful response. | "To make the code run faster." (Does not address accuracy/convergence trade-offs.) |

---

#### P5 Key Concepts to Look For

**Full Credit (2 points) - Student should mention:**

1. **Convergence behavior differs by regime:**
   - Series expansions converge rapidly for small T (few terms needed)
   - Series converge slowly or diverge for large T (many terms, numerical issues)
   - Asymptotic expansions are accurate for large T (first few terms suffice)
   - Asymptotic series diverge for small T

2. **Numerical stability concerns:**
   - Series for large T: alternating large terms lead to catastrophic cancellation
   - Asymptotic for small T: division by small numbers causes instability

3. **Computational efficiency:**
   - Using the optimal method for each regime minimizes computation
   - Avoids unnecessary iterations/terms

**Partial Credit (1 point):**
- Recognizes that different methods are suited to different ranges
- Does not explain the mathematical or computational reasons

**No Credit (0 points):**
- States that it is "just how it is done"
- Confuses with unrelated concepts
- No response

---

### Item P6: Convergence Tolerance (2 points)

**Question:** What does it mean for an iterative calculation to "converge to a tolerance of 10^-6"? Explain in one or two sentences.

**Learning Outcome Alignment:** LO3, LO4 (SCF convergence, parameter-outcome connections)

---

#### P6 Detailed Rubric

| Points | Criteria | Example Responses |
|--------|----------|-------------------|
| **2** | Correctly explains that consecutive iterations (or the quantity of interest vs. its limit) differ by less than 10^-6. May specify energy, density matrix, or other convergence metric. | "The calculation stops when the energy difference between consecutive iterations is less than 10^-6 Hartree. This means the answer has stabilized to at least six decimal places." |
| **1** | Mentions the tolerance value but does not clearly specify what is being compared (consecutive values? value vs. limit?) or is vague about the comparison. | "The calculation is accurate to 10^-6." (Correct direction but does not explain the mechanism of checking convergence.) |
| **0** | Incorrect or no meaningful response. | "The calculation uses 6 iterations." (Confuses tolerance with iteration count.) |

---

#### P6 Key Concepts to Look For

**Full Credit (2 points) - Student should explain:**

1. **What is compared:**
   - Energy: |E_{n+1} - E_n| < 10^-6
   - Density: ||P_{n+1} - P_n|| < 10^-6
   - Or generally: difference between successive approximations

2. **What the threshold means:**
   - When the change is smaller than 10^-6, we consider the result "converged"
   - The value has stabilized to within that precision
   - Further iterations would not meaningfully change the result

**Partial Credit (1 point):**
- Mentions "tolerance" or "10^-6" in the context of accuracy
- Does not clearly explain the comparison mechanism

**No Credit (0 points):**
- Confuses tolerance with iteration count
- Confuses tolerance with precision (significant figures)
- No response

---

## Post-Check Detailed Answer Key

### Item Q1: Rys Quadrature (1 point)

**Question:** In Rys quadrature for molecular integrals, increasing the quadrature order from n=5 to n=8...

**Correct Answer: A** - Increases accuracy by allowing exact integration of higher-order polynomial moments

**Detailed Justification:**

Rys quadrature with n points:
- Exactly integrates polynomials up to degree 2n-1
- n=5 exactly integrates up to degree 9
- n=8 exactly integrates up to degree 15
- Higher angular momentum integrals require higher polynomial accuracy
- The observed error decreases as order increases (students see this in the Rys module)

**Learning Outcome Alignment:** LO2 (Rys quadrature order-accuracy relationship)

---

#### Q1 Distractor Analysis

| Option | Why Students Choose It | Why It Is Wrong |
|--------|------------------------|-----------------|
| **B** - "More points = more round-off error" | Over-generalizing floating-point concerns. May recall warnings about accumulated numerical error. | At these modest orders (n=5-8), round-off error is negligible compared to truncation error. IEEE-754 double precision has ~15 decimal digits; we need only ~10-12 for chemistry accuracy. |
| **C** - "No effect - domain [0,1] is fixed" | Confusing the integration domain with the quadrature accuracy. The domain is indeed fixed, but accuracy depends on ORDER. | The domain [0,1] is a property of the Rys polynomial formulation. Accuracy depends on how well the quadrature approximates the integrand, which is determined by order, not domain. |
| **D** - "Only for small T" | Incomplete understanding. May have observed that accuracy behavior differs with T but did not understand the general relationship. | Higher order improves accuracy for ALL T values, not just small ones. The relationship between order and accuracy is fundamental to Gaussian quadrature theory. |

---

### Item Q2: DIIS Acceleration (1 point)

**Question:** DIIS (Direct Inversion in the Iterative Subspace) improves SCF convergence by...

**Correct Answer: B** - Extrapolating an optimal Fock matrix from previous iterations

**Detailed Justification:**

DIIS algorithm:
1. Store Fock matrices F_1, F_2, ..., F_k from previous iterations
2. Compute error vectors e_i = F_i * D_i * S - S * D_i * F_i
3. Solve B * c = [0, 0, ..., -1]^T for coefficients c_i
4. Form extrapolated Fock: F_extrap = sum(c_i * F_i)
5. Use F_extrap for next iteration

The extrapolation minimizes the residual norm, effectively predicting a Fock matrix closer to the converged solution.

**Learning Outcome Alignment:** LO3 (SCF convergence and DIIS acceleration)

---

#### Q2 Distractor Analysis

| Option | Why Students Choose It | Why It Is Wrong |
|--------|------------------------|-----------------|
| **A** - "Faster individual iterations" | Confusing "acceleration" (fewer iterations to converge) with "speedup" (faster computation per iteration). | DIIS actually adds computation per iteration (storing matrices, solving linear system). The benefit is FEWER total iterations, not faster ones. |
| **C** - "Fewer basis functions" | Confusing DIIS with basis set truncation or other approximations. May not understand what DIIS actually does. | DIIS does not change the basis set at all. It operates on the Fock/density matrices, which are of fixed dimension determined by the basis. |
| **D** - "Changes molecular geometry" | Confusing SCF (electronic structure) with geometry optimization. May not distinguish between different types of "optimization." | DIIS is purely an electronic structure acceleration technique. Geometry optimization uses gradients and is a separate outer loop. |

---

### Item Q3: Boys Function Asymptotics (1 point)

**Question:** For large values of T, the Boys function F_m(T) approaches...

**Correct Answer: D** - Zero

**Detailed Justification:**

Asymptotic behavior of F_m(T):
- F_m(T) = integral(0 to 1) of t^(2m) * exp(-T * t^2) dt
- As T increases, the exponential weight exp(-T * t^2) decays rapidly
- The dominant contribution shifts to t near 0
- Asymptotic formula: F_m(T) approximately equals (2m-1)!! / (2^(m+1)) * sqrt(pi / T^(2m+1))
- As T -> infinity, F_m(T) -> 0 for all m >= 0
- Students directly observe this decay in the Boys module visualization

**Learning Outcome Alignment:** LO1 (Boys function behavior and numerical regimes)

---

#### Q3 Distractor Analysis

| Option | Why Students Choose It | Why It Is Wrong |
|--------|------------------------|-----------------|
| **A** - "Infinity" | Incorrectly extrapolating. May think that "going to large T" means "values get large." | The exponential weight exp(-T * t^2) SUPPRESSES the integrand for large T. The integral shrinks, not grows. |
| **B** - "One" | Confusing with F_0(0) = 1. May misremember the limiting behavior. | F_m(0) = 1/(2m+1), with F_0(0) = 1. But for LARGE T, all F_m(T) -> 0, not 1. |
| **C** - "1/(2m+1)" | This is F_m(0), not F_m(large T). Students may confuse which limit is being asked about. | F_m(0) = 1/(2m+1) is the small-T limit (no exponential suppression). The large-T behavior is fundamentally different: F_m(T) -> 0. |

---

### Item Q4: Method Selection (1 point)

**Question:** IQCP uses a Taylor series for F_m(T) when T is small, but switches to a recurrence relation (erf + upward recurrence) when T exceeds a threshold. This is primarily because...

**Correct Answer: A** - The Taylor series converges too slowly (requires too many terms) for large T

**Detailed Justification:**

Series behavior for Boys function:
- Taylor series: F_m(T) = exp(-T) * sum(k=0 to inf) of T^k / (2m + 2k + 1)!!
- For small T: exp(-T) is near 1, terms decrease rapidly, few needed
- For large T: exp(-T) is tiny, but T^k grows rapidly; need many terms before exp(-T) * T^k decreases
- This leads to catastrophic cancellation: subtracting large nearly-equal numbers
- The recurrence method using erf(sqrt(T)) + upward recurrence is efficient for moderate-to-large T

**IQCP Implementation Note:** IQCP follows libcint in using only TWO computational methods (series and recurrence), with m-dependent turnover points. There is no separate asymptotic regime.

**Learning Outcome Alignment:** LO1, LO4 (Boys function methods, parameter-outcome connections)

---

#### Q4 Distractor Analysis

| Option | Why Students Choose It | Why It Is Wrong |
|--------|------------------------|-----------------|
| **B** - "Recurrence always more accurate" | May overgeneralize from the large-T case. Does not understand that series can be more accurate for small T. | The recurrence method is NOT always more accurate. For very small T, the series expansion can be more efficient and equally accurate. Each method has its appropriate domain determined by the m-dependent turnover point. |
| **C** - "Taylor series is faster" | May think series = fast, other methods = slow. Does not understand the slow convergence issue. | The Taylor series requires MANY terms for large T (slow convergence), making it SLOWER than the recurrence method, which computes erf once and then applies a simple recurrence. |
| **D** - "Integer overflow" | Technical misconception. May think T values cause computational errors. | Integer overflow is not the issue. The problem is mathematical: series convergence, not data type limitations. Modern codes use double precision which handles the relevant T range. |

---

### Item Q5: Root Count Rule (1 point)

**Question:** In Rys quadrature for electron repulsion integrals, the number of quadrature roots needed depends on the shell quartet. For a (pp|pp) shell quartet (L=4), how many roots are required according to the formula n_r = floor(L/2) + 1?

**Correct Answer: B** - 3

**Detailed Justification:**

For a (pp|pp) shell quartet:
- Each p shell has angular momentum l=1
- Total angular momentum L = l_A + l_B + l_C + l_D = 1 + 1 + 1 + 1 = 4
- Applying the formula: n_r = floor(4/2) + 1 = 2 + 1 = 3

The formula n_r = floor(L/2) + 1 ensures that the Rys quadrature can exactly integrate polynomials up to degree 2n_r - 1 = L, which is the highest polynomial degree appearing in the electron repulsion integral integrand for that shell quartet.

**Learning Outcome Alignment:** LO2 (Rys quadrature order-accuracy relationship)

---

#### Q5 Distractor Analysis

| Option | Why Students Choose It | Why It Is Wrong |
|--------|------------------------|-----------------|
| **A** - "2" | May compute L/2 = 2 but forget the +1 | The formula is n_r = floor(L/2) + 1, not just floor(L/2). The +1 ensures at least one quadrature point even for L=0. |
| **C** - "4" | May think n_r = L | The formula divides by 2 before adding 1. Using n_r = L would significantly over-estimate the needed quadrature points. |
| **D** - "5" | May add L+1 or make other arithmetic error | Neither L+1=5 nor floor(L/2)+2=4 is the correct formula. |

---

### Item Q6: Quadrature Order Selection (2 points)

**Question:** Based on your observations in the Rys module, explain how the optimal quadrature order depends on both the parameter T and the target accuracy. (2-3 sentences)

**Learning Outcome Alignment:** LO2, LO4 (Rys quadrature, parameter-outcome connections)

---

#### Q6 Detailed Rubric

| Points | Criteria | Example Responses |
|--------|----------|-------------------|
| **2** | Correctly explains BOTH relationships: (1) Higher accuracy targets require higher quadrature order, AND (2) The relationship between T and optimal order (larger T generally allows smaller order for same accuracy, or order-T interaction). May include specific observations from the Rys module. | "Higher target accuracy requires more quadrature points because each additional point allows exact integration of higher-degree polynomial terms. For larger T values, the integrand becomes more localized near t=0, which can actually require FEWER points for the same accuracy since the function is simpler in that region." |
| **1** | Correctly explains ONE relationship (either accuracy-order OR T-order) but not both, OR explains both but incompletely or with minor errors. | "More points give better accuracy" (correct but incomplete - does not address T dependence). OR "Larger T needs more points" (partially addresses T but may be reversed for some cases). |
| **0** | Incorrect or no meaningful response. Does not demonstrate understanding from the lab activity. | "The order does not matter much." |

---

#### Q6 Key Concepts to Look For

**Full Credit (2 points) - Student should explain:**

1. **Accuracy-order relationship:**
   - Higher order = exact integration of higher-degree polynomials
   - Tighter tolerance targets require more quadrature points
   - Observed in Rys module: error curves show decreasing error with increasing n

2. **T-order relationship:**
   - For large T, exp(-T * t^2) localizes integrand near t = 0
   - Localized integrands may need fewer points
   - For small T, the full [0,1] domain contributes, potentially needing more points
   - OR: the relationship is complex and depends on specific accuracy targets

**Partial Credit (1 point):**
- Correctly identifies one relationship
- Vague or incomplete explanation of both

**No Credit (0 points):**
- Completely incorrect understanding
- No reference to observations from the Rys module

---

### Item Q7: DIIS Effect on Convergence (2 points)

**Question:** Describe in 2-3 sentences how DIIS changes the SCF convergence behavior compared to standard (non-accelerated) iteration. What specific differences did you observe?

**Learning Outcome Alignment:** LO3 (SCF convergence and DIIS acceleration)

---

#### Q7 Detailed Rubric

| Points | Criteria | Example Responses |
|--------|----------|-------------------|
| **2** | Correctly describes DIIS effects including: (1) Reduces total number of iterations, AND (2) Changes convergence pattern (smoother, avoids oscillations, faster approach to minimum). May include specific iteration counts or observations from the SCF module. | "With DIIS enabled, H2O converged in about 8 iterations instead of 15 without DIIS. The energy values approached the final answer more smoothly without the oscillations I observed in standard SCF, where the energy bounced above and below the converged value." |
| **1** | Mentions one aspect of DIIS improvement (fewer iterations OR smoother convergence) but not both, OR description is vague without specific observations. | "DIIS made it converge faster." (Correct but lacks detail about convergence pattern or specific observations.) |
| **0** | Incorrect or no meaningful response. Does not demonstrate understanding from the lab activity. | "DIIS changes the final energy." (Incorrect - DIIS does not change the converged result, only how fast we get there.) |

---

#### Q7 Key Concepts to Look For

**Full Credit (2 points) - Student should describe:**

1. **Reduced iteration count:**
   - Specific numbers if possible (e.g., "15 iterations without DIIS, 8 with DIIS")
   - General statement about "fewer iterations needed"

2. **Changed convergence pattern:**
   - Standard SCF may oscillate or converge slowly (linear convergence)
   - DIIS produces smoother, more monotonic convergence
   - May mention "avoiding oscillations" or "faster approach"

3. **Specific observations from the lab:**
   - Reference to the H2O preset or other systems
   - Mention of observing the iteration plots/tables

**Partial Credit (1 point):**
- Only mentions "faster" without explaining how
- Does not include specific observations

**No Credit (0 points):**
- Confuses DIIS with other concepts
- Says DIIS changes the final answer (it does not)
- No response

---

## Learning Outcome Alignment Summary

| Learning Outcome | Description | Pre-Check Items | Post-Check Items |
|------------------|-------------|-----------------|------------------|
| **LO1** | Understand Boys function behavior and numerical methods | P3 (Gaussian integral), P5 (numerical methods) | Q3 (asymptotics), Q4 (method selection) |
| **LO2** | Connect quadrature order to integration accuracy | P1 (quadrature concepts) | Q1 (Rys order), Q5 (optimal order) |
| **LO3** | Understand SCF convergence and DIIS acceleration | P2 (iterative methods), P6 (tolerance) | Q2 (DIIS mechanism), Q6 (DIIS effect) |
| **LO4** | Make parameter-outcome connections | P5 (method rationale), P6 (tolerance meaning) | Q4 (method selection), Q5 (order-T-accuracy) |

---

## Common Misconceptions Reference

| Misconception | Appears In | Correct Understanding |
|---------------|------------|----------------------|
| "More computation always means more error" | P1, Q1 | Round-off error is negligible at typical quadrature orders; truncation error dominates and decreases with more points. |
| "Convergence means faster computation" | P2 | Convergence describes the mathematical behavior of a sequence approaching a limit, not computational speed. |
| "All integrals to infinity diverge" | P3 | Exponentially decaying integrands like exp(-x^2) give finite integrals. |
| "Hartree-Fock is exact" | P4 | HF neglects electron correlation and gives a variational upper bound. |
| "One formula works everywhere" | P5, Q4 | Different numerical methods have different convergence properties in different parameter ranges. The m-dependent turnover in IQCP shows this sophistication. |
| "Tolerance = number of iterations" | P6 | Tolerance is the threshold for the change between iterations, not the iteration count. |
| "DIIS changes the final answer" | Q2, Q6 | DIIS changes how fast we reach the answer, not what the answer is. |
| "F_m(large T) = 1/(2m+1)" | Q3 | That is F_m(0). For large T, F_m(T) approaches zero. |

---

## Score Interpretation Guide

### Individual Item Analysis

| Score Pattern | Interpretation | Recommended Action |
|---------------|----------------|-------------------|
| High pre, high post | Student had strong prior knowledge | May benefit from advanced extensions |
| Low pre, high post | Strong learning gain | Activity working as intended |
| Low pre, low post | Learning difficulty | Review pre-requisites; provide additional support |
| High pre, low post | Possible confusion from activity | Check for misconceptions introduced |

### Class-Level Analysis

| Metric | Calculation | Target |
|--------|-------------|--------|
| Pre-check mean | Sum of pre scores / n | N/A (baseline) |
| Post-check mean | Sum of post scores / n | > 6/8 (75%) |
| Normalized gain | (Post - Pre) / (8 - Pre) | > 0.3 (medium gain) |
| Item difficulty | (Correct responses) / n | 0.3 - 0.8 (optimal) |

### Normalized Gain Interpretation

| Gain (g) | Classification | Typical Causes |
|----------|----------------|----------------|
| g > 0.7 | High | Excellent activity design; engaged students |
| 0.3 < g < 0.7 | Medium | Typical for well-designed interactive activities |
| g < 0.3 | Low | Activity issues; assessment misalignment; implementation problems |

---

## Administration Checklist

Before administering the concept check:

- [ ] Print sufficient copies of pre-check and post-check forms
- [ ] Ensure students cannot access IQCP during assessment
- [ ] Plan timing: 10-12 minutes per phase
- [ ] Prepare answer key for quick grading
- [ ] Have backup plan if technology fails during lab (affects post-check validity)

After administration:

- [ ] Grade using this scoring key
- [ ] Calculate class-level statistics
- [ ] Identify items with low post-check accuracy for follow-up instruction
- [ ] Archive responses for assessment validity documentation

---

*IQCP Lab Pack #1 v1.0 | Concept Check Scoring Key | https://iqcp.dev*
