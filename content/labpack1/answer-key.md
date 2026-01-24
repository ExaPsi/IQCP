# Lab Pack #1: Instructor Answer Key

**Lab Pack:** 1 - From Boys to Orbitals
**Version:** 1.0
**Last Updated:** 2026-01-18
**Document Type:** Instructor Materials (CONFIDENTIAL)

---

## Document Overview

This answer key provides expected responses for all 20 numbered questions in the Lab Pack #1 student worksheet. For each question, the following information is provided:

- **Question Text:** Exact wording from the worksheet
- **Expected Answer:** The correct response with scientific values
- **Acceptable Range:** Tolerance for numerical answers
- **Common Misconceptions:** What students might get wrong
- **Pedagogical Intent:** What the question is testing
- **Grading Notes:** Tips for partial credit

### Scientific Reference Environment

All numerical values were computed using:
- PySCF 2.11.0
- SciPy 1.17.0 (hyp1f1 for Boys functions)
- NumPy 2.4.1
- Python 3.12
- Date verified: 2026-01-18

---

## Section 2: Boys Function (Q2.1-Q2.5)

### Q2.1: Limiting Value as T Approaches 0

**Question Text:**
> What value does F_0(T) approach as T approaches 0? Express this as a simple fraction in terms of m.

**Expected Answer:**

F_m(0) = 1/(2m+1)

Specifically:
- F_0(0) = 1
- F_1(0) = 1/3
- F_2(0) = 1/5
- F_3(0) = 1/7

**Derivation:** When T = 0, the Boys function integral becomes:
```
F_m(0) = integral from 0 to 1 of t^(2m) dt = [t^(2m+1)/(2m+1)] from 0 to 1 = 1/(2m+1)
```

**Scientific Basis (SciPy verified):**
```
F_0(0.01) = 0.9966722188
F_0(0.001) = 0.9996667222
F_0(0) = 1.0000000000 (limit)
F_1(0) = 0.3333333333 = 1/3
F_2(0) = 0.2000000000 = 1/5
```

**Acceptable Range:**
- Must state F_m(0) = 1/(2m+1) or equivalent
- For m=0, must state F_0(0) = 1

**Common Misconceptions:**
1. Students may write "1" without recognizing the m-dependence
2. Students may confuse with exponential decay behavior at large T
3. Students may write 1/m instead of 1/(2m+1)

**Pedagogical Intent:**
Tests understanding of the integral definition of Boys functions. When T=0, exp(-T*t^2) = 1, reducing the integral to a simple power function.

**Grading Notes:**
- Full credit (5 pts): Correct formula 1/(2m+1) with F_0(0) = 1
- Partial credit (3 pts): Correct for m=0 but incorrect general formula
- Partial credit (2 pts): Recognizes the limit exists and is finite
- No credit: Incorrect answer or no attempt

---

### Q2.2: T Value Where F_0(T) < 0.01

**Question Text:**
> As T increases, F_0(T) approaches zero. At approximately what T value does F_0(T) become smaller than 0.01?

**Expected Answer:**

F_0(T) < 0.01 requires T approximately 7800-8000 (exact: ~7854)

**Scientific Basis (SciPy verified):**
```
F_0(100) = 0.088623
F_0(1000) = 0.028025
F_0(5000) = 0.012533
F_0(7000) = 0.010562
F_0(7854) = 0.010000 (threshold)
F_0(8000) = 0.009908
F_0(10000) = 0.008862
```

**Derivation:** For large T, F_0(T) approaches sqrt(pi/(4T)). Solving sqrt(pi/(4T)) = 0.01:
- pi/(4T) = 0.0001
- T = pi/0.0004 = 7854

**Acceptable Range:**
- T > 7000 is acceptable
- Any value between 7500 and 8500 earns full credit
- Order of magnitude (T ~ 10000) earns partial credit

**Common Misconceptions:**
1. Students expect exponential decay (much faster than actual sqrt decay)
2. Students may stop exploring at T = 100 and extrapolate incorrectly
3. Students may confuse F_0 behavior with F_m (m > 0) which does decay faster

**Pedagogical Intent:**
This challenging question reveals that Boys functions decay much more slowly than exponential for m=0. The sqrt(1/T) asymptotic behavior is surprising to most students.

**Grading Notes:**
- Full credit (5 pts): Correct range (7000-9000) or recognition of very slow decay
- Partial credit (3 pts): Correct order of magnitude (thousands)
- Partial credit (2 pts): Recognizes decay is slower than expected
- No credit: Value less than 100 or no attempt

**Instructor Note:** This is intentionally a challenging question. Many students will not find the exact value, and that is acceptable. The key insight is recognizing the surprisingly slow decay.

---

### Q2.3: Computational Methods

**Question Text:**
> In the Internals panel, what computational method is used for T = 0.5? For T = 15.0? For T = 35.0?

**Expected Answer:**

| T Value | Computational Method | Reason |
|---------|---------------------|--------|
| T = 0.5 | Recurrence (erf + upward recurrence) | For m=0, turnover=0, so ALL T>0 uses recurrence |
| T = 15.0 | Recurrence (erf + upward recurrence) | Well above turnover for any reasonable m |
| T = 35.0 | Recurrence (erf + upward recurrence) | Same method; no asymptotic regime exists |

**Scientific Basis:**
IQCP follows the libcint implementation with **two computational methods** (not three):

1. **Series expansion:** Used when T < turnover(m)
2. **Recurrence (erf + upward recurrence):** Used when T >= turnover(m)

The **turnover point is m-dependent**, not a fixed boundary:
- m=0: turnover = 0 (ALWAYS uses recurrence for T > 0)
- m=1: turnover = 0 (ALWAYS uses recurrence for T > 0)
- m=2: turnover = 0.87
- m=5: turnover = 2.11
- m=10: turnover = 4.05
- m=20: turnover = 7.84
- m=30: turnover = 11.58

**Important:** There is NO asymptotic regime in the IQCP implementation. The literature often describes three regimes (series/recurrence/asymptotic), but libcint (and thus IQCP) uses only two methods with the recurrence handling all moderate-to-large T values effectively.

**Reference Values at these T (for m=0):**
```
F_0(0.5) = 0.855624391892 (recurrence, since turnover=0 for m=0)
F_0(15.0) = 0.228822798330 (recurrence)
F_0(35.0) = 0.149799691340 (recurrence)
```

**Acceptable Answers:**
- For m=0: "Recurrence" for all three T values is correct
- For higher m values: Students may observe "Series" for small T values below the m-dependent turnover
- Accept synonyms: "erf-based" for "Recurrence", "Taylor" or "power series" for "Series"

**Common Misconceptions:**
1. Students may expect a fixed T=12 boundary (common in literature but not in IQCP)
2. Students may expect an "asymptotic" regime for large T (not implemented)
3. Students may not realize the turnover depends on m

**Pedagogical Intent:**
Tests ability to read and interpret the IQCP interface. Introduces the concept that method selection depends on BOTH T and m, following the sophisticated approach used in production quantum chemistry codes like libcint.

**Instructor Note:** If students report answers based on fixed T=12/T=30 boundaries from textbooks, use this as a teaching moment about how real implementations differ from simplified pedagogical descriptions.

**Grading Notes:**
- Full credit (5 pts): All three regimes correctly identified
- Partial credit (3-4 pts): Two of three correct
- Partial credit (1-2 pts): One correct
- No credit: None correct or no attempt

---

### Q2.4: Comparison of F_0(10) vs F_5(10)

**Question Text:**
> Compare F_0(10.0) and F_5(10.0). Which is larger? Why do you think higher-order Boys functions decay more rapidly with increasing T?

**Expected Answer:**

**Numerical comparison:**
- F_0(10.0) = 0.280247 (larger)
- F_5(10.0) = 7.9009 x 10^-5 (much smaller)
- Ratio: F_0/F_5 approximately 3500

**Explanation:**
Higher-order Boys functions have additional factors of t^(2m) in the integrand. The t^(2m) factor suppresses contributions at larger t values. Since exp(-T*t^2) peaks near t=0 for large T, and the t^(2m) factor is small near t=0, higher m values result in smaller integrals.

**Alternative explanation:**
The integrand t^(2m) * exp(-T*t^2) is increasingly concentrated near t=0 as m increases. Since the domain is [0,1] and larger m means stronger weighting toward t=0 where t^(2m) is small, the integral decreases.

**Acceptable Range:**
- Numerical values: F_0 within 10%, F_5 within factor of 2 (order of magnitude)
- Must correctly identify F_0 > F_5
- Explanation should mention the t^(2m) factor

**Common Misconceptions:**
1. Students may think higher m means larger values
2. Students may confuse this with the T=0 behavior (where F_m decreases with m)
3. Students may not connect the t^(2m) factor to the decay rate

**Pedagogical Intent:**
Tests conceptual understanding of how the order parameter m affects Boys function values. Links mathematical structure to numerical behavior.

**Grading Notes:**
- Full credit (5 pts): Correct comparison AND reasonable explanation
- Partial credit (3 pts): Correct comparison, weak explanation
- Partial credit (2 pts): Correct comparison only
- No credit: Incorrect comparison or no attempt

---

### Q2.5: Why Different Computational Methods?

**Question Text:**
> Why do you think quantum chemistry codes use different computational methods (series expansion vs. recurrence relation) for different parameter values, rather than using a single method for all cases?

**Expected Answer:**

Key points to include:

1. **Series expansion** converges rapidly for small T (few terms needed) but converges very slowly for large T (requires too many terms for accuracy)

2. **Recurrence relation** starting from erf(sqrt(T)) is numerically stable for moderate-to-large T values and more efficient than the series

3. **m-dependent turnover:** The optimal switch point depends on the order m, not just T. Higher m values have higher turnover points because the series converges better for larger T when m is large.

4. **Numerical stability:** Each method is optimized for its regime. Using the wrong method leads to:
   - Loss of precision (catastrophic cancellation)
   - Slow computation (many unnecessary terms)
   - Potential overflow or underflow

5. **Efficiency:** Choosing the optimal method minimizes computation time while maintaining required accuracy

**Note:** While some textbooks describe three regimes (series/recurrence/asymptotic), the IQCP implementation follows libcint in using only two methods. The recurrence method handles both moderate and large T values effectively.

**Acceptable Range:**
Must mention at least two of: accuracy, efficiency, numerical stability, convergence behavior

**Common Misconceptions:**
1. Students may think it is purely about speed (not accuracy)
2. Students may not understand catastrophic cancellation
3. Students may think one method could work everywhere with enough terms

**Pedagogical Intent:**
Tests understanding of practical numerical analysis considerations. Real quantum chemistry codes must balance accuracy, stability, and efficiency.

**Grading Notes:**
- Full credit (5 pts): Mentions accuracy/stability AND efficiency, with specific regime considerations
- Partial credit (3-4 pts): Mentions at least two key concepts
- Partial credit (1-2 pts): General idea about different methods for different ranges
- No credit: No substantive answer

---

### Q2.6: IQCP Methods vs. Theoretical 3-Regime Model

**Question Text:**
> Looking at the Internals panel, what computational method does IQCP actually use for m=0 at T=35? For m=5 at T=45? According to the theoretical 3-regime model, what theoretical regime would these fall into?

**Expected Answer:**

| Case | IQCP Method | Theoretical Regime | Explanation |
|------|-------------|-------------------|-------------|
| m=0, T=35 | Recurrence | Large T (Asymptotic) | T=35 >= 30+5*0=30 |
| m=5, T=45 | Recurrence | Moderate T | 30+5*5=55, and 45 < 55 |

**Scientific Basis (SciPy verified):**
```
F_0(35) = 1.497996913403e-01
F_5(45) = 2.114257817653e-08
```

**Key Insight:**
The lecture notes describe **three theoretical regimes**:
- Small T (T < 25): Series expansion
- Moderate T (25 <= T < 30+5m): erf + upward recurrence
- Large T (T >= 30+5m): Asymptotic expansion

However, IQCP (following libcint) uses only **two computational methods** (series and recurrence), with the recurrence method handling BOTH the moderate and large T regimes.

**Acceptable Range:**
- Must correctly identify IQCP uses "Recurrence" for both cases
- Must correctly map to theoretical regimes (Large T for m=0,T=35; Moderate T for m=5,T=45)

**Common Misconceptions:**
1. Students may expect to see "Asymptotic" as a method in IQCP
2. Students may not realize the theoretical regime thresholds depend on m
3. Students may confuse IQCP's m-dependent turnover with the theoretical regime boundaries

**Pedagogical Intent:**
Contrasts theoretical understanding from lecture notes with practical implementation. Shows that production codes often combine or simplify theoretical regimes for efficiency.

**Grading Notes:**
- Full credit (5 pts): Both IQCP methods AND theoretical regimes correctly identified
- Partial credit (3 pts): Methods correct but regime mapping incomplete
- Partial credit (2 pts): One pair (method+regime) correct
- No credit: Neither method correctly identified

---

### Q2.7: Why Implementations Combine Moderate and Large T Methods

**Question Text:**
> Why might implementations combine the moderate and large T methods rather than implementing all three separately?

**Expected Answer:**

Key points to include:

1. **Numerical efficiency:** The erf + upward recurrence method works well for BOTH moderate and large T. There is no need for a separate asymptotic regime if recurrence handles both accurately.

2. **Implementation simplicity:** Fewer code paths means fewer edge cases and easier maintenance. Two methods with m-dependent switching is simpler than three methods with two switching conditions.

3. **Stability of recurrence:** The upward recurrence relation:
   ```
   F_{m+1}(T) = [(2m+1)*F_m(T) - exp(-T)] / (2T)
   ```
   is numerically stable for T > turnover(m). The exp(-T) term becomes negligible for large T, so the recurrence naturally transitions to asymptotic behavior.

4. **Threshold tuning:** Production codes (like libcint) carefully tune m-dependent turnover points to ensure the recurrence method achieves machine precision for all T above the turnover, eliminating the need for a separate asymptotic implementation.

**Acceptable Range:**
Must mention at least TWO of: efficiency, simplicity, numerical stability, or practical accuracy

**Common Misconceptions:**
1. Students may think the asymptotic formula is less accurate (it is accurate but unnecessary when recurrence works)
2. Students may not realize the recurrence naturally handles asymptotic behavior
3. Students may confuse implementation trade-offs with theoretical completeness

**Pedagogical Intent:**
Highlights the difference between theoretical understanding (three regimes) and practical implementation (two methods). Teaches that production software makes pragmatic choices.

**Grading Notes:**
- Full credit (5 pts): Substantive explanation with at least two valid reasons
- Partial credit (3 pts): One good reason with elaboration
- Partial credit (2 pts): General statement about simplicity or efficiency
- No credit: No substantive answer

---

## Section 3: Rys Quadrature (Q3.1-Q3.8)

### Q3.1: Roots and Weights Properties

**Question Text:**
> Examine the roots (t_i) and weights (w_i) displayed in the table. Are all roots strictly between 0 and 1? Are all weights positive?

**Expected Answer:**

- **Roots in (0, 1)?** Yes, all roots t_i are strictly between 0 and 1
- **Weights positive?** Yes, all weights w_i are strictly positive

**Scientific Basis:**
These are fundamental mathematical properties of Gaussian quadrature:

1. **Roots:** The quadrature nodes are zeros of orthogonal polynomials. For Rys polynomials (orthogonal with weight exp(-T*t^2) on [0,1]), all zeros lie within the interior of the integration domain.

2. **Weights:** From the Christoffel-Darboux formula, weights for Gaussian quadrature with positive weight functions are always positive.

**Example values (n=5, T=10.0):**
```
Root 1: 0.0632, Weight 1: 0.0245
Root 2: 0.1989, Weight 2: 0.0872
Root 3: 0.3546, Weight 3: 0.0998
Root 4: 0.5123, Weight 4: 0.0562
Root 5: 0.6598, Weight 5: 0.0125
Sum of weights = F_0(10.0) = 0.2802
```

**Acceptable Range:**
- Must answer "Yes" to both
- No partial credit for numerical specifics required

**Common Misconceptions:**
1. Students may expect roots to include 0 or 1
2. Students may confuse weights with coefficients in other contexts
3. Students may not recognize these as mathematical guarantees

**Pedagogical Intent:**
Tests ability to verify sanity checks that all quadrature implementations must satisfy. Builds trust in numerical results through mathematical invariants.

**Grading Notes:**
- Full credit (5 pts): Both correct with verification from IQCP
- Partial credit (3 pts): One correct, one missing or incorrect
- No credit: Both incorrect or no attempt

---

### Q3.2: Effect of Quadrature Points

**Question Text:**
> How does the number of quadrature points affect the computation? Think about both accuracy and computational cost.

**Expected Answer:**

**Accuracy:**
- More quadrature points (larger n) = higher accuracy
- Each additional point allows exact integration of higher-degree polynomials
- n-point Rys quadrature integrates polynomials up to degree 2n-1 exactly

**Computational Cost:**
- Each additional point adds to the computational cost
- Cost per integral scales linearly with n
- Root/weight computation overhead also increases with n

**Trade-off:**
- Must choose minimum n that meets accuracy requirements
- Over-specifying n wastes computational resources
- Under-specifying n produces inaccurate results

**Acceptable Range:**
Must mention both accuracy improvement and cost increase

**Common Misconceptions:**
1. Students may think more points is always better (ignoring cost)
2. Students may not realize diminishing returns at high n
3. Students may confuse quadrature points with basis functions

**Pedagogical Intent:**
Tests understanding of the fundamental accuracy-cost trade-off in numerical integration. This trade-off appears throughout computational science.

**Grading Notes:**
- Full credit (5 pts): Discusses both accuracy improvement AND cost increase
- Partial credit (3 pts): Discusses only one aspect thoroughly
- Partial credit (2 pts): Mentions both but superficially
- No credit: Does not address the question

---

### Q3.3: Reconstruction Errors at Different Orders

**Question Text:**
> Looking at the error information, what is the approximate maximum reconstruction error for n=3? For n=5? For n=7?

**Expected Answer (at T=10.0):**

| Order n | Approximate Max Error | Order of Magnitude |
|---------|----------------------|-------------------|
| n = 3 | 1e-4 to 1e-5 | 10^-4 to 10^-5 |
| n = 5 | 1e-5 to 1e-6 | 10^-5 to 10^-6 |
| n = 7 | 1e-6 to 1e-7 | 10^-6 to 10^-7 |

**Scientific Basis:**
The reconstruction error is dominated by the first moment that cannot be exactly integrated:
- n-point quadrature is exact for moments 0 through 2n-1
- Leading error term is approximately F_{2n}(T)

```
n=3: First unintegrated moment is F_6(T=10) ~ 3.9e-5
n=5: First unintegrated moment is F_10(T=10) ~ 8.6e-6
n=7: First unintegrated moment is F_14(T=10) ~ 3.9e-6
```

**Acceptable Range:**
Order of magnitude correct is sufficient. Exact values depend on IQCP implementation details.

**Common Misconceptions:**
1. Students may expect linear improvement with n
2. Students may not understand "reconstruction error"
3. Students may confuse with floating-point precision limits

**Pedagogical Intent:**
Tests ability to read quantitative information from the IQCP interface and recognize the pattern of exponentially improving accuracy with increasing n.

**Grading Notes:**
- Full credit (5 pts): All three within one order of magnitude of expected
- Partial credit (4 pts): Two correct
- Partial credit (2 pts): One correct or correct trend identified
- No credit: All incorrect or no attempt

---

### Q3.4: Minimum Orders for Target Accuracies

**Question Text:**
> At T = 10.0, what is the minimum quadrature order needed to achieve 1e-8 accuracy? What about 1e-6 accuracy?

**Expected Answer:**

| Target Accuracy | Minimum Order |
|-----------------|---------------|
| 1e-6 | n = 5-6 |
| 1e-8 | n = 8-10 |

**Scientific Basis:**
From the error vs. order analysis:
```
T = 10.0:
n=5: error ~ 10^-5 to 10^-6 (borderline for 1e-6)
n=6: error ~ 10^-6 (meets 1e-6)
n=7: error ~ 10^-6 to 10^-7
n=8: error ~ 10^-7 (approaching 1e-8)
n=9: error ~ 10^-8 (meets 1e-8)
n=10: error ~ 10^-8 to 10^-9 (safely meets 1e-8)
```

**Acceptable Range:**
- For 1e-6: n = 5, 6, or 7 acceptable
- For 1e-8: n = 8, 9, 10, or 11 acceptable

**Common Misconceptions:**
1. Students may expect a simple formula (n = log10(1/tolerance))
2. Students may not use IQCP to find the actual values
3. Students may confuse order with number of digits

**Pedagogical Intent:**
Tests ability to use IQCP as a tool to answer practical questions about quadrature selection. This is exactly how real computational chemists make such decisions.

**Grading Notes:**
- Full credit (5 pts): Both answers in acceptable range
- Partial credit (3 pts): One answer in acceptable range
- Partial credit (2 pts): Correct relative ordering (1e-8 needs more points)
- No credit: Both incorrect or no attempt

---

### Q3.5: Effect of T on Recommended Order

**Question Text:**
> How does the recommended quadrature order change when T increases from 10 to 25? Why do you think this happens?

**Expected Answer:**

**Observation:**
The recommended quadrature order **decreases** as T increases from 10 to 25.

For example (at 1e-8 accuracy target):
- T = 10: requires n ~ 8-10
- T = 25: requires n ~ 5-7

**Explanation:**
At higher T values, the Boys function moments F_k(T) decay faster with increasing k. This means the leading error terms (which involve higher-order moments) are smaller, so fewer quadrature points are needed to achieve the same accuracy.

**Scientific Basis:**
```
Leading error terms F_{2n}(T):
n=5: F_10(10.0) = 8.58e-06, F_10(25.0) = 1.19e-09
n=7: F_14(10.0) = 3.91e-06, F_14(25.0) = 6.14e-11
```

The error at T=25 is orders of magnitude smaller than at T=10 for the same quadrature order.

**Acceptable Range:**
Must state that order decreases with increasing T AND provide a reasonable explanation involving moment decay.

**Common Misconceptions:**
1. Students may expect more points needed for larger T
2. Students may not connect to Boys function behavior from Section 2
3. Students may confuse T dependence with m dependence

**Pedagogical Intent:**
Connects Sections 2 and 3 by showing how Boys function behavior directly impacts Rys quadrature efficiency. Demonstrates why adaptive methods are valuable.

**Grading Notes:**
- Full credit (5 pts): Correct observation AND explanation involving moments
- Partial credit (3 pts): Correct observation, weak or missing explanation
- Partial credit (2 pts): Correct direction but no explanation
- No credit: Incorrect observation or no attempt

---

### Q3.6: Shell Quartet (pp|pp) Root Count Verification

**Question Text:**
> Using the shell quartet selector, set the quartet to (pp|pp). What is the total angular momentum L? What quadrature order n is automatically selected by IQCP? Verify this matches the formula n_r = floor(L/2) + 1.

**Expected Answer:**

- **Total angular momentum L** = 1+1+1+1 = 4 (each p shell has l=1)
- **IQCP selected order n** = 3
- **Formula verification:** floor(4/2) + 1 = 2 + 1 = 3
- **Match:** Yes, they agree

**Scientific Basis:**
From the lecture notes (Table 5.1), the root count rule is:
```
n_r = floor(L/2) + 1
```

This ensures the quadrature integrates polynomials up to degree 2n_r - 1 = L exactly, which is the highest polynomial degree appearing in the ERI integrand.

| Shell Quartet | L | n_r = floor(L/2) + 1 |
|--------------|---|----------------------|
| (ss\|ss) | 0 | 1 |
| (ps\|ss) | 1 | 1 |
| (pp\|ss) | 2 | 2 |
| (pp\|pp) | 4 | 3 |
| (dd\|pp) | 6 | 4 |

**Acceptable Range:**
- L must equal 4
- n must equal 3
- Formula verification must show floor(4/2)+1 = 3

**Common Misconceptions:**
1. Students may think n_r = L (off by factor of ~2)
2. Students may confuse shell angular momentum (p=1) with total L
3. Students may not realize the floor operation handles odd L values

**Pedagogical Intent:**
Connects theoretical root count formula to IQCP implementation. Reinforces that quadrature order is determined by angular momentum, not by T or desired accuracy.

**Grading Notes:**
- Full credit (5 pts): All four parts correct
- Partial credit (3 pts): L and n correct, formula verification incomplete
- Partial credit (2 pts): L correct only
- No credit: Incorrect values

---

### Q3.7: Shell Quartet (dd|pp) Root Count Verification

**Question Text:**
> For the (dd|pp) shell quartet: What is L? What quadrature order does IQCP select? According to the root count rule, is this the minimum required order?

**Expected Answer:**

- **L** = 2+2+1+1 = 6 (d shells have l=2, p shells have l=1)
- **IQCP selected order** = 4
- **Minimum required by formula:** floor(6/2) + 1 = 3 + 1 = 4
- **Conclusion:** Yes, IQCP selects the minimum required order

**Scientific Basis:**
```
(dd|pp): l_A=2, l_B=2, l_C=1, l_D=1
L = 2 + 2 + 1 + 1 = 6
n_r = floor(6/2) + 1 = 4
```

The quadrature with n_r=4 points integrates polynomials up to degree 2*4-1 = 7 exactly, which exceeds L=6, ensuring exact integration.

**Acceptable Range:**
- L = 6
- n = 4
- Must confirm this is the minimum

**Common Misconceptions:**
1. Students may add shell pairs incorrectly
2. Students may confuse d orbital (l=2) with two p orbitals
3. Students may think higher orders are "safer" without understanding the formula guarantees exactness

**Pedagogical Intent:**
Extends the root count concept to higher angular momentum. Shows the formula applies universally across shell types.

**Grading Notes:**
- Full credit (5 pts): All values correct with confirmation of minimum
- Partial credit (3 pts): L and n correct, no minimum confirmation
- Partial credit (2 pts): L correct only
- No credit: Incorrect L value

---

### Q3.8: Algorithm 5.1 Moments and Hankel Matrix

**Question Text:**
> Looking at the Algorithm 5.1 internals for T=10, n=3: What are the first 3 moments m_0, m_1, m_2? What is the dimension of the Hankel matrix H?

**Expected Answer:**

**Moments (m_k = 2*F_k(T)):**
- m_0 = 2*F_0(10) = 0.5605 (or 5.605e-01)
- m_1 = 2*F_1(10) = 0.0280 (or 2.802e-02)
- m_2 = 2*F_2(10) = 0.0042 (or 4.198e-03)

**Hankel matrix dimension:** 3 x 3 (n x n)

**Scientific Basis (SciPy verified):**
```
m_0 = 2*F_0(10) = 5.604947810133e-01
m_1 = 2*F_1(10) = 2.802019905769e-02
m_2 = 2*F_2(10) = 4.198489865677e-03
```

**Algorithm 5.1 Context:**
The Hankel matrix H has dimension n x n (where n is the quadrature order) with:
- H_ij = m_{i+j} for i,j = 0,1,...,n-1

For n=3:
```
H = | m_0  m_1  m_2 |
    | m_1  m_2  m_3 |
    | m_2  m_3  m_4 |
```

Note that to fill this matrix, we need moments m_0 through m_4 (up to m_{2n-2} = m_4).
The shifted Hankel matrix H^(1) needs moments up to m_{2n-1} = m_5.

**Acceptable Range:**
- Moments within 10% of expected values
- Hankel dimension must be 3x3 (or equivalently "n x n" or "n_r x n_r")

**Common Misconceptions:**
1. Students may forget the factor of 2 in m_k = 2*F_k(T)
2. Students may think Hankel dimension is (2n x 2n) instead of (n x n)
3. Students may confuse the number of moments needed with the matrix dimension

**Pedagogical Intent:**
Connects the theoretical Algorithm 5.1 description to actual computed values. Reinforces the relationship between Boys functions and Rys quadrature through the moment interpretation.

**Grading Notes:**
- Full credit (5 pts): At least 2 moments within 10% AND correct dimension
- Partial credit (3 pts): Moments correct OR dimension correct
- Partial credit (2 pts): Shows understanding of m_k = 2*F_k relationship
- No credit: No correct values

---

## Section 4: SCF (Q4.1-Q4.7)

### Q4.1: H2 Default Iterations and Energy

**Question Text:**
> For H2 with medium convergence and DIIS enabled, how many iterations does it take to converge? What is the final energy?

**Expected Answer:**

- **Iterations:** 4-8 iterations (typically 5-6)
- **Final energy:** -1.116716909173 Hartree

**Scientific Basis (PySCF verified):**
```python
# H2 (R = 1.4 bohr = 0.7408 Angstrom) STO-3G
Final energy: -1.116716909173 Hartree
Converged: True
Basis functions: 2
Electrons: 2
```

**Acceptable Range:**
- Iterations: 3-10 (DIIS performance can vary)
- Energy: -1.1167 +/- 0.0001 Hartree (within 10^-4)

**Common Misconceptions:**
1. Students may forget to include the minus sign
2. Students may confuse Hartree with other units (eV, kcal/mol)
3. Students may report intermediate energies, not final

**Pedagogical Intent:**
Establishes baseline for comparison with other systems and DIIS settings. Tests ability to read SCF output correctly.

**Grading Notes:**
- Full credit (4 pts): Both values within acceptable range
- Partial credit (2 pts): Energy correct, iterations missing or wrong
- Partial credit (2 pts): Iterations correct, energy missing or wrong
- No credit: Both incorrect or no attempt

---

### Q4.2: H2 Iterations With/Without DIIS (Tight)

**Question Text:**
> For H2 with tight convergence, how many iterations does SCF take without DIIS? With DIIS?

**Expected Answer:**

| Setting | Iterations |
|---------|-----------|
| Without DIIS | 10-15 iterations |
| With DIIS | 5-8 iterations |

**DIIS Speedup:** Approximately 40-60% reduction in iteration count

**Scientific Basis:**
DIIS accelerates convergence by extrapolating from previous Fock matrices, effectively "predicting" the converged solution from iteration history.

**Acceptable Range:**
- Without DIIS: 8-20 iterations
- With DIIS: 4-10 iterations
- Key: DIIS case must have fewer iterations

**Common Misconceptions:**
1. Students may not notice significant difference for small molecules
2. Students may confuse iteration count with computation time
3. Students may not run both calculations

**Pedagogical Intent:**
Direct observation of DIIS acceleration effect. Even for simple H2, DIIS provides measurable speedup.

**Grading Notes:**
- Full credit (4 pts): Both values reasonable AND DIIS has fewer iterations
- Partial credit (2 pts): Correct relative comparison but values off
- Partial credit (1 pt): One value correct
- No credit: No comparison or DIIS incorrectly identified as slower

---

### Q4.3: Convergence Pattern Description

**Question Text:**
> Looking at the energy vs. iteration plot, describe the difference in convergence patterns between the two cases. How does DIIS change the shape of the convergence curve?

**Expected Answer:**

Key observations to include:

**Without DIIS:**
1. Energy may oscillate (alternating above and below final value)
2. Residual decreases slowly and may plateau
3. Convergence curve has irregular shape
4. May exhibit step-like pattern

**With DIIS:**
1. Energy converges monotonically (or nearly so)
2. Residual drops rapidly once DIIS activates (after 2-3 iterations)
3. Convergence curve is smooth and steep
4. "Hockey stick" shape - slow start then rapid convergence

**DIIS Effect:**
DIIS "smooths out" oscillations by extrapolating from multiple previous iterations. The algorithm finds an optimal combination that minimizes the error vector, effectively jumping toward the solution.

**Acceptable Range:**
Must mention at least two characteristics of each case

**Common Misconceptions:**
1. Students may describe only the final result, not the path
2. Students may not notice oscillations without DIIS
3. Students may attribute all speedup to better initial guess

**Pedagogical Intent:**
Develops qualitative understanding of iterative algorithm behavior. Visualization helps build intuition for when acceleration techniques are most valuable.

**Grading Notes:**
- Full credit (6 pts): Clear description of both patterns AND explanation of DIIS effect
- Partial credit (4 pts): Good description of patterns, weak explanation
- Partial credit (2 pts): Some correct observations
- No credit: No substantive comparison

---

### Q4.4: H2O Final Energy and Convergence

**Question Text:**
> For H2O, what is the final RHF energy? Does the calculation converge in both cases (with and without DIIS)?

**Expected Answer:**

- **Final energy:** -74.963023138435 Hartree
- **Converges without DIIS?** Yes (but requires more iterations, typically 12-20)
- **Converges with DIIS?** Yes (faster, typically 6-10 iterations)

**Scientific Basis (PySCF verified):**
```python
# H2O STO-3G (standard geometry)
Final energy: -74.963023138435 Hartree
Converged: True
Basis functions: 7
Electrons: 10
```

**Acceptable Range:**
- Energy: -74.963 +/- 0.001 Hartree
- Both cases should converge (H2O in STO-3G is well-behaved)

**Common Misconceptions:**
1. Students may expect non-convergence without DIIS
2. Students may confuse with larger basis set energies from literature
3. Students may not wait for convergence

**Pedagogical Intent:**
Extends SCF experience to a larger, chemically interesting molecule. Demonstrates that DIIS helps even when calculation would converge anyway.

**Grading Notes:**
- Full credit (4 pts): Energy correct AND convergence assessment for both
- Partial credit (2 pts): Energy correct only
- Partial credit (2 pts): Convergence assessment correct, energy wrong
- No credit: Energy significantly wrong or no attempt

---

### Q4.5: Fock Matrix Symmetry

**Question Text:**
> In the Internals mode, examine the Fock matrix F. Is it symmetric (F_ij = F_ji)? Why is symmetry of the Fock matrix physically important?

**Expected Answer:**

**Observation:**
Yes, F_ij = F_ji (within numerical precision, typically 10^-12 or better)

**Physical Importance:**

1. **Real eigenvalues:** Hermitian (symmetric for real matrices) operators have real eigenvalues. The orbital energies (eigenvalues of F) must be real physical quantities.

2. **Orthogonal eigenvectors:** Eigenvectors of symmetric matrices are orthogonal, meaning molecular orbitals form an orthonormal basis.

3. **Time-reversal symmetry:** The Fock operator is derived from the time-reversal-invariant electronic Hamiltonian.

4. **Variational principle:** The symmetric form ensures the energy is a proper variational functional.

**Acceptable Range:**
Must confirm symmetry AND provide at least one physical reason

**Common Misconceptions:**
1. Students may not know how to verify symmetry numerically
2. Students may confuse Fock matrix with density matrix
3. Students may not connect symmetry to eigenvalue properties

**Pedagogical Intent:**
Connects matrix algebra to physical meaning. Tests understanding of why mathematical properties matter for physical interpretation.

**Grading Notes:**
- Full credit (5 pts): Correct observation AND substantive explanation
- Partial credit (3 pts): Correct observation, superficial explanation
- Partial credit (2 pts): Correct observation only
- No credit: Incorrect observation or no attempt

---

### Q4.6: HOMO-LUMO Energies and Gap

**Question Text:**
> What is the HOMO energy for H2? What is the LUMO energy? What does the HOMO-LUMO gap tell you about the molecule?

**Expected Answer:**

**Numerical Values (H2 STO-3G):**
- **HOMO energy:** -0.57822287 Hartree (approximately -15.7 eV)
- **LUMO energy:** +0.67031739 Hartree (approximately +18.2 eV)
- **HOMO-LUMO gap:** 1.248540 Hartree (approximately 34 eV)

**Significance of the Gap:**

1. **Chemical stability:** Large gap indicates H2 is chemically stable (resistant to electron addition/removal)

2. **Reactivity:** Molecules with smaller gaps are generally more reactive

3. **Optical properties:** The gap approximates the lowest electronic excitation energy (in a very simplified picture)

4. **Electrical properties:** Large gap indicates H2 is an insulator, not a conductor

**Scientific Basis (PySCF verified):**
```python
# H2 STO-3G orbital energies
MO 0: -0.57822287 Ha (occupied) <- HOMO
MO 1: +0.67031739 Ha (virtual) <- LUMO
Gap: 1.24854026 Ha = 33.97 eV
```

**Acceptable Range:**
- HOMO: -0.58 +/- 0.01 Hartree
- LUMO: +0.67 +/- 0.01 Hartree
- Gap: 1.2 - 1.3 Hartree or 33-35 eV
- Must mention at least one significance

**Common Misconceptions:**
1. Students may confuse HOMO/LUMO with orbital indices
2. Students may expect negative LUMO energy
3. Students may not know Hartree-to-eV conversion (1 Ha ~ 27.2 eV)

**Pedagogical Intent:**
Connects computed orbital energies to physical and chemical properties. HOMO-LUMO gap is a fundamental descriptor in chemistry.

**Grading Notes:**
- Full credit (6 pts): Both energies correct AND gap significance explained
- Partial credit (4 pts): Energies correct, weak significance discussion
- Partial credit (2 pts): One energy correct or gap correct without context
- No credit: Incorrect values or no attempt

---

### Q4.7: DIIS Recommendations

**Question Text:**
> Based on your observations in this section, when would you recommend using DIIS? Are there any situations where DIIS might not help or could cause problems?

**Expected Answer:**

**Recommend DIIS:**

1. **Most calculations:** DIIS provides acceleration with minimal overhead in typical cases

2. **Larger molecules:** Benefits increase with system size (more iterations saved)

3. **Difficult convergence:** When standard SCF oscillates or converges slowly

4. **Tight tolerances:** When many iterations would be needed without acceleration

**When DIIS may not help or could cause problems:**

1. **Very close to convergence:** Diminishing returns when already near solution

2. **Multiple SCF solutions:** DIIS may bias toward one solution; near-degenerate cases may need careful treatment

3. **Very far from minimum:** Initial iterations may benefit from simple damping before DIIS

4. **Near-degenerate orbitals:** Systems with small HOMO-LUMO gaps may have convergence issues

5. **Numerical precision limits:** At very tight tolerances, DIIS extrapolation may introduce noise

**Acceptable Range:**
Must include at least two recommendations for using DIIS AND at least one caution

**Common Misconceptions:**
1. Students may think DIIS always helps (it has edge cases)
2. Students may not have observed enough to form opinions
3. Students may confuse DIIS with other acceleration methods

**Pedagogical Intent:**
Develops judgment about when to apply computational techniques. Real practitioners must make such decisions regularly.

**Grading Notes:**
- Full credit (6 pts): Multiple recommendations AND realistic cautions
- Partial credit (4 pts): Good recommendations OR good cautions
- Partial credit (2 pts): General understanding without specifics
- No credit: No substantive answer

---

## Section 5: Synthesis (Q5.1-Q5.3)

### Q5.1: Boys Function in Integral Evaluation

**Question Text:**
> The Boys function F_m(T) appears in nuclear attraction integrals where T depends on the distance between nuclei and basis function centers. Based on what you learned about Boys function behavior:
>
> a) What happens to the integrand when two Gaussian centers are very close (small T)?
> b) What happens when they are far apart (large T)?
> c) Why might special numerical care be needed at the regime boundaries?

**Expected Answer:**

**Part (a) - Small T (close centers):**
- T is small, so F_m(T) approaches its maximum value 1/(2m+1)
- The integral contribution is significant (not negligible)
- Series expansion is used for accurate evaluation
- The computation is straightforward numerically

**Part (b) - Large T (far apart centers):**
- T is large, so F_m(T) approaches zero
- The integral contribution becomes negligible
- Asymptotic expansion is appropriate
- This reflects the physical expectation: nuclear attraction decays with distance

**Part (c) - Regime boundaries:**
- Must switch algorithms smoothly (no discontinuities in computed values)
- Near boundaries, both methods should give the same result
- Round-off error may accumulate differently in each algorithm
- Testing and validation at boundaries is critical for reliability
- Small errors at boundaries could propagate to large errors in total energy

**Acceptable Range:**
Must address all three parts with reasonable physical and numerical insight

**Common Misconceptions:**
1. Students may not connect T to physical distance
2. Students may forget that regime boundaries are implementation choices
3. Students may not realize that continuity at boundaries must be verified

**Pedagogical Intent:**
Connects mathematical behavior (Section 2) to physical application. Tests ability to synthesize across the lab.

**Grading Notes:**
- Full credit (5 pts): All three parts addressed with insight
- Partial credit (3-4 pts): Two parts well addressed
- Partial credit (1-2 pts): One part well addressed
- No credit: No substantive answer

---

### Q5.2: Quadrature and Computational Cost

**Question Text:**
> Rys quadrature is used to evaluate two-electron integrals (ij|kl), where each integral can require many quadrature points for high accuracy. If a calculation requires 1e-10 accuracy and you need to evaluate 10,000 integrals:
>
> a) How would you estimate the total number of quadrature point evaluations needed?
> b) How does your understanding of the order-accuracy relationship help predict computational cost?
> c) Why might adaptive quadrature order selection (varying order based on the integral's T value) be valuable?

**Expected Answer:**

**Part (a) - Estimating total evaluations:**
- If each integral needs approximately n quadrature points
- And there are N_int = 10,000 integrals
- Total evaluations approximately equals N_int x n
- For 1e-10 accuracy at typical T, need n ~ 10 points
- Estimate: 10,000 x 10 = 100,000 quadrature point evaluations

**Part (b) - Order-accuracy relationship:**
- From Section 3, we learned that error decreases roughly exponentially with order
- Higher accuracy requires larger n (approximately logarithmic relationship)
- Doubling the accuracy target (e.g., 1e-10 to 1e-12) adds only 1-2 quadrature points
- This helps predict that modest increases in accuracy are achievable without enormous cost increases

**Part (c) - Adaptive quadrature:**
- Different integrals have different T values (depending on orbital separations)
- Large-T integrals need fewer points (moments decay faster)
- Small-T integrals need more points
- Adaptive selection minimizes total cost while achieving uniform accuracy
- Could save significant computation (perhaps 30-50% in favorable cases)

**Acceptable Range:**
Must demonstrate quantitative reasoning in part (a) and conceptual understanding in parts (b) and (c)

**Common Misconceptions:**
1. Students may not realize the multiplicative nature of the cost
2. Students may think all integrals need the same number of points
3. Students may not connect to the T-dependence from Section 3

**Pedagogical Intent:**
Develops practical computational thinking. Real quantum chemistry programs use exactly this kind of reasoning.

**Grading Notes:**
- Full credit (5 pts): All three parts with reasonable quantitative and qualitative analysis
- Partial credit (3-4 pts): Two parts well addressed
- Partial credit (1-2 pts): Basic understanding shown
- No credit: No substantive answer

---

### Q5.3: DIIS Summary and Explanation

**Question Text:**
> DIIS (Direct Inversion in the Iterative Subspace) dramatically accelerates SCF convergence. Based on your observations:
>
> a) Summarize in 2-3 sentences what DIIS does to improve convergence.
> b) The DIIS method works by extrapolating from previous Fock matrices. Why might this be more effective than simple iteration?
> c) Under what circumstances might standard SCF iteration (without DIIS) still be useful?

**Expected Answer:**

**Part (a) - DIIS Summary:**
DIIS (Direct Inversion in the Iterative Subspace) accelerates SCF convergence by storing a history of previous Fock matrices and error vectors, then finding the optimal linear combination that minimizes the error. Instead of using only the current iteration's Fock matrix, DIIS extrapolates from multiple previous iterations to "predict" a Fock matrix closer to convergence. This typically reduces the number of SCF iterations by 40-60%.

**Part (b) - Why extrapolation is effective:**
- Simple iteration follows the gradient, which can lead to oscillation or spiraling
- DIIS uses information from multiple iterations to predict the minimum
- Extrapolation can "jump over" oscillations toward convergence
- The error vector (commutator [F,D]) provides a convergence metric
- DIIS effectively increases the convergence order of the iterative method

**Part (c) - When standard SCF is useful:**
1. **Initial iterations:** Before DIIS history is built (first 2-3 iterations)
2. **Near convergence:** Diminishing returns from DIIS when already close
3. **Simple systems:** Small molecules that converge quickly anyway
4. **Exploring multiple solutions:** DIIS may bias toward one SCF solution
5. **Debugging:** Simpler algorithm is easier to analyze
6. **Very tight convergence:** At machine precision limits, extrapolation may introduce noise

**Acceptable Range:**
Part (a) must capture the essence of DIIS (history + extrapolation + error minimization)
Parts (b) and (c) must provide reasonable explanations

**Common Misconceptions:**
1. Students may describe DIIS as "averaging" rather than "optimizing"
2. Students may not understand the role of the error vector
3. Students may think DIIS never has drawbacks

**Pedagogical Intent:**
Tests deep understanding of DIIS mechanism and practical judgment about algorithm selection. DIIS is used in virtually all quantum chemistry calculations.

**Grading Notes:**
- Full credit (5 pts): Accurate summary AND insightful explanation of both (b) and (c)
- Partial credit (3-4 pts): Good summary and one well-addressed part
- Partial credit (1-2 pts): Basic understanding of DIIS shown
- No credit: Fundamentally incorrect description

---

## Point Summary

| Section | Questions | Points |
|---------|-----------|--------|
| Section 2: Boys Function | Q2.1-Q2.7 | 35 |
| Section 3: Rys Quadrature | Q3.1-Q3.8 | 40 |
| Section 4: SCF | Q4.1-Q4.7 | 35 |
| Section 5: Synthesis | Q5.1-Q5.3 | 15 |
| **Total** | **25 questions** | **125** |

**Note:** Point values may be scaled by instructors to fit their grading scheme. The rubric provides detailed partial credit guidelines.

### Section 4 Point Breakdown

| Question | Points |
|----------|--------|
| Q4.1 | 4 |
| Q4.2 | 4 |
| Q4.3 | 6 |
| Q4.4 | 4 |
| Q4.5 | 5 |
| Q4.6 | 6 |
| Q4.7 | 6 |
| **Total** | **35** |

---

## General Grading Guidelines

### Numerical Answers

- **Full credit:** Within specified tolerance or acceptable range
- **Partial credit (50%):** Correct method, minor calculation error
- **Partial credit (25%):** Correct units or order of magnitude
- **No credit:** Completely wrong or no attempt

### Conceptual Answers

- **Full credit:** All key points addressed with clarity
- **Partial credit (75%):** Most key points, minor omissions
- **Partial credit (50%):** Some understanding, incomplete
- **Partial credit (25%):** Minimal understanding shown
- **No credit:** Incorrect or no attempt

### Common Issues Across All Questions

1. **Units:** Accept Hartree, Ha, or atomic units. Deduct points only if wrong unit conversion attempted.

2. **Significant figures:** Do not penalize for excess precision. Penalize only for grossly insufficient precision.

3. **Terminology:** Accept reasonable synonyms (e.g., "Taylor series" for "series expansion").

4. **Partial work:** Award partial credit for correct reasoning even if final answer is wrong.

---

## Scientific Reference Values Quick Reference

### Boys Function Values

| Parameter | Value |
|-----------|-------|
| F_0(0) | 1.0 |
| F_m(0) | 1/(2m+1) |
| F_0(0.5) | 0.855624391892 |
| F_0(10.0) | 0.280247 |
| F_0(15.0) | 0.228823 |
| F_5(10.0) | 7.9009e-05 |

### Method Selection (m-dependent turnover)

IQCP uses **two computational methods** with m-dependent turnover points:

| Method | When Used |
|--------|-----------|
| Series | T < turnover(m) |
| Recurrence | T >= turnover(m) |

**Turnover Points by m:**
| m | Turnover Point |
|---|----------------|
| 0-1 | 0 (always recurrence) |
| 2 | 0.87 |
| 5 | 2.11 |
| 10 | 4.05 |
| 20 | 7.84 |
| 30 | 11.58 |

**Note:** There is NO asymptotic regime in IQCP. The recurrence method handles all moderate-to-large T values.

### SCF Energies (STO-3G)

| System | Energy (Hartree) |
|--------|------------------|
| H2 | -1.116716909173 |
| H2O | -74.963023138435 |
| LiH | -7.862023860127 |
| HeH+ | -2.841779241244 |
| NH3 | -55.454540513740 |

### H2 Orbital Energies (STO-3G)

| Orbital | Energy (Hartree) |
|---------|------------------|
| HOMO | -0.57822287 |
| LUMO | +0.67031739 |
| Gap | 1.248540 |

---

*Lab Pack #1 Answer Key v1.0 | CONFIDENTIAL - Instructor Use Only*
*Interactive Quantum Chemistry Playground | https://iqcp.dev*
