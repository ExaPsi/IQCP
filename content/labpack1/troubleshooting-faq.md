# Lab Pack #1: Troubleshooting FAQ for Instructors

**Lab Pack:** 1 - From Boys to Orbitals
**Version:** 1.0
**Last Updated:** 2026-01-18
**Document Type:** Instructor Support Materials

---

## Overview

This document provides solutions to common technical issues students may encounter while completing Lab Pack #1. Each FAQ entry includes the problem description, possible causes, step-by-step solutions, and guidance on when to escalate.

**Quick tip:** Most issues can be resolved by refreshing the page or using a different browser. Encourage students to try these steps first.

---

## Quick Reference Table

| Issue | Category | Quick Solution |
|-------|----------|----------------|
| Deep link shows wrong parameters | Deep Links | Clear cache and reload; manually enter parameters |
| Deep link does not open at all | Deep Links | Check for URL truncation; use local URL fallback |
| "Copy Link" button not working | Deep Links | Use browser URL bar; check clipboard permissions |
| WASM module not loading | Browser | Update browser; enable WebAssembly; try Chrome |
| Page is slow or unresponsive | Browser | Close other tabs; disable extensions; hard refresh |
| SCF calculation stuck | Computation | Wait 30 seconds; check max iterations; try Cancel |
| Values differ from answer key | Computation | Normal if within tolerance; check parameters |
| Values show NaN or undefined | Computation | Refresh page; check for invalid inputs |
| Export button missing | Export/Import | Switch to results view; ensure calculation complete |
| Artifact file not importing | Export/Import | Check file format; try different browser |
| Controls not responding | UI/Navigation | Click elsewhere first; check for pending calculation |
| Mode toggle not working | UI/Navigation | Wait for calculation; refresh if stuck |
| Display issues on mobile | UI/Navigation | Use desktop browser; rotate to landscape |

---

## Category 1: Deep Link Issues

### FAQ 1.1: Deep link shows different parameters than described in the worksheet

**Problem:** A student clicks a deep link but IQCP displays different parameter values than what the worksheet step describes.

**Possible Causes:**
1. Browser cached an old version of the application
2. Deep link state was overwritten by previous session
3. URL was modified or partially loaded

**Solution:**

1. **Clear browser cache and reload:**
   - Chrome: Press `Ctrl+Shift+R` (Windows/Linux) or `Cmd+Shift+R` (Mac)
   - Firefox: Press `Ctrl+Shift+R` or `Cmd+Shift+R`
   - Safari: Press `Cmd+Option+R`

2. **Try opening the link in a private/incognito window:**
   - Chrome: `Ctrl+Shift+N` or `Cmd+Shift+N`
   - Firefox: `Ctrl+Shift+P` or `Cmd+Shift+P`
   - This eliminates cached data issues

3. **Manually enter the parameters:**
   - If the link specifies "m=5, T=10.0", use the sliders/inputs in the Controls Panel to set these values
   - The worksheet always includes the target parameter values in the link description

**When to Escalate:** If multiple students experience this issue with the same link, the deep link encoding may need regeneration. Contact the course administrator.

---

### FAQ 1.2: Deep link does not open or shows an error page

**Problem:** Clicking a deep link results in a blank page, 404 error, or the wrong module loading.

**Possible Causes:**
1. URL was truncated when copied (common in email clients)
2. The `?run=` query parameter was cut off
3. Network connectivity issues
4. IQCP server is temporarily unavailable

**Solution:**

1. **Check the URL length:**
   - Deep links contain a `?run=...` parameter with encoded state
   - If the URL ends abruptly (e.g., `?run=N4Ig...` instead of a longer string), it was truncated
   - Copy the full link from the original worksheet PDF/document

2. **Use the local development fallback:**
   - If running IQCP locally: replace `https://iqcp.dev` with `http://localhost:5173`
   - Example: `http://localhost:5173/boys?run=N4Ig...`

3. **Navigate manually:**
   - Go to `https://iqcp.dev` directly
   - Click the appropriate module (Boys, Rys, or SCF)
   - Enter the parameters listed in the worksheet step

4. **Check network connectivity:**
   - Verify other websites load correctly
   - Try accessing `https://iqcp.dev` directly without parameters

**When to Escalate:** If `https://iqcp.dev` is completely unreachable for multiple students, there may be a server outage. Have students use the local version or pause the lab activity.

---

### FAQ 1.3: "Copy Link" button does not work or copy fails

**Problem:** Student clicks the "Copy Link" button but nothing appears in their clipboard, or they receive a permission error.

**Possible Causes:**
1. Browser blocked clipboard access
2. Page was opened from a file:// URL instead of http(s)
3. Browser does not support the Clipboard API
4. Focus was not on the page when clicking

**Solution:**

1. **Grant clipboard permissions:**
   - Chrome: Click the lock icon in the address bar > Site Settings > Clipboard > Allow
   - Firefox: The browser may prompt for permission; click "Allow"
   - Safari: System Preferences > Security & Privacy > Privacy > Accessibility (may need to add browser)

2. **Use the browser URL bar instead:**
   - The current URL in the browser address bar IS the shareable deep link
   - Students can select the URL and copy it directly (`Ctrl+C` / `Cmd+C`)

3. **Ensure page is served over HTTP(S):**
   - Clipboard API requires a secure context
   - Check that the URL starts with `https://` or `http://localhost`

4. **Try clicking inside the page first:**
   - Click anywhere on the IQCP interface
   - Then click the "Copy Link" button

**When to Escalate:** If clipboard functionality fails consistently, students can manually copy URLs from the address bar. This is not a blocking issue.

---

## Category 2: Browser and Compatibility Issues

### FAQ 2.1: WASM module not loading (WebAssembly error)

**Problem:** IQCP shows an error message about WebAssembly, the page displays a loading spinner indefinitely, or calculations do not run.

**Possible Causes:**
1. Browser does not support WebAssembly
2. WebAssembly is disabled in browser settings
3. Browser is severely outdated
4. Security software blocking WASM execution

**Solution:**

1. **Use a supported browser:**
   - Chrome 57+ (recommended)
   - Firefox 52+
   - Safari 11+
   - Edge 16+
   - Check version: enter `chrome://version` (Chrome) or `about:support` (Firefox)

2. **Update browser to latest version:**
   - Chrome: Menu > Help > About Google Chrome
   - Firefox: Menu > Help > About Firefox
   - Edge: Menu > Help and feedback > About Microsoft Edge

3. **Check WebAssembly support:**
   - Open browser DevTools (F12)
   - In Console, type: `typeof WebAssembly === 'object'`
   - Should return `true`

4. **Disable interfering extensions:**
   - Ad blockers or security extensions may block WASM
   - Try incognito/private mode (extensions usually disabled)
   - Or temporarily disable extensions one by one

5. **Check security software:**
   - Some antivirus programs block WebAssembly
   - Temporarily disable or add exception for `iqcp.dev`

**Browser Compatibility Matrix:**

| Browser | Minimum Version | WebAssembly Support |
|---------|-----------------|---------------------|
| Chrome | 57 | Full |
| Firefox | 52 | Full |
| Safari | 11 | Full |
| Edge | 16 | Full |
| Internet Explorer | N/A | Not Supported |

**When to Escalate:** If a student's browser meets requirements but WASM still fails, they may need IT support to check security software or network policies.

---

### FAQ 2.2: Page is slow, unresponsive, or browser becomes sluggish

**Problem:** IQCP runs slowly, slider updates are delayed, or the browser tab becomes unresponsive during calculations.

**Possible Causes:**
1. Insufficient system resources (RAM, CPU)
2. Too many browser tabs open
3. Browser extensions consuming resources
4. Large SCF calculation in progress
5. Browser DevTools open with heavy logging

**Solution:**

1. **Close unnecessary tabs and applications:**
   - IQCP's WASM computations use CPU resources
   - Close other computationally intensive tabs

2. **Perform a hard refresh:**
   - `Ctrl+Shift+R` (Windows/Linux) or `Cmd+Shift+R` (Mac)
   - This clears cached assets and reloads fresh

3. **Disable browser extensions:**
   - Open in private/incognito mode
   - Or disable extensions via browser settings

4. **Close DevTools if open:**
   - Developer tools can slow down rendering
   - Press F12 to close if open

5. **Wait for calculations to complete:**
   - SCF calculations on larger molecules (H2O, NH3) may take 5-15 seconds
   - The UI may be less responsive during computation
   - Look for the "Computing..." indicator

6. **Try a different browser:**
   - Chrome is generally fastest for WASM
   - Firefox is also well-optimized

**Performance Expectations:**

| Module | Typical Response Time |
|--------|----------------------|
| Boys Function | <100ms |
| Rys Quadrature | <200ms |
| SCF (H2) | 1-3 seconds |
| SCF (H2O) | 5-15 seconds |

**When to Escalate:** If performance issues persist on modern hardware with all optimizations, the student's computer may not meet minimum requirements. Consider pairing with another student or providing pre-computed artifacts.

---

## Category 3: Computation Issues

### FAQ 3.1: SCF calculation seems stuck or taking too long

**Problem:** After clicking "Run SCF" or changing parameters, the calculation appears to hang. The progress indicator spins indefinitely.

**Possible Causes:**
1. Calculation is still running (larger systems take longer)
2. SCF is not converging (oscillating or very slow convergence)
3. Max iterations set too low
4. Web Worker has crashed or stalled

**Solution:**

1. **Wait up to 30 seconds:**
   - Larger molecules (H2O, NH3) require more iterations
   - Each SCF iteration involves matrix operations
   - Normal calculation times:
     - H2: 2-5 seconds
     - HeH+: 2-5 seconds
     - H2O: 10-20 seconds

2. **Check the iteration counter:**
   - Look at the SCF iteration table or progress indicator
   - If iterations are incrementing, the calculation is running normally

3. **Click "Cancel" and retry:**
   - If available, click the Cancel button
   - Wait a moment, then click "Run SCF" again

4. **Increase max iterations:**
   - Default is usually 50 iterations
   - If SCF hasn't converged, try 100 or 150
   - Go to Controls Panel > Advanced Settings

5. **Check convergence settings:**
   - Tighter convergence (1e-10 vs 1e-6) requires more iterations
   - For exploration, "medium" convergence (1e-6) is usually sufficient

6. **Refresh the page:**
   - If truly stuck, `Ctrl+R` / `Cmd+R` will reset the application
   - Re-enter parameters and try again

**SCF Convergence Expectations:**

| System | Expected Iterations (no DIIS) | With DIIS |
|--------|------------------------------|-----------|
| H2 | 8-12 | 5-8 |
| HeH+ | 8-12 | 5-8 |
| LiH | 15-25 | 8-12 |
| H2O | 15-25 | 8-12 |

**When to Escalate:** If SCF consistently fails to converge for a preset system, there may be a bug. Have the student try a different preset or report the issue with the exact parameters used.

---

### FAQ 3.2: Numerical values differ from the answer key

**Problem:** Student's calculated values (Boys function, Rys roots, SCF energy) differ from the values in the answer key.

**Possible Causes:**
1. Normal floating-point differences between browsers
2. Parameters not exactly matching worksheet
3. Different convergence settings
4. Viewing intermediate vs. final results

**Solution:**

1. **Check if difference is within tolerance:**
   - Boys function: differences beyond 10th decimal place are normal
   - Rys roots/weights: differences beyond 8th decimal place are normal
   - SCF energy: differences beyond 6th decimal place (microHartree) are normal

   **Tolerance Guidelines:**
   | Quantity | Acceptable Difference |
   |----------|----------------------|
   | Boys F_m(T) | < 1e-10 |
   | Rys roots | < 1e-8 |
   | Rys weights | < 1e-8 |
   | SCF energy | < 1e-6 Ha |

2. **Verify parameters exactly match:**
   - Check m, T values for Boys
   - Check order, T for Rys
   - Check preset system, basis, DIIS settings for SCF

3. **Ensure calculation has completed:**
   - For SCF, check the "Converged" indicator
   - Intermediate energies will differ from final

4. **Check display precision:**
   - Switch between "Explain" and "Internals" modes
   - Internals mode may show more decimal places

**Example of Acceptable Differences:**

| Answer Key | Student Value | Difference | Acceptable? |
|------------|---------------|------------|-------------|
| 0.8556243918921488 | 0.8556243918921490 | 2e-16 | Yes |
| -1.116714349 | -1.116714348 | 1e-9 | Yes |
| -74.96590119 | -74.96589421 | 7e-6 | Marginal* |

*Marginal differences may indicate different convergence criteria.

**When to Escalate:** If values differ by more than the tolerance guidelines, verify the student's exact parameters. Significant discrepancies may indicate a calculation bug.

---

### FAQ 3.3: Values display as NaN, undefined, or Infinity

**Problem:** Instead of numerical results, IQCP displays "NaN", "undefined", "-Infinity", or similar error indicators.

**Possible Causes:**
1. Invalid input parameters (negative T, non-integer m)
2. Numerical overflow or underflow
3. Division by zero in intermediate calculation
4. WASM computation returned an error

**Solution:**

1. **Check input values:**
   - m must be a non-negative integer (0, 1, 2, ...)
   - T must be a non-negative real number (T >= 0)
   - For Rys, order must be a positive integer

2. **Reset to valid parameters:**
   - Use the "Reset" button if available
   - Or navigate to the module fresh: click "Boys", "Rys", or "SCF" in navigation

3. **Avoid extreme parameter values:**
   - Boys: T > 10000 may cause numerical issues
   - Rys: order > 15 may be unstable
   - Keep m < 20 for stability

4. **Refresh the page:**
   - `Ctrl+R` / `Cmd+R`
   - Sometimes WASM state can become corrupted

5. **Clear browser storage:**
   - Chrome: DevTools > Application > Storage > Clear site data
   - This resets any persisted state

**When to Escalate:** If NaN appears with normal parameters (e.g., m=0, T=5.0), this is a bug. Report with exact steps to reproduce.

---

## Category 4: Export and Import Issues

### FAQ 4.1: Export button not appearing or not working

**Problem:** Student cannot find the Export button, or clicking it has no effect.

**Possible Causes:**
1. No calculation results to export yet
2. Export button hidden in collapsed panel
3. Browser blocking file downloads
4. Calculation still in progress

**Solution:**

1. **Ensure a calculation has completed:**
   - The Export button may be disabled until results exist
   - For SCF, wait until "Converged" appears
   - For Boys/Rys, ensure the computed value is displayed

2. **Look in the correct location:**
   - Export button is typically in the Results area or a toolbar
   - May be labeled "Export Artifact", "Download", or show a download icon
   - Check both the Controls Panel and Results Panel

3. **Check browser download settings:**
   - Chrome: Settings > Downloads > Ask where to save
   - Firefox: Settings > General > Downloads
   - Ensure downloads are not blocked

4. **Check for download popup blockers:**
   - Allow popups for iqcp.dev
   - Some ad blockers prevent downloads

5. **Try right-click > Save As on results:**
   - If available, this bypasses the Export button

**When to Escalate:** If Export consistently fails, students can document their results via screenshots instead. This is not ideal but allows completion.

---

### FAQ 4.2: Imported artifact file shows errors or wrong data

**Problem:** After importing a previously exported artifact file, IQCP shows an error, displays wrong data, or does nothing.

**Possible Causes:**
1. File was corrupted during transfer
2. File is from a different IQCP version
3. Wrong file type (not a valid artifact JSON)
4. File encoding issues from email attachment

**Solution:**

1. **Verify file extension and format:**
   - Artifact files should end in `.json`
   - They should be readable JSON (try opening in text editor)
   - First line should contain `{` and file should end with `}`

2. **Check file integrity:**
   - Compare file size to original
   - Very small file (<100 bytes) likely incomplete
   - Re-download or re-receive the file

3. **Use correct Import function:**
   - Navigate to the matching module first (Boys, Rys, SCF)
   - Then use the Import button
   - Importing a Boys artifact in the Rys module will fail

4. **Try re-exporting from original browser:**
   - If student still has access, export again
   - Use direct file save, not copy-paste

5. **Check artifact schema version:**
   - Open JSON file and check `schema_version` field
   - Must match current IQCP version

**Artifact File Structure:**
```json
{
  "schema_version": "run_state_v1",
  "module": "boys",
  "timestamp": "2026-01-18T...",
  "params": { ... },
  "results": { ... }
}
```

**When to Escalate:** If artifact import/export is systematically failing, there may be a schema incompatibility. Students can submit screenshots as a fallback.

---

## Category 5: UI and Navigation Issues

### FAQ 5.1: Controls (sliders, inputs) not responding to changes

**Problem:** Moving sliders or changing input values has no effect on the display. Controls appear frozen.

**Possible Causes:**
1. A calculation is in progress
2. Focus is on another element
3. JavaScript error halted UI
4. Input validation rejecting values silently

**Solution:**

1. **Wait for current calculation to complete:**
   - During SCF iterations, controls may be disabled
   - Look for "Computing..." or spinner indicators

2. **Click elsewhere, then try again:**
   - Click on an empty area of the page
   - Then interact with the control

3. **Check for validation constraints:**
   - Sliders have min/max limits
   - Text inputs may reject invalid characters
   - m must be an integer; T must be a number

4. **Refresh the page:**
   - `Ctrl+R` / `Cmd+R` resets all state
   - Re-enter parameters manually

5. **Check browser console for errors:**
   - Press F12 to open DevTools
   - Go to Console tab
   - Red error messages indicate JavaScript problems
   - Report these to IT support

**When to Escalate:** If controls consistently fail to respond after page refresh, there may be a browser compatibility issue. Try a different browser.

---

### FAQ 5.2: Mode toggle (Explain/Internals) not switching views

**Problem:** Clicking the mode toggle between "Explain" and "Internals" does not change what is displayed.

**Possible Causes:**
1. View is updating but content is similar
2. Toggle state is not being persisted
3. Component rendering error

**Solution:**

1. **Look for subtle changes:**
   - "Explain" mode shows simplified, pedagogical explanations
   - "Internals" mode shows computation details (regime, coefficients, matrices)
   - Some sections may look similar; scroll to see differences

2. **Click the toggle multiple times:**
   - Toggle between modes 2-3 times
   - Ensure you are clicking the toggle, not nearby text

3. **Check current mode indicator:**
   - One mode should appear selected/highlighted
   - If both look the same, UI may not be updating

4. **Refresh and try again:**
   - `Ctrl+R` / `Cmd+R`
   - Navigate to the module
   - Click toggle immediately after page loads

**When to Escalate:** If mode switching is confirmed broken, students can complete the lab using whichever mode is visible. Note this in grading.

---

### FAQ 5.3: Display issues on mobile devices or tablets

**Problem:** IQCP is difficult to use on phones or tablets. Elements are cut off, too small, or overlap.

**Possible Causes:**
1. IQCP is designed for desktop browsers
2. Mobile viewport is too narrow
3. Touch events not properly handled
4. Plotly charts not resizing

**Solution:**

1. **Use a desktop or laptop computer:**
   - IQCP is optimized for screens 1024px or wider
   - Mobile is not officially supported for Lab Pack activities

2. **If mobile is necessary, try landscape mode:**
   - Rotate device to landscape orientation
   - This provides more horizontal space

3. **Use "Request Desktop Site":**
   - Chrome (Android): Menu > Desktop site
   - Safari (iOS): Press-hold refresh > Request Desktop Website

4. **Zoom out:**
   - Pinch to zoom out for a wider view
   - May make text smaller but shows more content

5. **Disable browser UI elements:**
   - Use full-screen mode if available
   - Hides address bar and tabs, giving more space

**Recommended Device Setup:**

| Device Type | Support Level | Notes |
|-------------|---------------|-------|
| Desktop/Laptop | Full | Recommended |
| Tablet (landscape) | Partial | Usable with adjustments |
| Tablet (portrait) | Limited | May have layout issues |
| Phone | Not Supported | Use desktop instead |

**When to Escalate:** Mobile support is not a priority. Students should use a computer for this lab.

---

## Additional Troubleshooting Tips

### General Debugging Steps

When any issue occurs, try these steps in order:

1. **Refresh the page** (`Ctrl+R` / `Cmd+R`)
2. **Hard refresh** (`Ctrl+Shift+R` / `Cmd+Shift+R`)
3. **Try incognito/private mode**
4. **Try a different browser**
5. **Clear browser cache and cookies for iqcp.dev**
6. **Check browser console (F12) for error messages**
7. **Verify internet connectivity**
8. **Restart the browser**

### Collecting Bug Reports

If students encounter issues that are not resolved by this FAQ, collect the following information for bug reports:

1. **Browser and version** (e.g., Chrome 120.0.6099.109)
2. **Operating system** (e.g., Windows 11, macOS Sonoma, Ubuntu 22.04)
3. **Exact URL** when the issue occurred
4. **Exact steps** to reproduce the issue
5. **Screenshot** of the problem
6. **Browser console errors** (F12 > Console > copy red text)

### Known Limitations

| Limitation | Workaround |
|------------|------------|
| No offline mode | Requires internet connection |
| No mobile support | Use desktop browser |
| Deep links may be long | Copy full URL carefully |
| SCF limited to preset systems | Cannot add custom molecules |
| Results not auto-saved | Export artifacts before closing |

---

## Contact and Support

If issues cannot be resolved using this FAQ:

1. **During class:** Instructor can assist with common issues
2. **Technical issues:** Report to course IT support
3. **Bug reports:** Submit via course management system with details above
4. **Feature requests:** Note for future development

---

## Document Revision History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2026-01-18 | Initial release |

---

**End of Troubleshooting FAQ**
