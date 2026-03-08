---
layout: two-cols
layoutClass: gap-8
---

# Problem

Calculating of the energy of 1-D harmonic oscillator

Hamiltonian:
$$
\hat{H} = - \frac{1}{2} \frac{\mathrm{d}^2}{\mathrm{d}x^2} + \frac{1}{2}x^2
$$

Trial wave function:
$$
\psi(x) = \exp \left( -\frac{1}{2}x^2 \right)
$$

Energy:
$$
E = \frac{\langle \psi_0 | \hat{H} | \psi_0 \rangle}{\langle \psi_0 | \psi_0 \rangle} = \frac{1}{2}
$$

::right::

# Results

<div class="chart-wrap">
  <p class="legend-inline">
    <span style="color:#1f77b4;">■</span> FDM
    <span style="color:#ff7f0e;">■</span> QTT+TCI
  </p>

```mermaid {scale: 0.62}
---
config:
  themeVariables:
    xyChart:
      plotColorPalette: "#1f77b4, #ff7f0e"
---
xychart-beta
  title "Runtime Comparison"
  x-axis "c" [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
  y-axis "Time (ms)" 0 --> 420
  line "FDM" [0.022541, 0.034369, 0.074828, 0.093436, 0.180841, 0.361011, 0.725843, 1.261, 3.093, 6.205, 13.965, 31.262, 59.420, 178.922, 400.168]
  line "QTT+TCI" [1.645, 2.415, 3.593, 5.042, 6.342, 11.272, 11.710, 10.414, 13.150, 13.704, 13.348, 14.062, 24.675, 16.491, 15.956]
```
</div>

For $c = 16, n = 2^c = 65536$, 

- 0.4999999272383 ← FDM
- 0.4999999272384 ← QTT

<style>
.slidev-layout.two-cols .chart-wrap {
  position: relative;
}

.slidev-layout.two-cols .mermaid {
  max-width: 100%;
  overflow: hidden;
}

.slidev-layout.two-cols .mermaid svg {
  width: 100% !important;
  height: auto !important;
}

.slidev-layout.two-cols .mermaid text {
  font-size: 10px !important;
}

/* Thicken series lines in Mermaid xychart */
.slidev-layout.two-cols .mermaid svg g.plot path,
.slidev-layout.two-cols .mermaid svg .line-plot-0,
.slidev-layout.two-cols .mermaid svg .line-plot-1 {
  stroke-width: 3px !important;
}

.slidev-layout.two-cols .legend-inline {
  position: absolute;
  top: 2.15rem;
  right: 0.7rem;
  z-index: 20;
  margin: 0;
  padding: 0.12rem 0.42rem;
  font-size: 0.68rem;
  border-radius: 9999px;
  background: rgba(255, 255, 255, 0.85);
  box-shadow: 0 1px 4px rgba(0, 0, 0, 0.12);
}

.slidev-layout.two-cols .legend-sep {
  margin: 0 0.35rem;
  opacity: 0.6;
}

.chart-wrap {
  position: relative;
  /* height: 0; */
  /* margin-top: -14rem; */
  /* overflow: visible; */
}

.legend-inline {
    position: absolute;
    top: 1.5em;
    left: 3.5em;
}
</style>
