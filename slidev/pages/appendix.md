# Appendix: Method

## FDM

$$
\pmb{H}
= - \frac{\hbar^2}{2m} \cdot \frac{1}{\Delta x^2}
  \left(\begin{array}{ccccccc}
    -2 & 1 & 0 & \ldots \\
    1 & -2 & 1 & \ldots \\
    0 & 1 & -2 & \ldots \\
    \vdots & \vdots & \vdots & \ddots \\
  \end{array}\right)
  +
  \left(\begin{array}{ccccccc}
    V(-50) & 0 & 0 & \ldots \\
    0 & V(-50+\Delta x) & 0 & \ldots \\
    0 & 0 & V(-50+2\Delta x) & \ldots \\
    \vdots & \vdots & \vdots & \ddots \\
  \end{array}\right),
\\
\Delta x = 0.1,
~~~\pmb{\psi}
=
\left(\begin{array}{c}
  \psi(-50) \\
  \psi(-50+\Delta x) \\
  \psi(-50+2\Delta x) \\
  \vdots \\
\end{array}\right),
~~~E = \frac{\langle\psi|\hat{H}|\psi\rangle}{\langle\psi|\psi\rangle}
\simeq \frac{\pmb{\psi}^\ast \pmb{H} \pmb{\psi}}{\pmb{\psi}^\ast \pmb{\psi}}.
$$

## QTT

https://arxiv.org/pdf/2505.17046 by Lucas Arenstein