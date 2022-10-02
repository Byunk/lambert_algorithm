# lambert_algorithm
Lambert Problem solver by Gauss', Battin's, Analytic Gradient methods

## User guide
```
lambert_gauss(T, lambda)
lambert_battin(T, lambda)
```
${T = \sqrt {\frac {8\mu}{s^3}}(t_2 - t_1)}$
and $\lambda$ represent a semi perimeter

output $x = {\sin}^2{\frac{1}{2}E}$
where $E = \frac{1}{2}(E_2-E_1)$


```
lambert_analytic_gradient(r1, r2, phi, tf, muC)
```
$r_1, r_2$ are scalar, and $\phi$ is true anomaly difference (deg)

$t_F$ is transfer time, and $\mu$ represents standard gravitational parameter

output x is a flight-path angle


## References

### Gauss method
[theory of the motion of the heavenly bodies moving about the sun in conic sections](https://www.biodiversitylibrary.org/item/58729#page/9/mode/1up)

### Introducing free parameters
An Elegant Lambert Algorithm (https://doi.org/10.2514/3.19910)

Richard H. Battin and Robin M. Vaughan

[An Improvement of Gauss' method for Solving Lambert's Problem](https://dspace.mit.edu/handle/1721.1/15587)

Vaughan, Robin M

### Analytic Gradient method
Lambert Algorithm Using Analytic Gradients (https://doi.org/10.2514/1.62091)

Jaemyung Ahn and Sang-Il Lee
