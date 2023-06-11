# BSplinesTrajOpt

 ### BSplines-based robot trajectory optimization solver.

<p align="center">
  <img src="data/torus.png" width="250" title="BSplines curve over torus manifold">
</p>

## Description

Direct-collocation techniques with B-Splines and geometric algorithms are used to generate minimum-effort control movements. When humanoid robots are considered, the computation of such movements require to solve the robot inverse dynamics, the center of mass projection and centroidal momentum, for equilibrium purposes, and the partial derivatives of all these features. This is why Lie groups theory plays an important role here to formulate recursive, exact and fast algorithms in order to retrieve first and second order information.

## Getting Started

### Dependencies

Project is created with:
* [CMake](https://cmake.org/) version: 3.10.2
* [gcc](https://gcc.gnu.org/) version: 7.5.0
* [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page) version: 3.3.4
* [CoppeliaSim](https://coppeliarobotics.com/) version: 4.0.0

### Installing

* Install [Ipopt](https://coin-or.github.io/Ipopt/INSTALL.html) at /usr/local/bin and /usr/local/include directories.

* To install this project,
```
$ ./install.sh
```
### Executing unit test

* Test analytical derivatives time for all dynamics objects,
```
$ ./unitest/test_analytic_derivatives_time.sh
```

* Test numerical derivatives time for all dynamics objects,
```
$ ./unitest/test_numeric_derivatives_time.sh
```

* Test derivatives accuracy for center of mass, centroidal momentum and inverse dynamics,
```
$ ./unitest/test_com_centroidal_dynamics_accuracy.sh
```

* Test derivatives accuracy for cost function and constraints,
```
$ ./unitest/test_cost_constraint_accuracy.sh
```

* Test time improvement when exploiting symmetry of second-order derivatives,
```
$ ./unitest/test_symmetry_exploited.sh
```

* Test trajectory optimization experiments when center-of-mass constraint is enabled,
```
$ ./unitest/test_com_experiments.sh
```

* Test trajectory optimization experiments when centroidal momentum constraint is enabled,
```
$ ./unitest/test_centroidal_experiments.sh
```
## Help

Any advise for common problems or issues please contact the authors.

## Authors

[Alvaro Paz](https://www.linkedin.com/in/alvaro-paz-anaya/)

[Gustavo Arechavaleta](https://sites.google.com/site/gustavoarechavaleta/)

## Version History
* 0.1
    * Initial Release

## License

This project is licensed under the [MIT](https://choosealicense.com/licenses/mit/) License
