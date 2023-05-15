# Schnakenberg-Model-Analysis
Let us explore the techniques and instruments presented in this paper by examining them within a simpler system that allows for analytical solutions and result comparisons. We will focus on the Schnakenberg system, a partial differential equation system described as follows:

In this system, we will analyze the unique steady state and derive stability conditions by solving it analytically. The Schnakenberg system is particularly convenient for this purpose. It consists of two equations governing the dynamics of the variables $u$ and $v$, represented as:

$$
  \frac{\partial u}{\partial t} = \Delta u + \eta(a - u + u^2v) = \Delta u + \eta f(u,v),
$$
  \frac{\partial v}{\partial t} = d\Delta v + \eta(b - u^2v) = d \Delta v + \eta g(u,v),

In this context, both species are uniformly generated throughout the domain. The variable $u$ undergoes linear decay, while $v$ undergoes nonlinear and autocatalytic conversion to $u$. The diffusion rate $d$ influences the relative speed of dispersal between the two species, and the parameter $\eta$ determines the balance between diffusion and the chemical reaction.


## Hopf bifurcation
![Alt text](figures/timestepHopfP0.jpg?raw=true "Hopf")


## Turing Bifurcation.
![Alt text](figures/timestepTuringP1.jpg?raw=true "Turing")
