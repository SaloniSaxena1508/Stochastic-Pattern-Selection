# Stochastic Pattern Selection

This code was written for my PhD thesis in statistical physics. The problem I solved with the code is known as wavenumber selection in pattern formation.

## Non-technical introduction to pattern formation and pattern selection
Pattern formation occurs in nature when systems arrange themselves in regular patterns due to changes in their environment. PICTURE.
How often a pattern repeats itself is measured in terms of something known as a "wavenumber". Typically, a system can adopt a large range of wavenumbers, however, in practice it is observed that most systems are overwhelmingly biased toward one particular wavenumber. This is known as "pattern selection". What makes this selected wavenumber so special and how is it chosen? Here we answer that question. Our hypothesis is that random fluctuations in the surroundings drive the system toward a pattern with a specific "selected" wavenumber.

The aim of this code is to take a differential equation which is a good proxy for a real system. Random fluctuations are introduced by using Python's built-in random number generator. The differential equation is solved and the behavior of the wavenumbers is studied.

## Definitions
1. Model equation: $\frac{\partial u}{\partial t} = -\alpha u - \frac{\partial^2 u}{\partial t^2} - \frac{\partial^4 u}{\partial t^2} + \left( \frac{\partial u}{\partial x} \right)^2$

