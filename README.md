# homuFoam
My OpenFOAM developments

## scripts


## homuPorousInterFoam
A better `porousInterFoam` with phase saturation.
That is, the total volume of liquid phase by 
```
porosity * VOF
```
is conserved provided that the solution converges sufficiently.

Currently, it has drag force model based on Ergun's equation.
The source terms by `fvOptions` can also be used.
It will be interesting to add a capillary effect in the porous medium.

