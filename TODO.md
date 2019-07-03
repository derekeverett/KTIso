## To Do for ITA

### $v_z$ propagation

right now trying to use a table of values of v_z with roots separated cubicly in v_z

take $v_z{,i}$ evenly spaced in [0,1], and then take $v_{z,i}^3$ as the quadrature abcissa

still not sure what to do for the quadrature weights? left / right differences?

Also, if only using $v_z$ in [0,1] rather than [-1,1], need to make sure that boundary conditions for $v_z$ propagation are set properly using the known values, by assuming symmetry.

 ### Questions

If I use nonuniform grid spacing in $v_z$, how small does $dt$ need to be to resolve the maccormack step in $v_z$. Maybe I can use a uniform grid in $v_z$ but take a small $dt$ s.t. $dt < dv_z / 4$ ?

