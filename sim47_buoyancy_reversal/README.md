# Buoyancy reversal

Simulation based on  
Grabowski 1993 "Cumulus entrainment, fine-scale mixing, and buoyancy reversal"  
<https://doi.org/10.1002/qj.49711951305>

To override default parameters defined in `gen`,
assign new values in `par.py`.

To start the simulation, setup Aphros environment and run
```
make cleanrun
```

To create images and video using ParaView
```
cd vis
./run
```

Expected result is in
[vis/precomputed/sim47_buoyancy_reversal.mp4](vis/precomputed/sim47_buoyancy_reversal.mp4)

## Slides

Slides describing an earlier prototype for buoyancy reversal.

<https://pkarnakov.github.io/slides/buoyancy_reversal>

## Examples

* Visible effect of buoyancy reversal
  (smaller initial velocity, stronger buoyancy)

```
# par.py
np = 2
tmax = 2

reversal = False

buo_bc_y = [0]
buo_bc_values = [10]
buo_chis_y = [0]
buo_chis_values = [0.5]
buo_bs_y = [0]
buo_bs_values = [-10 if reversal else 0]

psi_unit_u = 0.5
```
