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
