# Simulations

## Make skies

```
cd sky_models/
mkdir data
spython scripts/make_pygdsm.py
spython scripts/make_const.py
```

## Make Beams
```
cd beams
mkdir data
spython scripts/make_beam.py
```

## Run simulation
```
spython scripts/run_sim.py
```

This will currently fail unless you have a hacked verion of pyuvsim. When it would also fail, but a bit later. :)
