# SEOP-and-NMR


This repository contains theory an simulations of Spin Exchange Optical Pumping (SEOP) and Nuclear Magnetic Resonance (NMR) physics. 



## TODO

### Repository structure
- [x] Add requirements.txt --- 24/05/2021
- [x] Add a utils.py file --- 24/05/2021

### Theory
- [x] Write solution to open loop dual species NMRG dynamics. --- 30/05/2021
- [ ] Write solution to closed loop dual species NMRG dynamics.


### Simulations
- [x] Single specie without magnetic noise $B_{noise}=0, \omega_r>0$. --- 2/06/2021
- [x] Single specie dynamic-range simulation with and without magnetic noise $B_{noise}\leq 0, \omega_r>0$ . --- 2/06/2021
- [x] Dual species with magnetic noise $B_{noise}>0, \omega_r>0$. --- 6/06/2021
- [ ] compute PSD of $B_{noise}>0, \omega_r>0$ and find the connection to T_2.
- [ ] Dual species with magnetic noise $B_{noise}>0, \omega_r>0 and smooth magetic drift. Explore for which drift size the simulation breaks (there is no noise cancellation).




## List of Notebooks

| #   | Subject                                         | Colab             | Nbviwer               |
|:----:|------------------------------------------------|:-----------------:|:---------------------:|
| 0    | Gaussian White Noise simulation                | [![Open In Collab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/RoyElkabetz/SEOP-and-NMR/blob/main/src/How_to_generate_White_Gaussian_Noise.ipynb)        | [![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/RoyElkabetz/SEOP-and-NMR/blob/main/src/How_to_generate_White_Gaussian_Noise.ipynb)|
| 1   | Single species NMR simulation without magnetic noise                   | [![Open In Collab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/RoyElkabetz/SEOP-and-NMR/blob/main/src/single_specie_experiment.ipynb)        | [![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/RoyElkabetz/SEOP-and-NMR/blob/main/src/single_specie_experiment.ipynb)|
| 1   | Single species NMR simulation with magnetic noise                   | [![Open In Collab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/RoyElkabetz/SEOP-and-NMR/blob/main/src/single_specie_experiment.ipynb)        | [![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/RoyElkabetz/SEOP-and-NMR/blob/main/src/single_specie_experiment.ipynb)|
| 2   | Dual species NMR simulation                   | [![Open In Collab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/RoyElkabetz/SEOP-and-NMR/blob/main/src/dual_specie_experiment.ipynb)        | [![nbviewer](https://raw.githubusercontent.com/jupyter/design/master/logos/Badges/nbviewer_badge.svg)](https://nbviewer.jupyter.org/github/RoyElkabetz/SEOP-and-NMR/blob/main/src/dual_specie_experiment.ipynb)|

