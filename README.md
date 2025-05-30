# Drum Synthesis

## Structure

- `src/`
  - `membrane_fdtd.m`: 2D drum membrane simulation using the FDTD method.
  - `membrane_modal.m`: 2D drum membrane simulation using modal synthesis.
  - `wave_3d_fdtd.m`: 3D wave propagation simulation using the FDTD method.
- `output/`
  - Contains generated plots and simulation results (e.g., `fdtd_wave_frame_t0.003s.png`).
- `computations.pdf`: Theory.

## Requirements

- MATLAB (tested with R2020b or later)

## Description of Methods
### FDTD (Finite Difference Time Domain)
Regarding the membrane, the current scheme is a 2D FDTD method that simulates the vibration of a drum membrane with a damping term. Next step is to include the excitation term coming from the cavity under the membrane. Concerning the cavity itself, current state of the simulation is the 3D wave equation solved using FDTD method and Neumann boundary conditions on every side of the cube. Next step is to include the excitation term coming from the membrane above and update the boundary conditions on this side.

### Modal Synthesis
For the membrane, the scheme used in the FDTD method is fully calculated using modal analysis. Next step is to do the same for the cavity.


## References

- [S. Bilbao, *Timpani Drum Synthesis in 3d on gpgpus*](https://dafx12.york.ac.uk/papers/dafx12_submission_36.pdf)
- [S. Bilbao, *Large Scale Physical Modeling Sound Synthesis*](https://dafx.de/paper-archive/2019/DAFx2019_paper_22.pdf)
- [S. Bilbao, *Numerical Sound Synthesis*](https://onlinelibrary.wiley.com/doi/book/10.1002/9780470749012)

---

For more details, see `computations.pdf`.
