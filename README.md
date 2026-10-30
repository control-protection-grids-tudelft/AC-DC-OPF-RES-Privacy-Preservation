# AC/DC OPF with RES-Side Privacy Preservation: Distributed vs. Projection Method

Code and data accompanying the MEDPOWER 2026 contribution (Paphos, Cyprus, 3–6 Nov. 2026):

> H. Li, S. Singh, A. Lekić, "AC/DC OPF with Renewable Energy Source-Side Privacy Preservation: Distributed Method vs. Projection Method," MEDPOWER 2026.

Intelligent Electrical Power Grids (IEPG) group, Delft University of Technology. Slides: [`MEDPOWER2026-HLi-No14.pdf`](MEDPOWER2026-HLi-No14.pdf).

## Overview

A TSO dispatching a hybrid AC/MTDC grid with offshore RES plants needs the RES internal operational constraints to reach the global optimum, but RES owners may not want to disclose them. This repository compares two privacy-preserving approaches that exchange only boundary (PCC) information:

- **Distributed method (ADMM):** RES-side and TSO-side subproblems are solved alternately, with dual updates, until primal and dual residuals fall below thresholds. Requires boundary-information exchange at every iteration.
- **Projection method (ECM):** each RES projects its feasible region onto the boundary P–Q plane using the epsilon-constrained method (sweep P at the PCC, find min/max Q at each step). The resulting vertex set is converted to linear constraints and shared once with the TSO.

Both are benchmarked against centralized optimization (full information) and against a centralized model that ignores RES internal constraints.

### Model

- AC grid: SOCP-relaxed AC power flow
- MTDC grid: SOCP-relaxed DC power flow
- VSC: SOCP-relaxed converter model with losses and current limits
- RES feeders: LinDistFlow (radial), with voltage limits and per-unit P/S limits
- Objective: minimize quadratic generation cost

### Test system

5-terminal MTDC network. Terminals 1–3 connect to buses 4, 6 and 8 of the IEEE 9-bus AC grid; terminals 4 and 5 connect to two RES plants, each a 30-bus radial feeder.

## Repository structure

```
codes/
  Centralized.m                     Centralized benchmark (full RES information)
  Centralized_without_RES_detail.m  Centralized OPF ignoring RES internal constraints
  Distributed.m                     ADMM-based distributed OPF
  Projection.m                      ECM-based P–Q projection + TSO-side OPF
  plot_PQ_res.m                     Computes and plots the RES P–Q feasible region
  lindistflow_up.m / _down.m        Max / min boundary Q for a given PCC active power
  case5_dc.m                        5-terminal MTDC network data
  case2_vsc.m                       VSC model data
  case32_res.m                      RES radial feeder data
figures/
  Figure5Plot.py                    Plots the P–Q regions (Fig. 5)
  P_*.mat, Q_*.mat                  P–Q region data for sweep resolutions 0.1, 0.2, 0.5
```

## Requirements

- MATLAB
- [MATPOWER](https://matpower.org) (for `case9`, `makeYbus`, `define_constants`)
- [YALMIP](https://yalmip.github.io)
- Gurobi (solver used in the scripts; can be changed in `sdpsettings`)
- Python 3 with `numpy`, `scipy`, `matplotlib` (figure only)

## Usage

Add MATPOWER and YALMIP to the MATLAB path, then run from `codes/`:

```matlab
Centralized                     % benchmark
Centralized_without_RES_detail  % impact of neglecting RES constraints
Distributed                     % ADMM
Projection                      % ECM
plot_PQ_res                     % RES P–Q region only
```

To reproduce the P–Q region figure:

```bash
cd figures
python Figure5Plot.py
```

## Key results (test case)

| | Centralized | ADMM | ECM |
|---|---|---|---|
| Solver time [s] | 0.077 | 2.117 | 1.860 |

- Neglecting RES internal constraints changes the generation cost by about 10% (3683 vs. 4100 $/h).
- ADMM converges in 26 iterations.
- ADMM and ECM both give generator outputs and costs very close to the centralized benchmark (≈4100 vs. 4099 $/h).
- ECM is faster here, but this is case-dependent: ADMM time depends on iteration count, ECM time and accuracy on the P-sweep resolution.

## Citation

If you use this code, please cite the MEDPOWER 2026 paper above.

## Contact

Haixiao Li, IEPG, TU Delft
