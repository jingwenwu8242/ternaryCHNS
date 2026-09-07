# ternaryCHNS: Section 4.3, Fig. 9

This directory contains the C programs and MATLAB scripts for the **Droplet on a Concave Substrate** example.

## Cases

| Panel | Contact angle | C program | CMake target | Data and MATLAB directory |
|---|---:|---|---|---|
| (a), `fig2a` | 30 degrees | `chnsangle5b.c` | `ex2_theta30` | `data5b/` |
| (b), `fig2b` | 60 degrees | `chnsangle6b.c` | `ex2_theta60` | `data6b/` |
| (c), `fig2c` | 90 degrees | `chnsangle7b.c` | `ex2_theta90` | `data7b/` |
| (d), `fig2d` | 120 degrees | `chnsangle8b.c` | `ex2_theta120` | `data8b/` |

The four programs have the same parameters, initial condition, and boundary conditions. Only the contact angle and output directory differ.

## Main parameters

- Grid: `nr = nz = 256`
- Time step: `dt = 1e-6`
- Number of time steps: 14,000
- `Re = 1`, `We = 1`, `Pe = 0.01`
- `rho_1 = rho_2 = 1`, `eta_1 = eta_2 = 1`
- Gravity and initial velocity: zero
- Initial droplet radius: `0.2`
- Initial droplet center: `(r,z) = (0,0.25)`

## Build and run

```text
cmake -S . -B build
cmake --build build --config Release
```

On Windows with a Visual Studio generator:

```text
.\build\Release\ex2_theta30.exe
.\build\Release\ex2_theta60.exe
.\build\Release\ex2_theta90.exe
.\build\Release\ex2_theta120.exe
```

## MATLAB figures

After a simulation finishes, open its `data5b/`, `data6b/`, `data7b/`, or `data8b/` directory in MATLAB and run `show_giure3.m`.

The scripts use snapshots 3, 21, and 101 and generate `fig2a.pdf`, `fig2b.pdf`, `fig2c.pdf`, or `fig2d.pdf`. Generated simulation data are not included in this package.
