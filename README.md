# Lotka-Volterra Agents — Multi-Species Ecological Dynamics in C

The **Lotka-Volterra (LV) equations** are a system of coupled ordinary differential equations that model how populations of interacting species change over time. Originally formulated by Alfred Lotka (1925) and Vito Volterra (1926) for predator-prey pairs, the generalized N-species form describes competition, mutualism, and parasitism across an entire ecosystem. This library implements the generalized LV system in portable C with **Euler** and **RK4** numerical integration, analytical equilibrium computation, perturbation testing, and ecological diversity metrics.

## Why It Matters

Ecological dynamics are not just about foxes and rabbits. The same equations govern **agent ecosystems**: competing AI strategies in a marketplace, microservice populations sharing infrastructure, neural architecture variants competing for accuracy. When two agent types compete for the same resource (attention, compute, users), their population dynamics follow the same equations as two barnacle species on a rocky shore.

Understanding whether an agent ecosystem will reach **coexistence** (all species survive), **competitive exclusion** (one species dominates), or **oscillatory cycles** (populations boom and bust) is critical for fleet design. This library lets you simulate these dynamics on-device — no Python runtime, no external solver, just C and `libm`. It runs on embedded devices, in WebAssembly, and anywhere a C compiler exists.

## How It Works

### The Generalized Lotka-Volterra Equation

For N species with populations N₁, N₂, ..., Nₙ:

```
dNᵢ/dt = rᵢ · Nᵢ · (1 - Σⱼ αᵢⱼ · Nⱼ / Kᵢ)
```

Where:
- **rᵢ** = intrinsic growth rate of species i (exponential growth if alone)
- **Kᵢ** = carrying capacity of species i (environmental limit)
- **αᵢⱼ** = interaction coefficient (effect of species j on species i)
- **αᵢᵢ = 1** by convention (self-competition)

The term `rᵢ · Nᵢ` is exponential growth. The term `(1 - Σⱼ αᵢⱼ · Nⱼ / Kᵢ)` is the **logistic brake**: as populations grow, competition for resources slows growth.

**Interpretation of α:**
- αᵢⱼ < 1: species j has less impact on species i than i's own conspecifics (weak competition)
- αᵢⱼ = 1: species j and i are equivalent competitors
- αᵢⱼ > 1: species j displaces species i (strong competition / predation)

### Numerical Integration: Euler vs RK4

The library provides two integrators:

**Forward Euler** — first-order, simplest:

```
N(t + dt) = N(t) + dt · f(N(t))
```

Local truncation error: **O(dt²)**. Global error: **O(dt)**. Stable for dt < 1/rₘₐₓ.

**Classical RK4** — fourth-order Runge-Kutta:

```
k₁ = f(N)
k₂ = f(N + ½dt·k₁)
k₃ = f(N + ½dt·k₂)
k₄ = f(N + dt·k₃)
N(t + dt) = N(t) + (dt/6)(k₁ + 2k₂ + 2k₃ + k₄)
```

Local truncation error: **O(dt⁵)**. Global error: **O(dt⁴)**. The library's tests confirm RK4 converges to equilibrium at least as fast as Euler for the same step count.

**Big-O per step:**
| Integrator | Evaluations of f | Time |
|---|---|---|
| Euler | 1 | O(N²) |
| RK4 | 4 | O(N²) |

Each evaluation of f requires computing the competition sum Σⱼ αᵢⱼ · Nⱼ for each species — an O(N²) matrix-vector product. For N = 32 (the library maximum), this is 1,024 multiply-adds per evaluation, trivial for modern hardware.

### Equilibrium Analysis

For the 2-species competitive case, setting dN₁/dt = dN₂/dt = 0 yields a linear system:

```
N₁* + α₁₂ · N₂* = K₁
α₂₁ · N₁* + N₂* = K₂
```

Solving by Cramer's rule:

```
N₁* = (K₁ - α₁₂·K₂) / (1 - α₁₂·α₂₁)
N₂* = (K₂ - α₂₁·K₁) / (1 - α₁₂·α₂₁)
```

**Coexistence is feasible** if and only if N₁* > 0 AND N₂* > 0, which requires:

```
α₁₂ < K₁/K₂  AND  α₂₁ < K₂/K₁
```

If either inequality is violated, one species excludes the other — this is the **competitive exclusion principle** (Gause's law).

### Perturbation Testing

The `lv_perturbation_test()` function:

1. Records current populations (assumed to be at or near equilibrium)
2. Applies a multiplicative shock: `Nᵢ → Nᵢ × (1 + magnitude)`
3. Simulates recovery using RK4
4. Measures steps until all populations return within tolerance

Recovery time quantifies **ecological resilience**: systems near a bifurcation point recover slowly (critical slowing down), a known early-warning signal for regime shifts.

### Shannon Diversity Index

```
H = -Σᵢ pᵢ · log₂(pᵢ)
```

Where pᵢ = Nᵢ / Σⱼ Nⱼ. Maximum H = log₂(N_species) when all species are equally abundant. H = 0 when one species dominates completely.

## Quick Start

```bash
# Build
gcc -o test_lv tests/test_lv.c src/lv_agents.c -lm -Wall -O2

# Run all 15 tests
./test_lv
```

```c
#include "lv_agents.h"

int main(void) {
    lv_system_t sys;
    lv_init(&sys, 2);  // 2-species system

    double r[] = {1.0, 1.0};           // equal growth rates
    double k[] = {100.0, 100.0};       // equal carrying capacities
    double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES] = {{0}};
    alpha[0][0] = 1.0; alpha[0][1] = 0.5;  // weak competition
    alpha[1][0] = 0.5; alpha[1][1] = 1.0;

    lv_set_params(&sys, r, k, alpha);
    double pop[] = {20.0, 20.0};
    lv_set_population(&sys, pop);

    // Simulate 10,000 steps with dt=0.01
    lv_simulate(&sys, 10000, 0.01);

    // Check equilibrium: both species coexist at ~66.67
    printf("N1=%.1f  N2=%.1f  H=%.3f\n",
           sys.populations[0], sys.populations[1],
           lv_shannon_diversity(&sys));

    // Perturbation test: 50% shock
    perturbation_result_t pr = lv_perturbation_test(&sys, 0.5, 5000, 0.01, 0.1);
    printf("Recovered: %s, steps: %zu\n",
           pr.recovered ? "yes" : "no", pr.recovery_steps);

    return 0;
}
```

## API

### Core Types

```c
typedef struct {
    double populations[LV_MAX_SPECIES];      // current N_i
    double growth_rates[LV_MAX_SPECIES];      // r_i
    double carrying_capacities[LV_MAX_SPECIES]; // K_i
    double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES]; // interaction matrix
    size_t n_species;
} lv_system_t;
```

### Initialization

| Function | Description |
|---|---|
| `lv_init(sys, n)` | Initialize n-species system, set αᵢᵢ = 1 |
| `lv_set_params(sys, r, k, alpha)` | Set growth rates, capacities, interaction matrix |
| `lv_set_population(sys, pop)` | Set initial populations |

### Integration

| Function | Order | Error |
|---|---|---|
| `lv_euler_step(sys, dt)` | 1st | O(dt) global |
| `lv_rk4_step(sys, dt)` | 4th | O(dt⁴) global |
| `lv_simulate(sys, steps, dt)` | N RK4 steps | — |

### Analysis

| Function | Returns | Description |
|---|---|---|
| `lv_all_survive(sys, threshold)` | `int` | 1 if all Nᵢ ≥ threshold |
| `lv_total_population(sys)` | `double` | Σᵢ Nᵢ |
| `lv_shannon_diversity(sys)` | `double` | Shannon H index (bits) |
| `lv_resilience(sys, threshold)` | `double` | Fraction of species surviving |

### Equilibrium & Perturbation

| Function | Description |
|---|---|
| `lv_compute_equilibrium_2species(...)` | Analytical N₁*, N₂* for 2-species competitive LV |
| `lv_perturbation_test(sys, magnitude, max_recovery, dt, tolerance)` | Shock the system, measure recovery |

## Architecture Notes

This C library is the **on-device ecological engine** for the SuperInstance constellation. In the conservation law **γ + η = C**, Lotka-Volterra dynamics describe how γ (generation energy) distributes across agent species. When γ is partitioned among competing agent types, the LV equations predict which types will thrive, which will go extinct, and how the system responds to shocks.

The library runs on edge devices alongside `strategy-ecology-c`, which adds game-theoretic strategy interactions on top of the LV dynamics. Together they form the **ecological prediction layer** that lets a room anticipate which agent configurations will be stable before deploying them. See the [SuperInstance Architecture](https://github.com/SuperInstance/SuperInstance/blob/main/ARCHITECTURE.md).

## References

1. Lotka, A. J. *Elements of Physical Biology* (1925) — original formulation
2. Volterra, V. "Variazioni e fluttuazioni del numero d'individui in specie animali conviventi" (1926) — independent discovery, *Mem. R. Accad. Naz. dei Lincei* 6(2)
3. Case, T. J. *An Illustrated Guide to Theoretical Ecology* (2000) — accessible modern treatment with equilibrium analysis

## License

MIT
