# lotka-volterra-agents-c

C implementation of generalized Lotka-Volterra dynamics for multi-species agent ecology.

Cross-language companion to [lotka-volterra-agents](https://github.com/SuperInstance/lotka-volterra-agents) (Rust).

## API

- `lv_system_t` — N-species competitive LV system
- `lv_euler_step()` / `lv_rk4_step()` — Euler and RK4 integration
- `lv_simulate()` — run N steps
- `lv_compute_equilibrium_2species()` — analytical 2-species equilibrium
- `lv_perturbation_test()` — perturb and measure recovery time
- `lv_shannon_diversity()` / `lv_resilience()` — ecological metrics

## Build & Test

```bash
gcc -o test_lv tests/test_lv.c src/lv_agents.c -lm -Wall -O2
./test_lv
```

15 tests passing. No external dependencies.

## License

MIT
