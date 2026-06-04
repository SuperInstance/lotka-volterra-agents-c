/* lotka-volterra-agents-c: Implementation */
#include "lv_agents.h"
#include <math.h>
#include <string.h>

void lv_init(lv_system_t *sys, size_t n) {
    memset(sys, 0, sizeof(*sys));
    sys->n_species = n;
    for (size_t i = 0; i < n; i++) {
        sys->alpha[i][i] = 1.0; /* self-competition */
    }
}

void lv_set_params(lv_system_t *sys, const double *r, const double *k,
                   const double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES]) {
    for (size_t i = 0; i < sys->n_species; i++) {
        sys->growth_rates[i] = r[i];
        sys->carrying_capacities[i] = k[i];
        for (size_t j = 0; j < sys->n_species; j++) {
            sys->alpha[i][j] = alpha[i][j];
        }
    }
}

void lv_set_population(lv_system_t *sys, const double *pop) {
    for (size_t i = 0; i < sys->n_species; i++) {
        sys->populations[i] = pop[i];
    }
}

/* Compute dN/dt for species i */
static double dn_dt(const lv_system_t *sys, size_t i) {
    double competition = 0.0;
    for (size_t j = 0; j < sys->n_species; j++) {
        competition += sys->alpha[i][j] * sys->populations[j] / sys->carrying_capacities[i];
    }
    return sys->growth_rates[i] * sys->populations[i] * (1.0 - competition);
}

void lv_euler_step(lv_system_t *sys, double dt) {
    double delta[LV_MAX_SPECIES];
    for (size_t i = 0; i < sys->n_species; i++) {
        delta[i] = dn_dt(sys, i) * dt;
    }
    for (size_t i = 0; i < sys->n_species; i++) {
        sys->populations[i] += delta[i];
        if (sys->populations[i] < 0.01) sys->populations[i] = 0.01;
    }
}

void lv_rk4_step(lv_system_t *sys, double dt) {
    size_t n = sys->n_species;
    double k1[LV_MAX_SPECIES], k2[LV_MAX_SPECIES];
    double k3[LV_MAX_SPECIES], k4[LV_MAX_SPECIES];
    double orig[LV_MAX_SPECIES], tmp[LV_MAX_SPECIES];

    memcpy(orig, sys->populations, n * sizeof(double));

    /* k1 */
    for (size_t i = 0; i < n; i++) k1[i] = dn_dt(sys, i);

    /* k2 */
    for (size_t i = 0; i < n; i++) tmp[i] = orig[i] + 0.5 * dt * k1[i];
    lv_set_population(sys, tmp);
    for (size_t i = 0; i < n; i++) k2[i] = dn_dt(sys, i);

    /* k3 */
    for (size_t i = 0; i < n; i++) tmp[i] = orig[i] + 0.5 * dt * k2[i];
    lv_set_population(sys, tmp);
    for (size_t i = 0; i < n; i++) k3[i] = dn_dt(sys, i);

    /* k4 */
    for (size_t i = 0; i < n; i++) tmp[i] = orig[i] + dt * k3[i];
    lv_set_population(sys, tmp);
    for (size_t i = 0; i < n; i++) k4[i] = dn_dt(sys, i);

    /* Combine */
    for (size_t i = 0; i < n; i++) {
        sys->populations[i] = orig[i] + (dt / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        if (sys->populations[i] < 0.01) sys->populations[i] = 0.01;
    }
}

void lv_simulate(lv_system_t *sys, size_t steps, double dt) {
    for (size_t i = 0; i < steps; i++) {
        lv_rk4_step(sys, dt);
    }
}

int lv_all_survive(const lv_system_t *sys, double threshold) {
    for (size_t i = 0; i < sys->n_species; i++) {
        if (sys->populations[i] < threshold) return 0;
    }
    return 1;
}

double lv_total_population(const lv_system_t *sys) {
    double total = 0.0;
    for (size_t i = 0; i < sys->n_species; i++) {
        total += sys->populations[i];
    }
    return total;
}

double lv_shannon_diversity(const lv_system_t *sys) {
    double total = lv_total_population(sys);
    if (total <= 0.0) return 0.0;
    double h = 0.0;
    for (size_t i = 0; i < sys->n_species; i++) {
        if (sys->populations[i] > 0.0) {
            double p = sys->populations[i] / total;
            h -= p * log2(p);
        }
    }
    return h;
}

double lv_resilience(const lv_system_t *sys, double threshold) {
    size_t surviving = 0;
    for (size_t i = 0; i < sys->n_species; i++) {
        if (sys->populations[i] >= threshold) surviving++;
    }
    return (double)surviving / (double)sys->n_species;
}

int lv_compute_equilibrium_2species(double r1, double r2, double k1, double k2,
                                     double a12, double a21,
                                     double *n1_eq, double *n2_eq) {
    /* Coexistence equilibrium for 2-species competitive LV:
       N1* = K1(1 - a12*K2/K1) / (1 - a12*a21)  ... but the standard form:
       dN1/dt = 0 => N1 + a12*N2 = K1
       dN2/dt = 0 => a21*N1 + N2 = K2
       Solving: N1* = (K1 - a12*K2) / (1 - a12*a21)
                N2* = (K2 - a21*K1) / (1 - a21*a12) */
    double denom = 1.0 - a12 * a21;
    if (fabs(denom) < 1e-10) return 0; /* singular */

    *n1_eq = (k1 - a12 * k2) / denom;
    *n2_eq = (k2 - a21 * k1) / denom;

    /* Feasibility: both must be positive */
    return (*n1_eq > 0.0 && *n2_eq > 0.0);
}

perturbation_result_t lv_perturbation_test(lv_system_t *sys, double magnitude,
                                            size_t max_recovery, double dt, double tolerance) {
    perturbation_result_t result = {0};

    /* Record equilibrium */
    memcpy(result.pre_populations, sys->populations, sys->n_species * sizeof(double));

    /* Apply perturbation */
    for (size_t i = 0; i < sys->n_species; i++) {
        sys->populations[i] *= (1.0 + magnitude);
        if (sys->populations[i] < 0.01) sys->populations[i] = 0.01;
    }

    /* Record post-perturbation */
    memcpy(result.post_populations, sys->populations, sys->n_species * sizeof(double));

    /* Simulate recovery */
    for (size_t step = 0; step < max_recovery; step++) {
        lv_rk4_step(sys, dt);

        /* Check if recovered (close to pre-perturbation) */
        int recovered = 1;
        for (size_t i = 0; i < sys->n_species; i++) {
            double diff = fabs(sys->populations[i] - result.pre_populations[i]);
            double rel = diff / (result.pre_populations[i] > 0.01 ? result.pre_populations[i] : 0.01);
            if (rel > tolerance) {
                recovered = 0;
                break;
            }
        }
        if (recovered) {
            result.recovery_steps = step + 1;
            result.recovered = 1;
            return result;
        }
    }

    result.recovery_steps = max_recovery;
    result.recovered = 0;
    return result;
}
