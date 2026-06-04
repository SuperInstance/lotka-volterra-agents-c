/* lotka-volterra-agents-c: Generalized LV dynamics for multi-species agent ecology */
#ifndef LV_AGENTS_H
#define LV_AGENTS_H

#include <stddef.h>

#define LV_MAX_SPECIES 32
#define LV_MAX_STEPS 100000

/* N-species competitive Lotka-Volterra system */
typedef struct {
    double populations[LV_MAX_SPECIES];
    double growth_rates[LV_MAX_SPECIES];
    double carrying_capacities[LV_MAX_SPECIES];
    double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES]; /* interaction matrix */
    size_t n_species;
} lv_system_t;

void lv_init(lv_system_t *sys, size_t n);
void lv_set_params(lv_system_t *sys, const double *r, const double *k,
                   const double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES]);
void lv_set_population(lv_system_t *sys, const double *pop);

/* Integration */
void lv_euler_step(lv_system_t *sys, double dt);
void lv_rk4_step(lv_system_t *sys, double dt);
void lv_simulate(lv_system_t *sys, size_t steps, double dt);

/* Analysis */
int lv_all_survive(const lv_system_t *sys, double threshold);
double lv_total_population(const lv_system_t *sys);
double lv_shannon_diversity(const lv_system_t *sys);
double lv_resilience(const lv_system_t *sys, double threshold);

/* Equilibrium */
int lv_compute_equilibrium_2species(double r1, double r2, double k1, double k2,
                                     double a12, double a21,
                                     double *n1_eq, double *n2_eq);

/* Perturbation */
typedef struct {
    double pre_populations[LV_MAX_SPECIES];
    double post_populations[LV_MAX_SPECIES];
    double recovery_steps;
    int recovered;
} perturbation_result_t;

perturbation_result_t lv_perturbation_test(lv_system_t *sys, double magnitude,
                                            size_t max_recovery, double dt, double tolerance);

#endif /* LV_AGENTS_H */
