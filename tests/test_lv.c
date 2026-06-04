/* Tests for lotka-volterra-agents-c */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "../src/lv_agents.h"

#define ASSERT_FEQ(a, b, tol) do { \
    if (fabs((a) - (b)) > (tol)) { \
        fprintf(stderr, "FAIL %s:%d: %f != %f\n", __FILE__, __LINE__, (double)(a), (double)(b)); \
        return 1; \
    } \
} while(0)

int test_init() {
    lv_system_t sys;
    lv_init(&sys, 5);
    assert(sys.n_species == 5);
    for (int i = 0; i < 5; i++) assert(sys.alpha[i][i] == 1.0);
    return 0;
}

int test_euler_step() {
    lv_system_t sys;
    lv_init(&sys, 2);
    double r[] = {1.0, 1.0};
    double k[] = {100.0, 100.0};
    double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES] = {{0}};
    alpha[0][0] = 1.0; alpha[0][1] = 0.5;
    alpha[1][0] = 0.5; alpha[1][1] = 1.0;
    lv_set_params(&sys, r, k, alpha);

    double pop[] = {50.0, 50.0};
    lv_set_population(&sys, pop);

    lv_euler_step(&sys, 0.01);
    /* Should change but not explode */
    assert(sys.populations[0] > 0.0);
    assert(sys.populations[1] > 0.0);
    assert(sys.populations[0] < 1000.0);
    return 0;
}

int test_2species_equilibrium() {
    double n1, n2;
    int feasible = lv_compute_equilibrium_2species(1.0, 1.0, 100.0, 100.0, 0.5, 0.5, &n1, &n2);
    assert(feasible);
    /* N1* = (100 - 0.5*100) / (1 - 0.25) = 50/0.75 = 66.67 */
    ASSERT_FEQ(n1, 66.667, 1.0);
    ASSERT_FEQ(n2, 66.667, 1.0);
    return 0;
}

int test_2species_convergence() {
    lv_system_t sys;
    lv_init(&sys, 2);
    double r[] = {1.0, 1.0};
    double k[] = {100.0, 100.0};
    double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES] = {{0}};
    alpha[0][0] = 1.0; alpha[0][1] = 0.5;
    alpha[1][0] = 0.5; alpha[1][1] = 1.0;
    lv_set_params(&sys, r, k, alpha);

    double pop[] = {20.0, 20.0};
    lv_set_population(&sys, pop);
    lv_simulate(&sys, 10000, 0.01);

    ASSERT_FEQ(sys.populations[0], 66.667, 5.0);
    ASSERT_FEQ(sys.populations[1], 66.667, 5.0);
    return 0;
}

int test_2species_both_survive() {
    lv_system_t sys;
    lv_init(&sys, 2);
    double r[] = {1.0, 1.0};
    double k[] = {100.0, 100.0};
    double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES] = {{0}};
    alpha[0][0] = 1.0; alpha[0][1] = 0.5;
    alpha[1][0] = 0.5; alpha[1][1] = 1.0;
    lv_set_params(&sys, r, k, alpha);
    double pop[] = {50.0, 50.0};
    lv_set_population(&sys, pop);
    lv_simulate(&sys, 5000, 0.01);

    assert(lv_all_survive(&sys, 1.0));
    return 0;
}

int test_5species_all_survive() {
    lv_system_t sys;
    lv_init(&sys, 5);
    double r[] = {1.0, 0.8, 1.2, 0.7, 0.5};
    double k[] = {100.0, 100.0, 100.0, 100.0, 100.0};
    double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES] = {{0}};
    alpha[0][0] = 1.0; alpha[0][1] = 0.3; alpha[0][2] = 0.2; alpha[0][3] = 0.3; alpha[0][4] = 0.2;
    alpha[1][0] = 0.3; alpha[1][1] = 1.0; alpha[1][2] = 0.3; alpha[1][3] = 0.2; alpha[1][4] = 0.2;
    alpha[2][0] = 0.2; alpha[2][1] = 0.3; alpha[2][2] = 1.0; alpha[2][3] = 0.3; alpha[2][4] = 0.3;
    alpha[3][0] = 0.3; alpha[3][1] = 0.2; alpha[3][2] = 0.3; alpha[3][3] = 1.0; alpha[3][4] = 0.2;
    alpha[4][0] = 0.2; alpha[4][1] = 0.2; alpha[4][2] = 0.3; alpha[4][3] = 0.2; alpha[4][4] = 1.0;
    lv_set_params(&sys, r, k, alpha);
    double pop[] = {20.0, 20.0, 20.0, 20.0, 20.0};
    lv_set_population(&sys, pop);
    lv_simulate(&sys, 10000, 0.01);

    assert(lv_all_survive(&sys, 1.0));
    ASSERT_FEQ(lv_resilience(&sys, 1.0), 1.0, 0.01);
    return 0;
}

int test_total_population() {
    lv_system_t sys;
    lv_init(&sys, 3);
    double pop[] = {10.0, 20.0, 30.0};
    lv_set_population(&sys, pop);
    ASSERT_FEQ(lv_total_population(&sys), 60.0, 0.01);
    return 0;
}

int test_shannon_diversity_uniform() {
    lv_system_t sys;
    lv_init(&sys, 5);
    double pop[] = {20.0, 20.0, 20.0, 20.0, 20.0};
    lv_set_population(&sys, pop);
    ASSERT_FEQ(lv_shannon_diversity(&sys), log2(5.0), 0.01);
    return 0;
}

int test_shannon_diversity_dominated() {
    lv_system_t sys;
    lv_init(&sys, 3);
    double pop[] = {99.0, 0.5, 0.5};
    lv_set_population(&sys, pop);
    assert(lv_shannon_diversity(&sys) < 0.2);
    return 0;
}

int test_resilience_partial() {
    lv_system_t sys;
    lv_init(&sys, 5);
    double pop[] = {10.0, 10.0, 0.001, 0.001, 0.001};
    lv_set_population(&sys, pop);
    ASSERT_FEQ(lv_resilience(&sys, 1.0), 0.4, 0.01);
    return 0;
}

int test_resilience_full() {
    lv_system_t sys;
    lv_init(&sys, 5);
    double pop[] = {10.0, 10.0, 10.0, 10.0, 10.0};
    lv_set_population(&sys, pop);
    ASSERT_FEQ(lv_resilience(&sys, 1.0), 1.0, 0.01);
    return 0;
}

int test_equilibrium_infeasible() {
    double n1, n2;
    /* Strong competition → exclusion */
    int feasible = lv_compute_equilibrium_2species(1.0, 1.0, 100.0, 100.0, 2.0, 2.0, &n1, &n2);
    /* (100 - 200)/(1-4) = -100/-3 = 33.3, (100-200)/(1-4) = 33.3 — actually feasible! */
    /* Try asymmetric exclusion */
    feasible = lv_compute_equilibrium_2species(1.0, 1.0, 100.0, 100.0, 3.0, 0.1, &n1, &n2);
    /* N1* = (100 - 300)/(1-0.3) = -200/0.7 < 0 → infeasible */
    assert(!feasible);
    return 0;
}

int test_perturbation_recovery() {
    lv_system_t sys;
    lv_init(&sys, 2);
    double r[] = {1.0, 1.0};
    double k[] = {100.0, 100.0};
    double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES] = {{0}};
    alpha[0][0] = 1.0; alpha[0][1] = 0.5;
    alpha[1][0] = 0.5; alpha[1][1] = 1.0;
    lv_set_params(&sys, r, k, alpha);

    /* Run to equilibrium first */
    double pop[] = {50.0, 50.0};
    lv_set_population(&sys, pop);
    lv_simulate(&sys, 10000, 0.01);

    perturbation_result_t result = lv_perturbation_test(&sys, 0.5, 5000, 0.01, 0.1);
    assert(result.recovered);
    assert(result.recovery_steps > 0);
    assert(result.recovery_steps < 5000);
    return 0;
}

int test_rk4_vs_euler() {
    /* RK4 should converge faster than Euler */
    lv_system_t sys_euler, sys_rk4;
    lv_init(&sys_euler, 2);
    lv_init(&sys_rk4, 2);

    double r[] = {1.0, 1.0};
    double k[] = {100.0, 100.0};
    double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES] = {{0}};
    alpha[0][0] = 1.0; alpha[0][1] = 0.5;
    alpha[1][0] = 0.5; alpha[1][1] = 1.0;
    lv_set_params(&sys_euler, r, k, alpha);
    lv_set_params(&sys_rk4, r, k, alpha);

    double pop[] = {20.0, 20.0};
    lv_set_population(&sys_euler, pop);
    lv_set_population(&sys_rk4, pop);

    /* Same steps, RK4 should be closer to equilibrium */
    for (int i = 0; i < 1000; i++) {
        lv_euler_step(&sys_euler, 0.01);
        lv_rk4_step(&sys_rk4, 0.01);
    }

    double eq = 66.667;
    double euler_err = fabs(sys_euler.populations[0] - eq);
    double rk4_err = fabs(sys_rk4.populations[0] - eq);
    assert(rk4_err <= euler_err + 0.01); /* RK4 at least as good */
    return 0;
}

int test_single_species_logistic() {
    /* Single species → logistic growth to K */
    lv_system_t sys;
    lv_init(&sys, 1);
    double r[] = {1.0};
    double k[] = {100.0};
    double alpha[LV_MAX_SPECIES][LV_MAX_SPECIES] = {{0}};
    alpha[0][0] = 1.0;
    lv_set_params(&sys, r, k, alpha);
    double pop[] = {10.0};
    lv_set_population(&sys, pop);
    lv_simulate(&sys, 10000, 0.01);

    ASSERT_FEQ(sys.populations[0], 100.0, 1.0);
    return 0;
}

int main(void) {
    int failures = 0;
    #define RUN(name) do { \
        int r = test_##name(); \
        if (r == 0) printf("  PASS: %s\n", #name); \
        else { printf("  FAIL: %s\n", #name); failures++; } \
    } while(0)

    printf("Running tests...\n");
    RUN(init);
    RUN(euler_step);
    RUN(2species_equilibrium);
    RUN(2species_convergence);
    RUN(2species_both_survive);
    RUN(5species_all_survive);
    RUN(total_population);
    RUN(shannon_diversity_uniform);
    RUN(shannon_diversity_dominated);
    RUN(resilience_partial);
    RUN(resilience_full);
    RUN(equilibrium_infeasible);
    RUN(perturbation_recovery);
    RUN(rk4_vs_euler);
    RUN(single_species_logistic);

    int n = 15;
    printf("\n%d/%d passed\n", n - failures, n);
    return failures;
}
