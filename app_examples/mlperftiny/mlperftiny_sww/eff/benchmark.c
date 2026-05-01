// <<AUTOBENCH>> efficient_e1 efficient_e1x
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>

#include "fft4k.h"
#include "featurize.h"
#include "marvin.h"

int8_t expected_output_marvin[3] = {127, -128, -128};

int8_t _Alignas(4) output_arr_marvin[3];

const char *modelType = "sww";
#ifdef EFF_BLD_HAND_OPTIMIZED
#define model_func run_model_opt
#else
#define model_func run_model
#endif

int model_func(const int8_t input_0[1][30][1][40], int8_t output_0[1][3]);

#define DEBUG_PRINT
#ifdef DEBUG_PRINT
#define TENSOR_PRINT(x, ...)      \
    do {                          \
        printf(x, ##__VA_ARGS__); \
    } while (0)
#else
#define TENSOR_PRINT(x, ...)
#endif

bool validate_output(int8_t *tensor, int8_t *expected) {
    TENSOR_PRINT("Output tensor: [%d", tensor[0]);
    bool correct = tensor[0] == expected[0];
    for (int i = 1; i < 3; i++) {
        TENSOR_PRINT(", %d", tensor[i]);
        if (tensor[i] != expected[i]) {
            correct = false;
        }
    }

    TENSOR_PRINT("]\n");

    return correct;
}

int8_t _Alignas(4) feat[1200];

int benchmarkModel() {
    __effcc_enter_profile_region("model");
    for (int iter = 0; iter < __effcc_get_profile_iterations(); iter++) {
        featurize(input_data_marvin, feat);
        model_func((int8_t (*)[30][1][40])feat, &output_arr_marvin);
    }
    __effcc_exit_profile_region();

    // results are incorrect on hardware for now,
    // likely due to this wide integer issue:
    // https://github.com/EfficientComputer/riptools/issues/3116
    // results verified through simulation of optimized version through QEMU
    if (validate_output(output_arr_marvin, expected_output_marvin)) {
        printf("[mlperftiny_%s] PASS\n", modelType);
    } else {
        printf("[mlperftiny_%s] FAIL\n", modelType);
        return 1;
    }
    return 0;
}
