// <<AUTOBENCH>> efficient_e1x efficient_e1 efficient_e0
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>

int8_t _Alignas(4) input_arr[490] = {
    2,   85, 81, 81,  79, 81, 82,  81,  83, 84, 0,   82, 81, 82,  79,  81, 81,
    81,  83, 84, 0,   85, 83, 84,  79,  81, 81, 80,  82, 81, 0,   85,  83, 82,
    82,  81, 80, 81,  82, 79, 0,   83,  83, 85, 81,  83, 81, 80,  80,  82, 1,
    83,  80, 83, 81,  81, 82, 81,  82,  83, 1,  81,  81, 82, 80,  82,  81, 82,
    83,  83, 23, 72,  83, 81, 74,  80,  83, 82, 85,  81, 34, 65,  83,  81, 78,
    82,  81, 83, 86,  83, 43, 67,  87,  83, 76, 79,  78, 80, 83,  84,  39, 75,
    89,  79, 76, 77,  78, 81, 83,  84,  36, 75, 92,  76, 75, 77,  78,  81, 84,
    85,  33, 73, 93,  75, 76, 78,  76,  80, 84, 86,  36, 76, 91,  76,  77, 78,
    77,  80, 84, 85,  83, 71, 84,  80,  81, 83, 81,  84, 87, 85,  113, 87, 81,
    79,  79, 80, 80,  82, 83, 80,  112, 85, 76, 76,  80, 81, 79,  84,  84, 79,
    111, 85, 76, 76,  81, 82, 80,  84,  84, 80, 110, 85, 74, 77,  81,  83, 79,
    82,  84, 80, 108, 87, 74, 77,  80,  84, 79, 82,  82, 81, 106, 88,  74, 79,
    79,  84, 79, 82,  81, 82, 104, 87,  76, 80, 79,  85, 79, 82,  81,  82, 95,
    89,  77, 80, 79,  84, 80, 81,  82,  82, 79, 90,  81, 82, 79,  84,  79, 82,
    82,  80, 67, 87,  81, 79, 77,  83,  80, 83, 81,  80, 57, 89,  80,  77, 75,
    82,  80, 83, 82,  82, 50, 91,  78,  78, 77, 81,  79, 82, 82,  82,  39, 89,
    76,  77, 78, 83,  79, 82, 83,  80,  47, 81, 78,  77, 76, 84,  77,  82, 84,
    81,  56, 72, 80,  82, 81, 86,  79,  83, 83, 81,  63, 71, 81,  82,  80, 85,
    81,  83, 81, 81,  67, 70, 77,  81,  78, 83, 81,  83, 84, 82,  68,  71, 76,
    84,  79, 84, 82,  82, 81, 80,  65,  75, 79, 84,  80, 83, 80,  81,  82, 82,
    50,  76, 75, 82,  79, 82, 79,  80,  82, 81, 45,  78, 76, 81,  78,  84, 79,
    83,  83, 80, 43,  79, 74, 78,  78,  86, 80, 84,  81, 79, 40,  80,  74, 80,
    77,  85, 80, 83,  82, 80, 35,  82,  75, 79, 80,  85, 80, 81,  82,  81, 35,
    83,  73, 78, 79,  87, 80, 83,  82,  80, 34, 84,  72, 78, 78,  85,  79, 82,
    82,  79, 33, 84,  70, 77, 77,  84,  80, 84, 83,  79, 33, 87,  73,  77, 75,
    79,  79, 84, 81,  79, 36, 89,  75,  81, 75, 78,  79, 82, 82,  81,  34, 91,
    75,  79, 74, 78,  79, 82, 84,  82,  29, 88, 76,  78, 75, 81,  81,  81, 82,
    81,  25, 88, 77,  80, 76, 80,  80,  79, 80, 83,  25, 92, 78,  79,  75, 78,
    79,  79, 82, 82,  23, 94, 81,  77,  75, 76, 80,  81, 81, 82};
int8_t _Alignas(4) output_arr[12];

int8_t (*input)[49][10][1] = (int8_t (*)[49][10][1])input_arr;
int8_t (*output)[12] = &output_arr;

const char *modelType = "kws";

#ifdef EFF_BLD_HAND_OPTIMIZED
#define model_func run_model_opt
#else
#define model_func run_model
#endif

int model_func(int8_t input[1][49][10][1], int8_t output[1][12]);

int8_t expected_output[12] = {-128, -128, -128, -128, -128, -128,
                              -128, 127,  -128, -128, -128, -128};

#define DEBUG_PRINT 1

#ifdef DEBUG_PRINT
#define TENSOR_PRINT(x, ...)      \
    do {                          \
        printf(x, ##__VA_ARGS__); \
    } while (0)
#else
#define TENSOR_PRINT(x, ...)
#endif

bool validate_output(int8_t *tensor) {
    TENSOR_PRINT("Output tensor: [%d", tensor[0]);
    bool correct = tensor[0] == expected_output[0];
    for (int i = 1; i < 12; i++) {
        TENSOR_PRINT(", %d", tensor[i]);
        if (tensor[i] != expected_output[i]) {
            correct = false;
        }
    }

    TENSOR_PRINT("]\n");

    return correct;
}

int benchmarkModel() {
    __effcc_enter_profile_region("model");
    for (int iter = 0; iter < __effcc_get_profile_iterations(); iter++)
        model_func(input, output);
    __effcc_exit_profile_region();

    bool correct = validate_output((int8_t *)output);
    if (!correct) {
        printf("[mlperftiny_%s] FAIL\n", modelType);
        return 1;
    }

    printf("[mlperftiny_%s] PASS\n", modelType);
    return 0;
}
