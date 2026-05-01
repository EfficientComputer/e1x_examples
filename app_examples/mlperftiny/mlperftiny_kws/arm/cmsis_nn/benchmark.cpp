// <<AUTOBENCH>> ambiq_apollo4p ambiq_apollo510 nxp_lpc55s69 renesas_ra8m1
#define TF_LITE_STATIC_MEMORY
#define CMSIS_NN
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++11-narrowing"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-compare"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#include "tensorflow/lite/micro/kernels/micro_ops.h"
#include "tensorflow/lite/micro/micro_interpreter.h"
#include "tensorflow/lite/micro/micro_mutable_op_resolver.h"
#include "tensorflow/lite/schema/schema_generated.h"
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#pragma clang diagnostic pop
#pragma clang diagnostic pop

#include "model_data.h"

constexpr int kTensorArenaSize = 90 * 1024;
alignas(16) uint8_t tensor_arena[kTensorArenaSize];

alignas(16) int8_t input_arr[490] = {
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

bool validate_output(TfLiteTensor* tensor) {
    int max = 0;
    for (int i = 0; i < 12; i++) {
        if (tensor->data.int8[i] > tensor->data.int8[max]) {
            max = i;
        }
    }

    return max == 7;
}

extern "C" int benchmarkModel() {
    printf("Starting ML Perf Tiny KWS benchmark\n");

    const tflite::Model* model = tflite::GetModel(g_kws_model_data);
    if (model->version() != TFLITE_SCHEMA_VERSION) {
        printf("Model schema version mismatch!\n");
        return 1;
    }

    static tflite::MicroMutableOpResolver<7> resolver;

    resolver.AddAdd();
    resolver.AddFullyConnected();
    resolver.AddConv2D();
    resolver.AddDepthwiseConv2D();
    resolver.AddReshape();
    resolver.AddSoftmax();
    resolver.AddAveragePool2D();

    static tflite::MicroInterpreter interpreter(model, resolver, tensor_arena,
                                                kTensorArenaSize, nullptr);

    if (interpreter.AllocateTensors() != kTfLiteOk) {
        printf("Failed to allocate tensors\n");
        return 1;
    }

    printf("size: %d\n", interpreter.arena_used_bytes());

    TfLiteTensor* input = interpreter.input(0);
    for (int i = 0; i < 3072; i++) {
        input->data.int8[i] = input_arr[i];
    }

    TfLiteStatus status;
    __effcc_enter_profile_region("model");
    for (int iter = 0; iter < __effcc_get_profile_iterations(); iter++)
        status = interpreter.Invoke();
    __effcc_exit_profile_region();

    if (status != kTfLiteOk) {
        printf("Failed to run model: %d\n", status);
        return 1;
    }

    TfLiteTensor* output = interpreter.output(0);
    if (validate_output(output)) {
        printf("[mlperftiny_kws] PASS\n");
    } else {
        printf("[mlperftiny_kws] FAIL\n");
        return 1;
    }

    return 0;
}

#include <stdarg.h>

void MicroPrintf(const char* format, ...) {
    printf("Got something? %s\n", format);
}