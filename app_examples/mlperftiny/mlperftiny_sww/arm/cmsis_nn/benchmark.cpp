// <<AUTOBENCH>> ambiq_apollo4p ambiq_apollo510 nxp_lpc55s69 renesas_ra8m1
#include <stdio.h>
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

#include "dsp/transform_functions.h"
#include "model_data.h"
#include "fft4k.h"
#include "featurize.h"

#include "marvin.h"

int8_t expected_output_marvin[3] = {127, -128, -128};

int8_t _Alignas(16) feat[1200];

// Smallest multiple-of-1024-sized buffer that will hold the arena based on
// interpreter.arena_used_bytes() call
constexpr int kTensorArenaSize = 20 * 1024;
alignas(16) uint8_t tensor_arena[kTensorArenaSize];

bool validate_output(TfLiteTensor* tensor, int8_t* expected) {
    bool correct = true;
    for (int i = 0; i < 3; i++) {
        if (tensor->data.int8[i] != expected[i]) {
            correct = false;
        }
    }
    return correct;
}

extern "C" int benchmarkModel() {
    printf("Starting ML Perf Tiny SWW benchmark\n");
    const tflite::Model* model = tflite::GetModel(str_ww_ref_model_int8_tflite);
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

    TfLiteStatus status;

    __effcc_enter_profile_region("model");
    for (int iter = 0; iter < __effcc_get_profile_iterations(); iter++) {
        featurize(input_data_marvin, feat);
        for (int i = 0; i < 1200; i++) {
            input->data.int8[i] = feat[i];
        }
        status = interpreter.Invoke();
    }
    __effcc_exit_profile_region();

    if (status != kTfLiteOk) {
        printf("Failed to run model: %d\n", status);
        return 1;
    }

    TfLiteTensor* output = interpreter.output(0);
    if (validate_output(output, expected_output_marvin)) {
        printf("[mlperftiny_sww] PASS\n");
    } else {
        printf("[mlperftiny_sww] FAIL\n");
        return 1;
    }

    return 0;
}

#include <stdarg.h>

void MicroPrintf(const char* format, ...) {
    printf("Got something? %s\n", format);
}
