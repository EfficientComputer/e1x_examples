// <<AUTOBENCH>> renesas_ra8p1
/*
 * Runs a Vela-optimized TFLite model on the Ethos-U55 NPU using the
 * TFLite Micro MicroInterpreter with the EthosU custom operator kernel.
 */

/* Suppress warnings from TFLite Micro headers */
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wc++11-narrowing"
#pragma clang diagnostic ignored "-Wconversion"
#pragma clang diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic ignored "-Wsign-compare"
#pragma clang diagnostic ignored "-Wunused-parameter"
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"

#include "tensorflow/lite/micro/micro_interpreter.h"
#include "tensorflow/lite/micro/micro_mutable_op_resolver.h"
#include "tensorflow/lite/micro/kernels/ethosu.h"
#include "tensorflow/lite/schema/schema_generated.h"

#pragma clang diagnostic pop

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "str_ww_ref_model_vela.tflite.h"
#include "dsp/transform_functions.h"
#include "fft4k.h"
#include "featurize.h"

#include "marvin.h"

// Smallest multiple-of-1024-sized buffer that will hold the arena based on
// interpreter.arena_used_bytes() call
#define ARENA_SIZE_BYTES (26 * 1024)

int8_t expected_output_marvin[3] = {127, -128, -128};
int8_t feat[1200];

static uint8_t __attribute__((aligned(16))) g_arena[ARENA_SIZE_BYTES];

bool validate_output(TfLiteTensor *tensor, int8_t *expected) {
    bool correct = true;
    for (int i = 0; i < 3; i++) {
        if (tensor->data.int8[i] != expected[i]) {
            correct = false;
        }
    }
    return correct;
}

extern "C" int benchmarkModel() {
    /* Parse the FlatBuffer */
    const tflite::Model *model =
        tflite::GetModel(output_str_ww_ref_model_vela_tflite);
    if (model->version() != TFLITE_SCHEMA_VERSION) {
        printf("tflm_npu_run: schema mismatch (%u != %u)\r\n",
               (unsigned)model->version(), (unsigned)TFLITE_SCHEMA_VERSION);
        return -1;
    }

    /*
     * Op resolver with one slot: the EthosU custom op.
     * Vela-compiled models contain a single EthosU operator that offloads
     * the entire graph to the NPU command stream.
     * AddEthosU() is a no-op if Register_ETHOSU() returns nullptr (stub build).
     */
    tflite::MicroMutableOpResolver<1> resolver;
    if (resolver.AddEthosU() != kTfLiteOk) {
        printf("tflm_npu_run: AddEthosU() failed\r\n");
        return -2;
    }

    tflite::MicroInterpreter interpreter(model, resolver, g_arena,
                                         ARENA_SIZE_BYTES);

    /*
     * Build the interpreter on the stack; all internal allocations come
     * from 'arena' so the object itself stays small.
     */
    if (interpreter.AllocateTensors() != kTfLiteOk) {
        printf(
            "tflm_npu_run: AllocateTensors() failed - arena too small "
            "(provided %u bytes)\r\n",
            (unsigned)ARENA_SIZE_BYTES);
        return -3;
    }

    printf("size: %d\n", interpreter.arena_used_bytes());

    TfLiteTensor *input = interpreter.input(0);

    /* Run inference - dispatches to the NPU via ethosu_invoke_v3 */
    int status = 0;
    __effcc_enter_profile_region("model");
    for (int iter = 0; iter < __effcc_get_profile_iterations(); iter++) {
        featurize(input_data_marvin, input->data.int8);
        status = interpreter.Invoke();
    }
    __effcc_exit_profile_region();

    if (status != kTfLiteOk) {
        printf("tflm_npu_run: Invoke() failed\r\n");
        return -4;
    }

    TfLiteTensor *output = interpreter.output(0);
    if (validate_output(output, expected_output_marvin)) {
        printf("[mlperftiny_sww] PASS\n");
    } else {
        printf("[mlperftiny_sww] FAIL\n");
        return 1;
    }
    return 0;
}
