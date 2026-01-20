#include <eff.h>
#include <eff/profile.h>
#include <stdlib.h>
#include <stdio.h>
#include "image.h"

int stbi_write_jpg(int x, int y, int comp, const void  *data, int quality);

extern unsigned char jpeg_output_arr[16 * 1024];
extern int output_arr_length;

const int answer_length = 7556;

int main() {
    START_PROFILE_REGION("kernel");
    for (int iter = 0; iter < 10; iter++) {
        output_arr_length = 0;
        stbi_write_jpg(IMAGE_WIDTH, IMAGE_HEIGHT, 3, src_image, 50);
    }
    END_PROFILE_REGION();

    // See your image!
    #ifdef EFF_BLD_SIM
    FILE* f = fopen("test.jpg", "wb");
    fwrite(jpeg_output_arr, 1, output_arr_length, f);
    fflush(f);
    fclose(f);
    #endif

    // comparing the result - allowed to be up to 2% off (float vs fixed point)
    if (abs(output_arr_length - answer_length) > answer_length / 50) {
        printf("[jpeg] FAIL: Encoded jpeg doesn't match answer's size. Got %d expected %d\r\n", output_arr_length, answer_length);
        return 1;
    }

    printf("[jpeg] PASS\r\n");
    return 0;
}