#include <eff.h>
#include <stdio.h>
#include <stdlib.h>
#include "image.h"
#include "jpeg.h"

extern unsigned char jpeg_output_arr[16 * 1024];
extern int output_arr_length;

#define IMAGE_DIVIDER 10
const int answer_length = 1056;

int main() {
    stbi_write_jpg(IMAGE_WIDTH, IMAGE_HEIGHT / IMAGE_DIVIDER, 3, src_image,
                    50);

    // comparing the result - allowed to be up to 2% off (float vs fixed point)
    if (abs(output_arr_length - answer_length) > answer_length / 50) {
        printf(
            "[jpeg] FAIL: Encoded jpeg doesn't match answer's size. Got %d "
            "expected %d\r\n",
            output_arr_length, answer_length);
        return 1;
    }

    printf("[jpeg] PASS\r\n");
    return 0;
}
