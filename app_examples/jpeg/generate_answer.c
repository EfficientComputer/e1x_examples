// <<AUTOBENCH>> skip
// to compile, run `gcc -O3 -Wall generate_answer.c image.c -o generate_answer`

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "jpeg_ref.h"

#define FIO_CHECK(x)                 \
    if ((x) < 0)                     \
    {                                \
        printf("File IO failed.\n"); \
        return 1;                    \
    }

int main()
{
    stbi_write_jpg("answer.jpg", IMAGE_WIDTH, IMAGE_HEIGHT / 10, 3, src_image,
                   50);

    // Creating answer C file
    FILE *jpg = fopen("answer.jpg", "rb");
    FILE *outC = fopen("answer.c.inc", "w");
    if (!jpg || !outC)
    {
        printf("failed to open files.\n");
        return 1;
    }

    FIO_CHECK(
        fputs("#include <stdint.h>\n\nuint8_t encoded_jpg_arr[] = {\n", outC));

    int counter = 0;
    int length = 0;

    while (1)
    {
        uint8_t val;
        if (!fread(&val, 1, 1, jpg))
        {
            break;
        }

        char buf[16];
        sprintf(buf, "0x%02X, ", val);

        FIO_CHECK(fputs(buf, outC));
        length++;

        if (++counter >= 16)
        {
            FIO_CHECK(fputc('\n', outC));
            counter = 0;
        }
    }

    FIO_CHECK(fputs("};\n\nuint8_t* encoded_jpg = &encoded_jpg_arr[0];", outC));

    FIO_CHECK(fprintf(outC, "\nint answer_length = %d;\n", length));

    FIO_CHECK(fflush(outC));
    FIO_CHECK(fclose(outC));
    FIO_CHECK(fclose(jpg));
    return 0;
}
