/* Based on stb_image_write - v1.16 - public domain - http://nothings.org/stb
   writes out PNG/BMP/TGA/JPEG/HDR images to C stdio - Sean Barrett 2010-2015
   

   Modified by Daniel Cooper of Efficient Computer for application demonstration
*/

#include <eff.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#ifndef STBIW_ASSERT
#include <assert.h>
#define STBIW_ASSERT(x) assert(x)
#endif

#define STBIW_UCHAR(x) (unsigned char) ((x) & 0xff)

#define FIXED_POINT 32

typedef int32_t fixed_point;
#define FRAC_BITS 16
#define MUL_TYPE int64_t


// assumes x is quantized such that real_x = x / (1 << FRAC_BITS)
#define FIXED_POINT_MUL_CONST(x,real_y) ((fixed_point) (((MUL_TYPE) (x) * (MUL_TYPE) ((fixed_point)((real_y) * (float)(1<<FRAC_BITS)))) >> FRAC_BITS))
#define FIXED_POINT_MUL(x,y) ((fixed_point) (((MUL_TYPE) (x) * (MUL_TYPE) (y)) >> FRAC_BITS))

// -3 <= x <= 3
#define FIXED_POINT_DIV_CONST_SMALL_DIVIDEND(real_x,y) ((fixed_point) (((fixed_point) (real_x) << (FRAC_BITS * 2 - 3)) / (y)) << 3)
#define FIXED_POINT_LITERAL(real_x) ((fixed_point) ((real_x) * (float)(1 << FRAC_BITS)))

#define FIXED_POINT_ROUND_TO_INT(x) ((int) (((x) + (1 << (FRAC_BITS - 1))) >> FRAC_BITS))
#define FIXED_POINT_CAST_FROM_INT(x) ((fixed_point) (((fixed_point)(x)) << FRAC_BITS))

typedef void stbi_write_func(void *context, void *data, int size);

typedef struct
{
   stbi_write_func *func;
   void *context;
   unsigned char buffer[64];
   int buf_used;
} stbi__write_context;

// initialize a callback-based context
static void stbi__start_write_callbacks(stbi__write_context *s, stbi_write_func *c, void *context)
{
   s->func    = c;
   s->context = context;
}

#define KB *1024
unsigned char jpeg_output_arr[16 KB];
int output_arr_length = 0;
const int stbi__flip_vertically_on_write = 0;

int counter = 0;
static void stbi__stdio_write(void *context, void *data, int size)
{
    counter += size;
    for (int i = 0; i < size; i++) {
        jpeg_output_arr[output_arr_length++] = ((unsigned char*)data)[i];
    }
}

static void stbiw__fopen()
{
    // no file IO, just static arrays
    output_arr_length = 0;
}

static int stbi__start_write_file(stbi__write_context *s)
{
   stbiw__fopen();
   stbi__start_write_callbacks(s, stbi__stdio_write, NULL);
   return 1;
}

static void stbi__end_write_file(stbi__write_context *s)
{
   // no file IO, just static arrays
}

static void stbiw__putc(stbi__write_context *s, unsigned char c)
{
   counter += 1;
   jpeg_output_arr[output_arr_length++] = c;
}

typedef unsigned int stbiw_uint32;
typedef int stb_image_write_test[sizeof(stbiw_uint32)==4 ? 1 : -1];


/* ***************************************************************************
 *
 * JPEG writer
 *
 * This is based on Jon Olick's jo_jpeg.cpp:
 * public domain Simple, Minimalistic JPEG writer - http://www.jonolick.com/code.html
 */

static const unsigned char stbiw__jpg_ZigZag[] = { 0,1,5,6,14,15,27,28,2,4,7,13,16,26,29,42,3,8,12,17,25,30,41,43,9,11,18,
      24,31,40,44,53,10,19,23,32,39,45,52,54,20,22,33,38,46,51,55,60,21,34,37,47,50,56,59,61,35,36,48,49,57,58,62,63 };

static unsigned char stbiw__jpg_ZigZag_inv[64];

typedef struct {
   unsigned short bitBuf;
   unsigned short count;
} stbiw_bits;

static int stbiw__jpg_writeBits_direct(stbi__write_context *s, int *bitBufP, int *bitCntP, stbiw_bits bs, int pos) {
   int bitBuf = *bitBufP, bitCnt = *bitCntP;
   bitCnt += bs.count;
   bitBuf |= bs.bitBuf << (24 - bitCnt);
   int itCount = bitCnt >> 3;
   for(int i = 0; i < itCount; i++) {
      unsigned char c = (bitBuf >> 16) & 255;
      jpeg_output_arr[pos++] = c;
      if(c == 255) {
         jpeg_output_arr[pos++] = 0;
      }
      bitBuf <<= 8;
      bitCnt -= 8;
   }
   *bitBufP = bitBuf;
   *bitCntP = bitCnt;

   return pos;
}

__attribute__((always_inline)) static void stbiw__jpg_DCT(fixed_point *d0p, fixed_point *d1p, fixed_point *d2p, fixed_point *d3p, fixed_point *d4p, fixed_point *d5p, fixed_point *d6p, fixed_point *d7p) {
   fixed_point d0 = *d0p, d1 = *d1p, d2 = *d2p, d3 = *d3p, d4 = *d4p, d5 = *d5p, d6 = *d6p, d7 = *d7p;
   fixed_point z1, z2, z3, z4, z5, z11, z13;

   fixed_point tmp0 = d0 + d7;
   fixed_point tmp7 = d0 - d7;
   fixed_point tmp1 = d1 + d6;
   fixed_point tmp6 = d1 - d6;
   fixed_point tmp2 = d2 + d5;
   fixed_point tmp5 = d2 - d5;
   fixed_point tmp3 = d3 + d4;
   fixed_point tmp4 = d3 - d4;

   // Even part
   fixed_point tmp10 = tmp0 + tmp3;   // phase 2
   fixed_point tmp13 = tmp0 - tmp3;
   fixed_point tmp11 = tmp1 + tmp2;
   fixed_point tmp12 = tmp1 - tmp2;

   d0 = tmp10 + tmp11;       // phase 3
   d4 = tmp10 - tmp11;

   z1 = FIXED_POINT_MUL_CONST((tmp12 + tmp13), 0.707106781f); // c4
   d2 = tmp13 + z1;       // phase 5
   d6 = tmp13 - z1;

   // Odd part
   tmp10 = tmp4 + tmp5;       // phase 2
   tmp11 = tmp5 + tmp6;
   tmp12 = tmp6 + tmp7;

   // The rotator is modified from fig 4-8 to avoid extra negations.
   z5 = FIXED_POINT_MUL_CONST((tmp10 - tmp12), 0.382683433f); // c6
   z2 = FIXED_POINT_MUL_CONST(tmp10, 0.541196100f) + z5; // c2-c6
   z4 = FIXED_POINT_MUL_CONST(tmp12, 1.306562965f) + z5; // c2+c6
   z3 = FIXED_POINT_MUL_CONST(tmp11, 0.707106781f); // c4

   z11 = tmp7 + z3;      // phase 5
   z13 = tmp7 - z3;

   *d5p = z13 + z2;         // phase 6
   *d3p = z13 - z2;
   *d1p = z11 + z4;
   *d7p = z11 - z4;

   *d0p = d0;  *d2p = d2;  *d4p = d4;  *d6p = d6;
}

__attribute__((always_inline)) stbiw_bits stbiw__jpg_calcBits(int val) {
   int absVal = val < 0 ? -val : val;
   val = val < 0 ? val-1 : val;
   int bitCount = 1;
   while(absVal >>= 1) {
      bitCount++;
   }

   stbiw_bits ret = { val & ((1<<bitCount)-1), bitCount };
   return ret;
}

__attribute__((always_inline)) int stbiw__jpg_dct_DU(int32_t DU[64], fixed_point *CDU, int du_stride, fixed_point *fdtbl) {
   int lastNonZero0 = 0;
   int lastNonZero1 = 0;
   int lastNonZero2 = 0;
   int lastNonZero3 = 0;

   for (int outIdx = 0; outIdx < 64; outIdx+=4) {
      int j0 = stbiw__jpg_ZigZag_inv[outIdx];
      int i0 = (j0 & 0x7) + (j0 >> 3) * du_stride;
      int j1 = stbiw__jpg_ZigZag_inv[outIdx + 1];
      int i1 = (j1 & 0x7) + (j1 >> 3) * du_stride;
      int j2 = stbiw__jpg_ZigZag_inv[outIdx + 2];
      int i2 = (j2 & 0x7) + (j2 >> 3) * du_stride;
      int j3 = stbiw__jpg_ZigZag_inv[outIdx + 3];
      int i3 = (j3 & 0x7) + (j3 >> 3) * du_stride;

      fixed_point v0 = FIXED_POINT_MUL(CDU[i0], fdtbl[j0]);
      fixed_point v1 = FIXED_POINT_MUL(CDU[i1], fdtbl[j1]);
      fixed_point v2 = FIXED_POINT_MUL(CDU[i2], fdtbl[j2]);
      fixed_point v3 = FIXED_POINT_MUL(CDU[i3], fdtbl[j3]);
      int32_t outVal0 = FIXED_POINT_ROUND_TO_INT(v0);
      int32_t outVal1 = FIXED_POINT_ROUND_TO_INT(v1);
      int32_t outVal2 = FIXED_POINT_ROUND_TO_INT(v2);
      int32_t outVal3 = FIXED_POINT_ROUND_TO_INT(v3);
      DU[outIdx] = outVal0;
      DU[outIdx + 1] = outVal1;
      DU[outIdx + 2] = outVal2;
      DU[outIdx + 3] = outVal3;

      if (outVal0 != 0) {
         lastNonZero0 = outIdx;
      }
      if (outVal1 != 0) {
         lastNonZero1 = outIdx + 1;
      }
      if (outVal2 != 0) {
         lastNonZero2 = outIdx + 2;
      }
      if (outVal3 != 0) {
         lastNonZero3 = outIdx + 3;
      }
   }

   int last01 = lastNonZero0 > lastNonZero1 ? lastNonZero0 : lastNonZero1;
   int last23 = lastNonZero2 > lastNonZero3 ? lastNonZero2 : lastNonZero3;

   return last01 > last23 ? last01 : last23;
}

//__efficient__
void stbiw__jpg_encode_DC(stbi__write_context *s, int *bitBuf, int *bitCnt, int32_t DU[64], int DC, const unsigned short HTDC[256][2], const unsigned short HTAC[256][2], int* pos) {
   int bitBufLocal = *bitBuf;
   int bitCntLocal = *bitCnt;
   int posLocal = *pos;

   // Encode DC
   int diff = DU[0] - DC;
   int numDC = diff == 0 ? 1 : 2;

   stbiw_bits bitsDC = stbiw__jpg_calcBits(diff);
   
   for (int i = 0; numDC > i; i++) {
      stbiw_bits buf;
      if (numDC == 1) {
         buf.bitBuf = HTDC[0][0];
         buf.count = HTDC[0][1];
      } else if (i == 0) {
         buf.bitBuf = HTDC[bitsDC.count][0];
         buf.count = HTDC[bitsDC.count][1];
      } else {
         buf = bitsDC;
      }

      posLocal = stbiw__jpg_writeBits_direct(s, &bitBufLocal, &bitCntLocal, buf, posLocal);
   }

   *bitBuf = bitBufLocal;
   *bitCnt = bitCntLocal;
   *pos = posLocal;

   return;
}

//__efficient__
void stbiw__jpg_encode_AC(stbi__write_context *s, int *bitBuf, int *bitCnt, int32_t DU[64], int DC, const unsigned short HTDC[256][2], const unsigned short HTAC[256][2], int* pos, int end0pos) {
   const stbiw_bits M16zeroes = { HTAC[0xF0][0], HTAC[0xF0][1] };

   int bitBufLocal = *bitBuf;
   int bitCntLocal = *bitCnt;
   int posLocal = *pos;

   // Encode ACs
   
   int startpos = 1;
   for(int i = 1; i <= end0pos; ++i) {
      int64_t du_val = DU[i];
      if (du_val != 0) {
         stbiw_bits bits = stbiw__jpg_calcBits(du_val);

         int nrzeroes = i-startpos;
         int nrzeroesMasked = nrzeroes & 0xF;
         int lng = nrzeroes>>4;
         for (int nrmarker=0; nrmarker < lng + 2; ++nrmarker) {
            stbiw_bits buf;
            if (nrmarker < lng) {
               buf = M16zeroes;
            } else if (nrmarker == lng) {
               buf.bitBuf = HTAC[(nrzeroesMasked<<4)+bits.count][0];
               buf.count = HTAC[(nrzeroesMasked<<4)+bits.count][1];
            } else {
               buf = bits;
            }
            posLocal = stbiw__jpg_writeBits_direct(s, &bitBufLocal, &bitCntLocal, buf, posLocal);
         }

         startpos = i + 1;
      }
   }

   *bitBuf = bitBufLocal;
   *bitCnt = bitCntLocal;
   *pos = posLocal;
}

__attribute__((always_inline)) int stbiw__jpg_encode_DC_AC(stbi__write_context *s, int *bitBuf, int *bitCnt, int32_t DU[64], int DC, const unsigned short HTDC[256][2], const unsigned short HTAC[256][2], int* pos, int lastNonZero) {
   const stbiw_bits EOB = { HTAC[0x00][0], HTAC[0x00][1] };

   stbiw__jpg_encode_DC(s, bitBuf, bitCnt, DU, DC, HTDC, HTAC, pos);

   if (lastNonZero > 0) {
      stbiw__jpg_encode_AC(s, bitBuf, bitCnt, DU, DC, HTDC, HTAC, pos, lastNonZero);
   }

   if(lastNonZero == 0 || lastNonZero != 63) {
      *pos = stbiw__jpg_writeBits_direct(s, bitBuf, bitCnt, EOB, *pos);
   }

   return DU[0];
}

__efficient__
static void stbiw__subsample_UV(int subIterCount, fixed_point subUs[][64], fixed_point subVs[][64], fixed_point Us[][256], fixed_point Vs[][256]) {
   for (int j = 0; j < subIterCount; j++) {
      for(int yy = 0, pos = 0; yy < 8; ++yy) {
         for(int xx = 0; xx < 8; ++xx, ++pos) {
            int k = yy*32+xx*2;
            subUs[j][pos] = FIXED_POINT_MUL_CONST((Us[j][k+0] + Us[j][k+1] + Us[j][k+16] + Us[j][k+17]), 0.25f);
            subVs[j][pos] = FIXED_POINT_MUL_CONST((Vs[j][k+0] + Vs[j][k+1] + Vs[j][k+16] + Vs[j][k+17]), 0.25f);
         }
      }
   }
}

__efficient__
static void stbiw__convert_colors_16(int iterIdx,
                                  int subIterCount,
                                  int numIterW,
                                  int numIterH,
                                  int32_t magicNumber,
                                  int height,
                                  int width,
                                  int comp,
                                  int widthComp,
                                  const unsigned char *dataR,
                                  const unsigned char *dataG,
                                  const unsigned char *dataB,
                                  fixed_point Ys[][256],
                                  fixed_point Us[][256],
                                  fixed_point Vs[][256]) {
   for (int j = 0; j < subIterCount; j++) {
      // int x = ((j + iterIdx) % numIterW) * 16;
      // int y = ((j + iterIdx) / numIterW) * 16;
      int y = (((j + iterIdx) * magicNumber) >> 20) & 0xFFF0;
      int x = ((j + iterIdx)) * 16 - y * numIterW;

      __effcc_parallel(2) // <-- should be here!
      for(int row = y; row < y+16; row += 1) {
         // row >= height => use last input row
         int clamped_row = (row < height) ? row : height - 1;
         int base_p = (stbi__flip_vertically_on_write ? (height-1-clamped_row) : clamped_row)*widthComp;

         for(int col = x; col < x+16; ++col) {
            int pos = (row-y)*16 + (col - x);

            // if col >= width => use pixel from last input column
            int p = base_p + ((col < width) ? col : (width-1))*comp;
            fixed_point r = FIXED_POINT_CAST_FROM_INT(dataR[p]), g = FIXED_POINT_CAST_FROM_INT(dataG[p]), b = FIXED_POINT_CAST_FROM_INT(dataB[p]);
            Ys[j][pos] = FIXED_POINT_MUL_CONST(r, 0.29900f) + FIXED_POINT_MUL_CONST(g, 0.58700f) + FIXED_POINT_MUL_CONST(b, 0.11400f) - FIXED_POINT_LITERAL(128);
            Us[j][pos] = FIXED_POINT_MUL_CONST(r, -0.16874f) - FIXED_POINT_MUL_CONST(g, 0.33126f) + FIXED_POINT_MUL_CONST(b, 0.50000f);
            Vs[j][pos] = FIXED_POINT_MUL_CONST(r, 0.50000f) - FIXED_POINT_MUL_CONST(g, 0.41869f) - FIXED_POINT_MUL_CONST(b, 0.08131f);
         }
      }
   }
}

__efficient__
static void stbiw__convert_colors_8(int iterIdx,
                                  int subIterCount,
                                  int numIterW,
                                  int numIterH,
                                  int32_t magicNumber,
                                  int height,
                                  int width,
                                  int comp,
                                  int widthComp,
                                  const unsigned char *dataR,
                                  const unsigned char *dataG,
                                  const unsigned char *dataB,
                                  fixed_point Ys[][64],
                                  fixed_point Us[][64],
                                  fixed_point Vs[][64]) {
   for (int j = 0; j < subIterCount; j++) {
      // int x = ((j + iterIdx) % numIterW) * 8;
      // int y = ((j + iterIdx) / numIterW) * 8;
      int y = (((j + iterIdx) * magicNumber) >> 21) & 0xFFF8;
      int x = ((j + iterIdx)) * 8 - y * numIterW;

     // __effcc_parallel(2) // <-- should be here!
      for(int row = y; row < y+8; row += 1) {
         // row >= height => use last input row
         int clamped_row = (row < height) ? row : height - 1;
         int base_p = (stbi__flip_vertically_on_write ? (height-1-clamped_row) : clamped_row)*widthComp;

         for(int col = x; col < x+8; ++col) {
            int pos = (row-y)*8 + (col - x);

            // if col >= width => use pixel from last input column
            int p = base_p + ((col < width) ? col : (width-1))*comp;
            fixed_point r = FIXED_POINT_CAST_FROM_INT(dataR[p]), g = FIXED_POINT_CAST_FROM_INT(dataG[p]), b = FIXED_POINT_CAST_FROM_INT(dataB[p]);
            Ys[j][pos] = FIXED_POINT_MUL_CONST(r, 0.29900f) + FIXED_POINT_MUL_CONST(g, 0.58700f) + FIXED_POINT_MUL_CONST(b, 0.11400f) - FIXED_POINT_LITERAL(128);
            Us[j][pos] = FIXED_POINT_MUL_CONST(r, -0.16874f) - FIXED_POINT_MUL_CONST(g, 0.33126f) + FIXED_POINT_MUL_CONST(b, 0.50000f);
            Vs[j][pos] = FIXED_POINT_MUL_CONST(r, 0.50000f) - FIXED_POINT_MUL_CONST(g, 0.41869f) - FIXED_POINT_MUL_CONST(b, 0.08131f);
         }
      }
   }
}

static const int YQT[] = {16,11,10,16,24,40,51,61,12,12,14,19,26,58,60,55,14,13,16,24,40,57,69,56,14,17,22,29,51,87,80,62,18,22,
                             37,56,68,109,103,77,24,35,55,64,81,104,113,92,49,64,78,87,103,121,120,101,72,92,95,98,112,100,103,99};
static const int UVQT[] = {17,18,24,47,99,99,99,99,18,21,26,66,99,99,99,99,24,26,56,99,99,99,99,99,47,66,99,99,99,99,99,99,
                              99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99};

__efficient__
static void stbiw__jpg_prepare_UY_tables(int quality, unsigned char UVTable[64], unsigned char YTable[64]) {
   #ifdef EFF_BLD_HAND_OPTIMIZED
   #define DIV_100(x) (((x) * 2622) >> 18)
   #else
   #define DIV_100(x) ((x) / 100)
   #endif

   for(int i = 0; i < 64; i += 3) {
      int uvti, yti = DIV_100(YQT[i]*quality+50);
      YTable[stbiw__jpg_ZigZag[i]] = (unsigned char) (yti < 1 ? 1 : yti > 255 ? 255 : yti);
      uvti = DIV_100(UVQT[i]*quality+50);
      UVTable[stbiw__jpg_ZigZag[i]] = (unsigned char) (uvti < 1 ? 1 : uvti > 255 ? 255 : uvti);
   }

   for(int i = 1; i < 64; i += 3) {
      int uvti, yti = DIV_100(YQT[i]*quality+50);
      YTable[stbiw__jpg_ZigZag[i]] = (unsigned char) (yti < 1 ? 1 : yti > 255 ? 255 : yti);
      uvti = DIV_100(UVQT[i]*quality+50);
      UVTable[stbiw__jpg_ZigZag[i]] = (unsigned char) (uvti < 1 ? 1 : uvti > 255 ? 255 : uvti);
   }

   for(int i = 2; i < 64; i += 3) {
      int uvti, yti = DIV_100(YQT[i]*quality+50);
      YTable[stbiw__jpg_ZigZag[i]] = (unsigned char) (yti < 1 ? 1 : yti > 255 ? 255 : yti);
      uvti = DIV_100(UVQT[i]*quality+50);
      UVTable[stbiw__jpg_ZigZag[i]] = (unsigned char) (uvti < 1 ? 1 : uvti > 255 ? 255 : uvti);
   }
}

static const int32_t aasf_quant[] = {262143, 363604, 342507, 308248, 262143, 205965, 141871, 72325, 363604, 504333, 475071, 427553, 363604,
                                       285681, 196781, 100318, 342507, 475071, 447507, 402746, 342507, 269106, 185363, 94497, 308248, 427553,
                                       402746, 362462, 308248, 242189, 166823, 85045, 262143, 363604, 342507, 308248, 262143, 205965, 141871,
                                       72325, 205965, 285681, 269106, 242189, 205965, 161825, 111467, 56825, 141871, 196781, 185363, 166823,
                                       141871, 111467, 76780, 39142, 72325, 100318, 94497, 85045, 72325, 56825, 39142, 19954};

__efficient__
static void stbiw__jpg_prepare_fdtbl(unsigned char UVTable[64], unsigned char YTable[64], fixed_point fdtbl_Y[64], fixed_point fdtbl_UV[64], fixed_point zero) {
   for(int row = 0; row < 8; ++row) {
      for(int col = 0; col < 8; ++col) {
         int k = 8 * row + col;

         #define QUANT_FRAC 15

         #if QUANT_FRAC < FRAC_BITS
         #define CONST_SHIFT << (FRAC_BITS - QUANT_FRAC)
         #else
         #define CONST_SHIFT >> (QUANT_FRAC - FRAC_BITS)
         #endif

         fixed_point yVal = FIXED_POINT_CAST_FROM_INT(YTable [stbiw__jpg_ZigZag[k]]) + zero;
         fixed_point uvVal = FIXED_POINT_CAST_FROM_INT(UVTable[stbiw__jpg_ZigZag[k]]) + zero;

         fdtbl_Y[k]  = FIXED_POINT_DIV_CONST_SMALL_DIVIDEND(1, 
            FIXED_POINT_MUL(yVal, (int32_t)aasf_quant[row * 8 + col] CONST_SHIFT));
         fdtbl_UV[k] = FIXED_POINT_DIV_CONST_SMALL_DIVIDEND(1,
            FIXED_POINT_MUL(uvVal, (int32_t)aasf_quant[row * 8 + col] CONST_SHIFT));
      }
   }
}

__efficient__
static void stbiw__jpg_process_dct_rows_16(int subIterCount,
                                        fixed_point Ys[][256],
                                        fixed_point subUs[][64],
                                        fixed_point subVs[][64]) {
   for (int j = 0; j < subIterCount*6; j++) {
      int32_t div6 = ((uint32_t)j * 2796203) >> 24; // magic division!
      int32_t mod6 = j - 6 * div6;

      int duStride;
      fixed_point* CDU;

      if (mod6 < 4) {
         // Greyscale DUs
         CDU = Ys[div6] + (mod6 > 1 ? (mod6 == 2 ? 128 : 136) : (mod6 == 0 ? 0 : 8));
         duStride = 16;
      } else {
         // Color DUs
         CDU = mod6 == 4 ? subUs[div6] : subVs[div6];
         duStride = 8;
      }

      __effcc_parallel(1)
      for(int y = 0; y < duStride * 8; y+=duStride) {
         stbiw__jpg_DCT(CDU+y, CDU+y+1, CDU+y+2, CDU+y+3, CDU+y+4, CDU+y+5, CDU+y+6, CDU+y+7);
      }
   }
}

__efficient__
static void stbiw__jpg_process_dct_rows_8(int subIterCount,
                                        fixed_point Ys[][64],
                                        fixed_point subUs[][64],
                                        fixed_point subVs[][64]) {
   for (int j = 0; j < subIterCount*3; j++) {
      int32_t div3 = ((uint32_t)j * 5592406) >> 24; // magic division!
      int32_t mod3 = j - 3 * div3;

      int duStride = 8;
      fixed_point* CDU;

      CDU = mod3 == 0 ? Ys[div3] : (mod3 == 1 ? subUs[div3] : subVs[div3]);

      __effcc_parallel(1)
      for(int y = 0; y < duStride * 8; y+=duStride) {
         stbiw__jpg_DCT(CDU+y, CDU+y+1, CDU+y+2, CDU+y+3, CDU+y+4, CDU+y+5, CDU+y+6, CDU+y+7);
      }
   }
}

__efficient__
static void stbiw__jpg_process_dct_cols_16(int subIterCount,
                                        fixed_point Ys[][256],
                                        fixed_point subUs[][64],
                                        fixed_point subVs[][64]) {
   for (int j = 0; j < subIterCount*6; j++) {
      int32_t div6 = ((uint32_t)j * 2796203) >> 24; // magic division!
      int32_t mod6 = j - 6 * div6;

      int duStride;
      fixed_point* CDU;

      if (mod6 < 4) {
         // Greyscale DUs
         CDU = Ys[div6] + (mod6 > 1 ? (mod6 == 2 ? 128 : 136) : (mod6 == 0 ? 0 : 8));
         duStride = 16;
      } else {
         // Color DUs
         CDU = mod6 == 4 ? subUs[div6] : subVs[div6];
         duStride = 8;
      }

      __effcc_parallel(1)
      for(int x = 0; x < 8; x++) {
         stbiw__jpg_DCT(CDU+x, CDU+x+duStride, CDU+x+duStride*2, CDU+x+duStride*3, CDU+x+duStride*4, CDU+x+duStride*5, CDU+x+duStride*6, CDU+x+duStride*7);
      }
   }
}

__efficient__
static void stbiw__jpg_process_dct_cols_8(int subIterCount,
                                        fixed_point Ys[][64],
                                        fixed_point subUs[][64],
                                        fixed_point subVs[][64]) {
   for (int j = 0; j < subIterCount*3; j++) {
      int32_t div3 = ((uint32_t)j * 5592406) >> 24; // magic division!
      int32_t mod3 = j - 3 * div3;

      int duStride = 8;
      fixed_point* CDU;

      CDU = mod3 == 0 ? Ys[div3] : (mod3 == 1 ? subUs[div3] : subVs[div3]);

      __effcc_parallel(1)
      for(int x = 0; x < 8; x++) {
         stbiw__jpg_DCT(CDU+x, CDU+x+duStride, CDU+x+duStride*2, CDU+x+duStride*3, CDU+x+duStride*4, CDU+x+duStride*5, CDU+x+duStride*6, CDU+x+duStride*7);
      }
   }
}

__efficient__
void stbiw__jpg_complete_du_16(int subIterCount,
                                        fixed_point Ys[][256],
                                        fixed_point subUs[][64],
                                        fixed_point subVs[][64],
                                        fixed_point fdtbl_Y[64],
                                        fixed_point fdtbl_UV[64],
                                        int lastNonZeroIdxes[][6],
                                        int32_t DUs[][6][64]) {
   for (int j = 0; j < subIterCount*6; j++) {
      int32_t div6 = ((uint32_t)j * 2796203) >> 24; // magic division!
      int32_t mod6 = j - 6 * div6;

      int duStride;
      fixed_point* fdtbl;
      fixed_point* CDU;

      if (mod6 < 4) {
         // Greyscale DUs
         CDU = Ys[div6] + (mod6 > 1 ? (mod6 == 2 ? 128 : 136) : (mod6 == 0 ? 0 : 8));
         duStride = 16;
         fdtbl = fdtbl_Y;
      } else {
         // Color DUs
         CDU = mod6 == 4 ? subUs[div6] : subVs[div6];
         duStride = 8;
         fdtbl = fdtbl_UV;
      }

      int lastNonZero = stbiw__jpg_dct_DU(DUs[div6][mod6], CDU, duStride, fdtbl);
      lastNonZeroIdxes[div6][mod6] = lastNonZero;
   }
}

__efficient__
void stbiw__jpg_complete_du_8(int subIterCount,
                                        fixed_point Ys[][64],
                                        fixed_point subUs[][64],
                                        fixed_point subVs[][64],
                                        fixed_point fdtbl_Y[64],
                                        fixed_point fdtbl_UV[64],
                                        int lastNonZeroIdxes[][3],
                                        int32_t DUs[][3][64]) {
   for (int j = 0; j < subIterCount*3; j++) {
      int32_t div3 = ((uint32_t)j * 5592406) >> 24; // magic division!
      int32_t mod3 = j - 3 * div3;

      int duStride = 8;
      fixed_point* fdtbl;
      fixed_point* CDU;

      fdtbl = mod3 == 0 ? fdtbl_Y : fdtbl_UV;

      CDU = mod3 == 0 ? Ys[div3] : (mod3 == 1 ? subUs[div3] : subVs[div3]);

      int lastNonZero = stbiw__jpg_dct_DU(DUs[div3][mod3], CDU, duStride, fdtbl);
      lastNonZeroIdxes[div3][mod3] = lastNonZero;
   }
}

static const unsigned char std_dc_luminance_nrcodes[] = {0,0,1,5,1,1,1,1,1,1,0,0,0,0,0,0,0};
static const unsigned char std_dc_luminance_values[] = {0,1,2,3,4,5,6,7,8,9,10,11};
static const unsigned char std_ac_luminance_nrcodes[] = {0,0,2,1,3,3,2,4,3,5,5,4,4,0,0,1,0x7d};
static const unsigned char std_ac_luminance_values[] = {
   0x01,0x02,0x03,0x00,0x04,0x11,0x05,0x12,0x21,0x31,0x41,0x06,0x13,0x51,0x61,0x07,0x22,0x71,0x14,0x32,0x81,0x91,0xa1,0x08,
   0x23,0x42,0xb1,0xc1,0x15,0x52,0xd1,0xf0,0x24,0x33,0x62,0x72,0x82,0x09,0x0a,0x16,0x17,0x18,0x19,0x1a,0x25,0x26,0x27,0x28,
   0x29,0x2a,0x34,0x35,0x36,0x37,0x38,0x39,0x3a,0x43,0x44,0x45,0x46,0x47,0x48,0x49,0x4a,0x53,0x54,0x55,0x56,0x57,0x58,0x59,
   0x5a,0x63,0x64,0x65,0x66,0x67,0x68,0x69,0x6a,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7a,0x83,0x84,0x85,0x86,0x87,0x88,0x89,
   0x8a,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9a,0xa2,0xa3,0xa4,0xa5,0xa6,0xa7,0xa8,0xa9,0xaa,0xb2,0xb3,0xb4,0xb5,0xb6,
   0xb7,0xb8,0xb9,0xba,0xc2,0xc3,0xc4,0xc5,0xc6,0xc7,0xc8,0xc9,0xca,0xd2,0xd3,0xd4,0xd5,0xd6,0xd7,0xd8,0xd9,0xda,0xe1,0xe2,
   0xe3,0xe4,0xe5,0xe6,0xe7,0xe8,0xe9,0xea,0xf1,0xf2,0xf3,0xf4,0xf5,0xf6,0xf7,0xf8,0xf9,0xfa
};
static const unsigned char std_dc_chrominance_nrcodes[] = {0,0,3,1,1,1,1,1,1,1,1,1,0,0,0,0,0};
static const unsigned char std_dc_chrominance_values[] = {0,1,2,3,4,5,6,7,8,9,10,11};
static const unsigned char std_ac_chrominance_nrcodes[] = {0,0,2,1,2,4,4,3,4,7,5,4,4,0,1,2,0x77};
static const unsigned char std_ac_chrominance_values[] = {
   0x00,0x01,0x02,0x03,0x11,0x04,0x05,0x21,0x31,0x06,0x12,0x41,0x51,0x07,0x61,0x71,0x13,0x22,0x32,0x81,0x08,0x14,0x42,0x91,
   0xa1,0xb1,0xc1,0x09,0x23,0x33,0x52,0xf0,0x15,0x62,0x72,0xd1,0x0a,0x16,0x24,0x34,0xe1,0x25,0xf1,0x17,0x18,0x19,0x1a,0x26,
   0x27,0x28,0x29,0x2a,0x35,0x36,0x37,0x38,0x39,0x3a,0x43,0x44,0x45,0x46,0x47,0x48,0x49,0x4a,0x53,0x54,0x55,0x56,0x57,0x58,
   0x59,0x5a,0x63,0x64,0x65,0x66,0x67,0x68,0x69,0x6a,0x73,0x74,0x75,0x76,0x77,0x78,0x79,0x7a,0x82,0x83,0x84,0x85,0x86,0x87,
   0x88,0x89,0x8a,0x92,0x93,0x94,0x95,0x96,0x97,0x98,0x99,0x9a,0xa2,0xa3,0xa4,0xa5,0xa6,0xa7,0xa8,0xa9,0xaa,0xb2,0xb3,0xb4,
   0xb5,0xb6,0xb7,0xb8,0xb9,0xba,0xc2,0xc3,0xc4,0xc5,0xc6,0xc7,0xc8,0xc9,0xca,0xd2,0xd3,0xd4,0xd5,0xd6,0xd7,0xd8,0xd9,0xda,
   0xe2,0xe3,0xe4,0xe5,0xe6,0xe7,0xe8,0xe9,0xea,0xf2,0xf3,0xf4,0xf5,0xf6,0xf7,0xf8,0xf9,0xfa
};
// Huffman tables
static const unsigned short YDC_HT[256][2] = { {0,2},{2,3},{3,3},{4,3},{5,3},{6,3},{14,4},{30,5},{62,6},{126,7},{254,8},{510,9}};
static const unsigned short UVDC_HT[256][2] = { {0,2},{1,2},{2,2},{6,3},{14,4},{30,5},{62,6},{126,7},{254,8},{510,9},{1022,10},{2046,11}};
static const unsigned short YAC_HT[256][2] = {
   {10,4},{0,2},{1,2},{4,3},{11,4},{26,5},{120,7},{248,8},{1014,10},{65410,16},{65411,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {12,4},{27,5},{121,7},{502,9},{2038,11},{65412,16},{65413,16},{65414,16},{65415,16},{65416,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {28,5},{249,8},{1015,10},{4084,12},{65417,16},{65418,16},{65419,16},{65420,16},{65421,16},{65422,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {58,6},{503,9},{4085,12},{65423,16},{65424,16},{65425,16},{65426,16},{65427,16},{65428,16},{65429,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {59,6},{1016,10},{65430,16},{65431,16},{65432,16},{65433,16},{65434,16},{65435,16},{65436,16},{65437,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {122,7},{2039,11},{65438,16},{65439,16},{65440,16},{65441,16},{65442,16},{65443,16},{65444,16},{65445,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {123,7},{4086,12},{65446,16},{65447,16},{65448,16},{65449,16},{65450,16},{65451,16},{65452,16},{65453,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {250,8},{4087,12},{65454,16},{65455,16},{65456,16},{65457,16},{65458,16},{65459,16},{65460,16},{65461,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {504,9},{32704,15},{65462,16},{65463,16},{65464,16},{65465,16},{65466,16},{65467,16},{65468,16},{65469,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {505,9},{65470,16},{65471,16},{65472,16},{65473,16},{65474,16},{65475,16},{65476,16},{65477,16},{65478,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {506,9},{65479,16},{65480,16},{65481,16},{65482,16},{65483,16},{65484,16},{65485,16},{65486,16},{65487,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {1017,10},{65488,16},{65489,16},{65490,16},{65491,16},{65492,16},{65493,16},{65494,16},{65495,16},{65496,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {1018,10},{65497,16},{65498,16},{65499,16},{65500,16},{65501,16},{65502,16},{65503,16},{65504,16},{65505,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {2040,11},{65506,16},{65507,16},{65508,16},{65509,16},{65510,16},{65511,16},{65512,16},{65513,16},{65514,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {65515,16},{65516,16},{65517,16},{65518,16},{65519,16},{65520,16},{65521,16},{65522,16},{65523,16},{65524,16},{0,0},{0,0},{0,0},{0,0},{0,0},
   {2041,11},{65525,16},{65526,16},{65527,16},{65528,16},{65529,16},{65530,16},{65531,16},{65532,16},{65533,16},{65534,16},{0,0},{0,0},{0,0},{0,0},{0,0}
};
static const unsigned short UVAC_HT[256][2] = {
   {0,2},{1,2},{4,3},{10,4},{24,5},{25,5},{56,6},{120,7},{500,9},{1014,10},{4084,12},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {11,4},{57,6},{246,8},{501,9},{2038,11},{4085,12},{65416,16},{65417,16},{65418,16},{65419,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {26,5},{247,8},{1015,10},{4086,12},{32706,15},{65420,16},{65421,16},{65422,16},{65423,16},{65424,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {27,5},{248,8},{1016,10},{4087,12},{65425,16},{65426,16},{65427,16},{65428,16},{65429,16},{65430,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {58,6},{502,9},{65431,16},{65432,16},{65433,16},{65434,16},{65435,16},{65436,16},{65437,16},{65438,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {59,6},{1017,10},{65439,16},{65440,16},{65441,16},{65442,16},{65443,16},{65444,16},{65445,16},{65446,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {121,7},{2039,11},{65447,16},{65448,16},{65449,16},{65450,16},{65451,16},{65452,16},{65453,16},{65454,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {122,7},{2040,11},{65455,16},{65456,16},{65457,16},{65458,16},{65459,16},{65460,16},{65461,16},{65462,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {249,8},{65463,16},{65464,16},{65465,16},{65466,16},{65467,16},{65468,16},{65469,16},{65470,16},{65471,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {503,9},{65472,16},{65473,16},{65474,16},{65475,16},{65476,16},{65477,16},{65478,16},{65479,16},{65480,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {504,9},{65481,16},{65482,16},{65483,16},{65484,16},{65485,16},{65486,16},{65487,16},{65488,16},{65489,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {505,9},{65490,16},{65491,16},{65492,16},{65493,16},{65494,16},{65495,16},{65496,16},{65497,16},{65498,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {506,9},{65499,16},{65500,16},{65501,16},{65502,16},{65503,16},{65504,16},{65505,16},{65506,16},{65507,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {2041,11},{65508,16},{65509,16},{65510,16},{65511,16},{65512,16},{65513,16},{65514,16},{65515,16},{65516,16},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},
   {16352,14},{65517,16},{65518,16},{65519,16},{65520,16},{65521,16},{65522,16},{65523,16},{65524,16},{65525,16},{0,0},{0,0},{0,0},{0,0},{0,0},
   {1018,10},{32707,15},{65526,16},{65527,16},{65528,16},{65529,16},{65530,16},{65531,16},{65532,16},{65533,16},{65534,16},{0,0},{0,0},{0,0},{0,0},{0,0}
};

static int stbi_write_jpg_core(stbi__write_context *s, int width, int height, int comp, const void* data, int quality, fixed_point zero) {
   int row, col, i, k, subsample;
   fixed_point fdtbl_Y[64], fdtbl_UV[64];
   unsigned char YTable[64], UVTable[64];

   if(!data || !width || !height || comp > 4 || comp < 1) {
      return 0;
   }

   quality = quality ? quality : 90;
   subsample = quality <= 90 ? 1 : 0;
   quality = quality < 1 ? 1 : quality > 100 ? 100 : quality;
   quality = quality < 50 ? 5000 / quality : 200 - quality * 2;

   stbiw__jpg_prepare_UY_tables(quality, UVTable, YTable);

   stbiw__jpg_prepare_fdtbl(UVTable, YTable, fdtbl_Y, fdtbl_UV, zero);

   for (int i = 0; i < 64; i++) {
      stbiw__jpg_ZigZag_inv[stbiw__jpg_ZigZag[i]] = i;
   }

   // Write Headers
   {
      static const unsigned char head0[] = { 0xFF,0xD8,0xFF,0xE0,0,0x10,'J','F','I','F',0,1,1,0,0,1,0,1,0,0,0xFF,0xDB,0,0x84,0 };
      static const unsigned char head2[] = { 0xFF,0xDA,0,0xC,3,1,0,2,0x11,3,0x11,0,0x3F,0 };
      const unsigned char head1[] = { 0xFF,0xC0,0,0x11,8,(unsigned char)(height>>8),STBIW_UCHAR(height),(unsigned char)(width>>8),STBIW_UCHAR(width),
                                      3,1,(unsigned char)(subsample?0x22:0x11),0,2,0x11,1,3,0x11,1,0xFF,0xC4,0x01,0xA2,0 };
      stbi__stdio_write(s->context, (void*)head0, sizeof(head0));
      stbi__stdio_write(s->context, (void*)YTable, sizeof(YTable));
      stbiw__putc(s, 1);
      stbi__stdio_write(s->context, UVTable, sizeof(UVTable));
      stbi__stdio_write(s->context, (void*)head1, sizeof(head1));
      stbi__stdio_write(s->context, (void*)(std_dc_luminance_nrcodes+1), sizeof(std_dc_luminance_nrcodes)-1);
      stbi__stdio_write(s->context, (void*)std_dc_luminance_values, sizeof(std_dc_luminance_values));
      stbiw__putc(s, 0x10); // HTYACinfo
      stbi__stdio_write(s->context, (void*)(std_ac_luminance_nrcodes+1), sizeof(std_ac_luminance_nrcodes)-1);
      stbi__stdio_write(s->context, (void*)std_ac_luminance_values, sizeof(std_ac_luminance_values));
      stbiw__putc(s, 1); // HTUDCinfo
      stbi__stdio_write(s->context, (void*)(std_dc_chrominance_nrcodes+1), sizeof(std_dc_chrominance_nrcodes)-1);
      stbi__stdio_write(s->context, (void*)std_dc_chrominance_values, sizeof(std_dc_chrominance_values));
      stbiw__putc(s, 0x11); // HTUACinfo
      stbi__stdio_write(s->context, (void*)(std_ac_chrominance_nrcodes+1), sizeof(std_ac_chrominance_nrcodes)-1);
      stbi__stdio_write(s->context, (void*)std_ac_chrominance_values, sizeof(std_ac_chrominance_values));
      stbi__stdio_write(s->context, (void*)head2, sizeof(head2));
   }

   // Encode 8x8 macroblocks
   {
      const stbiw_bits fillBits = {0x7F, 7};
      int DCY=0, DCU=0, DCV=0;
      int bitBuf=0, bitCnt=0;
      // comp == 2 is grey+alpha (alpha is ignored)
      int ofsG = comp > 2 ? 1 : 0, ofsB = comp > 2 ? 2 : 0;
      const unsigned char *dataR = (const unsigned char *)data;
      const unsigned char *dataG = dataR + ofsG;
      const unsigned char *dataB = dataR + ofsB;
      int pos;
      if(subsample) {
         #define BLOCK_SIZE 8 // may not be larger than 42
         int numIterW = ((width + 15) / 16);
         int numIterH = ((height + 15) / 16);
         int numIter = numIterW * numIterH;

         for (int i = 0; i < numIter; i+=BLOCK_SIZE) {

            static fixed_point Ys[BLOCK_SIZE][256];
            static fixed_point Us[BLOCK_SIZE][256];
            static fixed_point Vs[BLOCK_SIZE][256];
            static fixed_point subUs[BLOCK_SIZE][64];
            static fixed_point subVs[BLOCK_SIZE][64];
            static int32_t DUs[BLOCK_SIZE][6][64];
            static int lastNonZeroIdxes[BLOCK_SIZE][6];

            int subIterCount = i + BLOCK_SIZE >= numIter ? numIter - i : BLOCK_SIZE;

            int32_t magicNumber = ((1 << 24) + numIterW - 1) / numIterW;

            // Converting color space from RGB to YUV
            stbiw__convert_colors_16(i,
                                  subIterCount,
                                  numIterW,
                                  numIterH,
                                  magicNumber,
                                  height,
                                  width,
                                  comp,
                                  width * comp,
                                  dataR,
                                  dataG,
                                  dataB,
                                  Ys,
                                  Us,
                                  Vs);

            // Prepparing sub-sampled U, Vs
            stbiw__subsample_UV(subIterCount, subUs, subVs, Us, Vs);

            // Calculating DCTs
            stbiw__jpg_process_dct_rows_16(subIterCount, Ys, subUs, subVs);
            stbiw__jpg_process_dct_cols_16(subIterCount, Ys, subUs, subVs);

            stbiw__jpg_complete_du_16(subIterCount, Ys, subUs, subVs, fdtbl_Y, fdtbl_UV, lastNonZeroIdxes, DUs);

            int filePos = output_arr_length;

            for (int j = 0; j < subIterCount; j++) {
               DCY = stbiw__jpg_encode_DC_AC(s, &bitBuf, &bitCnt, DUs[j][0], DCY, YDC_HT, YAC_HT, &filePos, lastNonZeroIdxes[j][0]);
               DCY = stbiw__jpg_encode_DC_AC(s, &bitBuf, &bitCnt, DUs[j][1], DCY, YDC_HT, YAC_HT, &filePos, lastNonZeroIdxes[j][1]);
               DCY = stbiw__jpg_encode_DC_AC(s, &bitBuf, &bitCnt, DUs[j][2], DCY, YDC_HT, YAC_HT, &filePos, lastNonZeroIdxes[j][2]);
               DCY = stbiw__jpg_encode_DC_AC(s, &bitBuf, &bitCnt, DUs[j][3], DCY, YDC_HT, YAC_HT, &filePos, lastNonZeroIdxes[j][3]);
               DCU = stbiw__jpg_encode_DC_AC(s, &bitBuf, &bitCnt, DUs[j][4], DCU, UVDC_HT, UVAC_HT, &filePos, lastNonZeroIdxes[j][4]);
               DCV = stbiw__jpg_encode_DC_AC(s, &bitBuf, &bitCnt, DUs[j][5], DCV, UVDC_HT, UVAC_HT, &filePos, lastNonZeroIdxes[j][5]);
            }

            output_arr_length = filePos;
         }
      } else {
         #define BLOCK_SIZE 8 // may not be larger than 42
         int numIterW = ((width + 7) / 8);
         int numIterH = ((height + 7) / 8);
         int numIter = numIterW * numIterH;

         for (int i = 0; i < numIter; i+=BLOCK_SIZE) {

            static fixed_point Ys[BLOCK_SIZE][64];
            static fixed_point Us[BLOCK_SIZE][64];
            static fixed_point Vs[BLOCK_SIZE][64];
            static int32_t DUs[BLOCK_SIZE][3][64];
            static int lastNonZeroIdxes[BLOCK_SIZE][3];

            int subIterCount = i + BLOCK_SIZE >= numIter ? numIter - i : BLOCK_SIZE;

            int32_t magicNumber = ((1 << 24) + numIterW - 1) / numIterW;

            // Converting color space from RGB to YUV
            stbiw__convert_colors_8(i,
                                  subIterCount,
                                  numIterW,
                                  numIterH,
                                  magicNumber,
                                  height,
                                  width,
                                  comp,
                                  width*comp,
                                  dataR,
                                  dataG,
                                  dataB,
                                  Ys,
                                  Us,
                                  Vs);

            // Calculating DCTs
            stbiw__jpg_process_dct_rows_8(subIterCount, Ys, Us, Vs);
            stbiw__jpg_process_dct_cols_8(subIterCount, Ys, Us, Vs);

            stbiw__jpg_complete_du_8(subIterCount, Ys, Us, Vs, fdtbl_Y, fdtbl_UV, lastNonZeroIdxes, DUs);

            int filePos = output_arr_length;

            for (int j = 0; j < subIterCount; j++) {
               DCY = stbiw__jpg_encode_DC_AC(s, &bitBuf, &bitCnt, DUs[j][0], DCY, YDC_HT, YAC_HT, &filePos, lastNonZeroIdxes[j][0]);
               DCU = stbiw__jpg_encode_DC_AC(s, &bitBuf, &bitCnt, DUs[j][1], DCU, UVDC_HT, UVAC_HT, &filePos, lastNonZeroIdxes[j][1]);
               DCV = stbiw__jpg_encode_DC_AC(s, &bitBuf, &bitCnt, DUs[j][2], DCV, UVDC_HT, UVAC_HT, &filePos, lastNonZeroIdxes[j][2]);
            }

            output_arr_length = filePos;
         }
      }

      // Do the bit alignment of the EOI marker
      output_arr_length = stbiw__jpg_writeBits_direct(s, &bitBuf, &bitCnt, fillBits, output_arr_length);
   }

   // EOI
   stbiw__putc(s, 0xFF);
   stbiw__putc(s, 0xD9);

   return 1;
}

int stbi_write_jpg(int x, int y, int comp, const void *data, int quality)
{
   stbi__write_context s = { 0 };
   if (stbi__start_write_file(&s)) {
      int r = stbi_write_jpg_core(&s, x, y, comp, data, quality, 0);
      stbi__end_write_file(&s);
      return r;
   } else
      return 0;
}
