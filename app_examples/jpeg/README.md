# JPEG Encode with Reference Validation (JPEG)

This example runs a JPEG codec over an embedded test image on the Electron E1x general-purpose processor. It compresses the image with an optimized encode path and validates the result against a trusted reference, a realistic image-processing workload with mixed control flow and arithmetic.

---

## 1. Overview

### What is JPEG?

JPEG is the most widely used lossy image compression format. Encoding an image works block by block: the pixels are converted from RGB into a luma and chroma color space, each 8x8 block is transformed with a discrete cosine transform (DCT), the transformed coefficients are quantized to discard visually unimportant detail, and the quantized values are entropy coded (Huffman) into a compact bitstream. This example runs an optimized fixed-point encoder over an embedded image and checks the output against a known-good reference.

### Mathematical Definition

The central step is the 2D DCT applied to each 8x8 pixel block, followed by quantization. For an 8x8 block of samples f(x, y), the transform coefficients F(u, v) are:

    F(u,v) = (1/4) * C(u) * C(v) * Σx Σy f(x,y) * cos((2x+1)*u*pi/16) * cos((2y+1)*v*pi/16)

where the sums run over x = 0..7 and y = 0..7, and C(w) = 1/sqrt(2) for w = 0 and 1 otherwise. Each coefficient is then quantized:

    Fq(u,v) = round( F(u,v) / Q(u,v) )

where Q(u,v) is the quantization table entry, scaled by the chosen quality setting. This example computes the DCT and quantization in fixed-point arithmetic.

---

## 2. Why This Kernel Matters

JPEG is deployed almost everywhere images are stored or moved, so its codec is a representative real-world workload:

- Cameras and imaging pipelines, where captured frames are compressed for storage or transmission
- Embedded and edge vision systems, where images are encoded before being sent upstream
- Machine learning preprocessing, where large image sets are decoded and re-encoded at scale
- Mobile and battery-powered devices, where image handling must stay within a tight energy budget

Because it combines block transforms, quantization, entropy coding, and irregular control flow, JPEG is a strong test of how well an architecture handles a mix of arithmetic and data-dependent branching.

---

## 3. Why EFF Hardware Performs Well

The Electron E1x runs programs on the Fabric architecture, a spatial dataflow design. Rather than repeatedly fetching, decoding, and scheduling instructions the way a traditional processor does, the effcc Compiler maps the kernel onto the Fabric as a dataflow graph. Operations fire as soon as their inputs are ready, and intermediate values flow directly between compute elements instead of moving through memory.

This fits JPEG encoding well:

- The 8x8 DCT is a fixed pattern of multiplies and adds that maps cleanly onto the dataflow graph and runs the same way for every block
- Block coefficients stay resident on the Fabric and flow directly from the transform into quantization and entropy coding instead of spilling to memory
- The regular, per-block pixel access pattern is known ahead of time, so addressing is built into the graph
- Fixed-point arithmetic keeps the transform and quantization efficient on integer hardware
- The steady stream of blocks overlaps computation with data movement, keeping throughput high while control overhead stays low

The result is efficient image compression at low energy, which is the metric that matters most for battery-powered and always-on devices.

---

## 4. Configurable Parameters

The test image is embedded in the app (a 240x160 RGB image in `image.c`, declared in `image.h`), so there is no external input file. These definitions and constants in `main.c` control the benchmark. Change them to re-run it or adjust what is encoded.

| Definition       | Default | Effect                                                                                                                                                                                                                                                                                                                                          |
| ---------------- | ------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `NUM_ITERATIONS` | `1`     | How many times the encode kernel runs. Increase it to average out measurement noise when benchmarking.                                                                                                                                                                                                                                          |
| `IMAGE_DIVIDER`  | `10`    | Divides the image height, so only the top `IMAGE_HEIGHT / IMAGE_DIVIDER` rows are encoded. Lowering it encodes more of the image and increases the work per run.                                                                                                                                                                                |
| `answer_length`  | `1056`  | The expected size in bytes of the encoded JPEG produced by the reference. The optimized output size is compared against this value with a 2 percent tolerance (to account for fixed-point versus floating-point rounding). If you change the image, the encoded region, or the quality, this must be updated to match the new reference result. |

The image dimensions are fixed by the embedded data: `IMAGE_WIDTH` (240) and `IMAGE_HEIGHT` (160) are defined in `image.h`. The encode quality is passed as an argument to the encoder (50 in this example). The reference size (`answer_length`) is produced offline by `generate_answer.c`, which runs the reference encoder on the same image. If you change the input or any encode setting, regenerate the reference so the comparison matches the new correct result.
