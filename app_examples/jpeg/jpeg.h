int stbi_write_jpg(int x, int y, int comp, const void *data, int quality);
typedef void stbi_write_func(void *context, void *data, int size);

typedef struct
{
    stbi_write_func *func;
    void *context;
    unsigned char buffer[64];
    int buf_used;
} stbi__write_context;
