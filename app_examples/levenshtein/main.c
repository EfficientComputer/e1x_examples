#include <stdio.h>
#include <stdlib.h>

#include "levenshtein.h"

#define NUM_ITERATIONS 1

char *s1 = "Efficient Launches E1x today";
int len1 = 28;
char *s2 =
    "Efficent lanched E1x today!"; // 2 insertions, 2 substitution, 1 deletion
int len2 = 27;

#define EDIT_DISTANCE 5

int main()
{
    int res;

    char *dynS1 = (char *)malloc(len1);
    char *dynS2 = (char *)malloc(len2);
    for (int i = 0; i < len1; i++)
    {
        dynS1[i] = s1[i];
    }
    for (int i = 0; i < len2; i++)
    {
        dynS2[i] = s2[i];
    }

    for (int iter = 0; iter < NUM_ITERATIONS; iter++)
        levenshtein(dynS1, len1, dynS2, len2, &res);

    if (res == EDIT_DISTANCE)
    {
        printf("[levenshtein] PASS\n");
    }
    else
    {
        printf("[levenshtein] FAIL: %d\n", res);
        return 1;
    }

    free(dynS1);
    free(dynS2);

    return 0;
}
