#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

bool doesFileExist(const char *fname) {
    FILE *file;
    if ((file = fopen(fname, "r")))
    {
        fclose(file);
        return true;
    }
    return false;
}
