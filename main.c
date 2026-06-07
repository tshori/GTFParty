#include <stdio.h>
#include <stdlib.h>

#include "gtf.h"

int main(int argc, char *argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: %s <file.gtf>\n", argv[0]);
        return 1;
    }

    FILE *fp = fopen(argv[1], "r");
    if (!fp) {
        perror(argv[1]);   /* prints: "filename: <system error message>" */
        return 1;
    }

    char line[65536];
    GtfRecord rec;
    int parsed = 0, skipped = 0, errors = 0;

    while (fgets(line, sizeof(line), fp)) {
        int ret = gtf_parse_line(line, &rec);
        if (ret == 1) { skipped++; continue; }
        if (ret == -1) { errors++; continue; }

        gtf_print(&rec);
        parsed++;
    }

    fclose(fp);

    /* Always write diagnostics to stderr so they don't pollute redirected output */
    fprintf(stderr, "parsed=%d  skipped=%d  errors=%d\n",
            parsed, skipped, errors);
    return errors > 0 ? 1 : 0;
}
