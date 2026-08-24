#include <stdio.h>
#include <string.h>

#include "iplib.h"

int
main(int argc, char **argv)
{
    const char *actual_version;
    const char *expected_version;

    if (argc != 2) {
        fprintf(stderr, "Usage: %s EXPECTED_VERSION\n", argv[0]);
        return 1;
    }

    expected_version = argv[1];
    actual_version = ip_get_version();

    if (actual_version == NULL) {
        fprintf(stderr, "ip_get_version() returned NULL\n");
        return 1;
    }

    if (strcmp(actual_version, expected_version) != 0) {
        fprintf(stderr, "Version mismatch: expected \"%s\", got \"%s\"\n",
                expected_version, actual_version);
        return 1;
    }

    printf("NCEPLIBS-ip version: %s\n", actual_version);

    return 0;
}
