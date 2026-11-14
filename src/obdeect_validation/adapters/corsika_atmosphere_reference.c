/* Small command-line adapter for the CORSIKA IACT/atmo reference implementation. */
#include "atmo.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* The tabulated profile path does not use CORSIKA's five-layer density routines. */
double rhof_(double *height) { (void)height; return 0.; }
double thick_(double *height) { (void)height; return 0.; }
double heigh_(double *thickness) { (void)thickness; return 0.; }

int main(int argc, char **argv) {
  if (argc != 2) return 2;
  const double observation_level = 2158.0e2; /* cm a.s.l. */
  atmnam_(argv[1], NULL);
  int profile = 99;
  atmset_(&profile, (double *)&observation_level);
  const double heights_m[] = {3000., 4000., 6000., 8000., 12000., 20000.};
  const double zenith_deg[] = {0., 20., 40., 60.};
  for (size_t ih = 0; ih < sizeof(heights_m) / sizeof(heights_m[0]); ++ih) {
    for (size_t iz = 0; iz < sizeof(zenith_deg) / sizeof(zenith_deg[0]); ++iz) {
      const double theta = zenith_deg[iz] * 3.14159265358979323846 / 180.;
      double zem = heights_m[ih] * 100.;
      double u = sin(theta), v = 0., w = cos(theta);
      double dx = 0., dy = 0., dt = 0.;
      raybnd_(&zem, &u, &v, &w, &dx, &dy, &dt);
      printf("%.0f %.0f %.17g %.17g %.17g %.17g\n", heights_m[ih], zenith_deg[iz], u, w, dx,
             dt);
    }
  }
  return 0;
}
