/* CORSIKA IACT/atmo raybnd() for the high-zenith refraction case.
 *
 * usage: corsika-refraction-high-zenith PROFILE OBS_LEVEL_M < rays
 * Each input line: "emission_altitude_cm u w" with u = sin(theta), w = cos(theta) (downward, > 0).
 * Each output line: "u_out w_out dx_cm dt_ns" (dx and dt start at 0, as in the existing adapter).
 */
#include "atmo.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* The tabulated profile path does not use CORSIKA's five-layer density routines. */
double rhof_(double *height) { (void)height; return 0.; }
double thick_(double *height) { (void)height; return 0.; }
double heigh_(double *thickness) { (void)thickness; return 0.; }

int main(int argc, char **argv) {
  if (argc != 3) return 2;
  double observation_level = atof(argv[2]) * 100.0; /* cm a.s.l. */
  atmnam_(argv[1], NULL);
  int profile = 99;
  atmset_(&profile, &observation_level);
  double zem, u, w, v;
  while (scanf("%lf %lf %lf", &zem, &u, &w) == 3) {
    v = 0.;
    double dx = 0., dy = 0., dt = 0.;
    raybnd_(&zem, &u, &v, &w, &dx, &dy, &dt);
    printf("%.17g %.17g %.17g %.17g\n", u, w, dx, dt);
  }
  return 0;
}
