#include "hmm.h"

#include <math.h>

#define INV_SQRT_2PI 0.3989422804014327


double gaussian_pdf(double x, double m,double s)
{
        double a = (x-m) / s;
        return INV_SQRT_2PI / s *exp(-0.5 * a * a);
}
