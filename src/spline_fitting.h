//
//  spline_fitting.h
//  tagdust2
//
//  Created by lassmann on 5/28/13.
//  Copyright (c) 2013 lassmann. All rights reserved.
//

#ifndef tagdust2_spline_fitting_h
#define tagdust2_spline_fitting_h


double seval (int n, double u, double x[], double y[],double b[], double c[], double d[],    int *last);

int spline (int n, int end1, int end2, double slope1, double slope2, double x[], double y[],  double b[], double c[], double d[], int *iflag);
void test_spline(void);

#endif
