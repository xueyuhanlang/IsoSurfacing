#pragma once

#include <cstdio>

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60

template <typename Real>
void print_progress_bar(Real percentage, bool verbose = true)
{
	if (verbose)
	{
		int val = (int)(percentage * 100);
		int lpad = (int)(percentage * PBWIDTH);
		int rpad = PBWIDTH - lpad;
		printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
		fflush(stdout);
	}
}
