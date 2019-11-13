#include "common.h"

double CURRENT_COMPARE_PRECISION = CONST_COMPARE_PRECISION;
streamsize CURRENT_DISPLAY_PRECISION = CONST_DISPLAY_PRECISION;

void set_compare_precision(double precision)
{
	CURRENT_COMPARE_PRECISION = precision;
}

double get_current_compare_precision()
{
	return CURRENT_COMPARE_PRECISION;
}

void reset_compare_precision()
{
	CURRENT_COMPARE_PRECISION = CONST_COMPARE_PRECISION;
}

void set_display_precision(streamsize precision)
{
	CURRENT_DISPLAY_PRECISION = precision;
}

streamsize get_current_display_precision()
{
	return CURRENT_DISPLAY_PRECISION;
}

void reset_display_precision()
{
	CURRENT_DISPLAY_PRECISION = CONST_DISPLAY_PRECISION;
}
