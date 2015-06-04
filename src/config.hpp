/*
 * Copyright (c) 2014, Blind Motion Project 
 * All rights reserved.
 */

#pragma once
#include <string>
#include "types.hpp"

//"Threshold for detecting when the phone was moved inside the vehicle. "
//"Shows maximal difference between vectors of acceleration when they are considered from the same interval."
static double FLAGS_block_diff_thres = 100.0;

//"Threshold for detecting when the phone was moved inside the vehicle. "
//"Shows maximal difference between time of neighboring vectors of acceleration "
//"when they are considered from the same interval."
static double FLAGS_block_time_thres = 3.0 * EXCEL_SECOND;

//"Type of the algorithm for detecting when the phone was moved inside the vehicle."
static bool FLAGS_adjacent = false;

//"Radius in excel time, used in acceleration data smoothing."
static double FLAGS_sm_radius = 0.5 * EXCEL_SECOND;

//"Shows how many values will be taken for smoothing, "
static double FLAGS_sm_range_part = 0.5;

//"Shows how many acceleration vectors will be taken for counting mean acceleration vector, "
//"very small and large ones will be thrown aside."
static double FLAGS_z_range_part = 0.3;

//"Sets maximal time length (in excel time) "
//"between two speed values for which speed derivative can be calculated."
static double FLAGS_speed_detection_thres = 3.0 * EXCEL_SECOND;

//"Output file name, if nothing passed, will be \"norm_<input>\""
static string FLAGS_output = "";

