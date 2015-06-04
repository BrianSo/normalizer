/*
 * Copyright (c) 2014, Blind Motion Project 
 * All rights reserved.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <string.h>

#include "types.hpp"
#include "plot.hpp"
#include "utils.hpp"
#include "config.hpp"

using namespace std;

vector<vector<string>> read_table(string const &filename) { // 2d array of string
    vector<vector<string>> table;
    ifstream input(filename);
    string line;
    int i = -1;
    while (getline(input, line)) {
        table.push_back(vector<string>());
        ++i;
        istringstream ss(line);
        string token;
        while (getline(ss, token, ';')) {
            table[i].push_back(token);
        }
    }
    return table;
}

void parse_data(vector<vector<string>> const &table,
        list &ta, list &xa, list &ya, list &za,
        list &tg, list &xg, list &yg, list &zg,
        list &t_geo, list &speed_geo) {
    for (int i = 0; i < table.size(); ++i) {
        if (table[i].size() == 0) {
            continue;
        }
        string type = table[i][0];
        if (type == "1") {
            double t = atof(table[i][2].c_str());
            ta.push_back(t);							//accelerometer excel time
            xa.push_back(atof(table[i][3].c_str()));	//accelerometer x
            ya.push_back(atof(table[i][4].c_str()));	//accelerometer y
            za.push_back(atof(table[i][5].c_str()));	//accelerometer z
        } else if (type == "4") {
            double t = atof(table[i][2].c_str());
            tg.push_back(t);							//gyroscope excel time
            xg.push_back(atof(table[i][3].c_str()));	//gyroscope x
            yg.push_back(atof(table[i][4].c_str()));	//gyroscope y
            zg.push_back(atof(table[i][5].c_str()));	//gyroscope z
        } else if (type == "geo") {
            t_geo.push_back(atof(table[i][2].c_str()));			//gps time
            speed_geo.push_back(atof(table[i][9].c_str()));		//gps speed
        }
    }
}

void replace_data(vector<vector<string>> &table,
        list const &ta, list const &xa,
        list const &ya, list const &za,
        list const &tg, list const &xg,
        list const &yg, list const &zg) {
    int ia = 0;
    int ig = 0;
    for (int i = 0; i < table.size(); ++i) {
        if (table[i].size() == 0) {
            continue;
        }
        string type = table[i][0];
        if (type == "1") {
            table[i][2] = to_string(ta[ia]);
            table[i][3] = to_string(xa[ia]);
            table[i][4] = to_string(ya[ia]);
            table[i][5] = to_string(za[ia]);
            ++ia;
        } else if (type == "4") {
            table[i][2] = to_string(tg[ig]);
            table[i][3] = to_string(xg[ig]);
            table[i][4] = to_string(yg[ig]);
            table[i][5] = to_string(zg[ig]);
            ++ig;
        }
    }
}

void write_data(string const &filename, vector<vector<string>> const &table) {
    ofstream output(filename);
    for (int i = 0; i < table.size(); ++i) {
        for (int j = 0; j < table[i].size(); ++j) {
            output << table[i][j] << ((j + 1 < table[i].size()) ? ';' : '\n');
        }
    }
    output.close();

}

void dumb_track_calculation(list const &ta, list const &xa, list const &ya,
        list const &tg, list const &zg, list &x, list &y, int start = 0) {
    double x_cur = 0;
    double y_cur = 0;
    double alpha_cur = 0;
    double vx_cur = 0;
    double vy_cur = 0;
    for (int i = 0; i < start; ++i) {
        x.push_back(0);
        y.push_back(0);
    }
    for (int i = start, j = 0; i < ta.size() - 1; ++i) {
        while (j < zg.size() - 1 && tg[j + 1] < ta[i]) {
            alpha_cur += zg[j] * (tg[j + 1] - tg[j]);
            cout << alpha_cur << ' ' << j << ' ' << zg[j] << ' ' << tg[j + 1] << ' ' << tg[j] << endl;
            ++j;
        }
        double dt = ta[i + 1] - ta[i];
        double ax_cur = xa[i] * cos(alpha_cur) + ya[i] * sin(alpha_cur);
        double ay_cur = -xa[i] * sin(alpha_cur) + ya[i] * cos(alpha_cur);
        vx_cur += ax_cur * dt;
        vy_cur += ay_cur * dt;
        x_cur += vx_cur * dt;
        y_cur += vy_cur * dt;
        x.push_back(x_cur);
        y.push_back(y_cur);
    }
}

// 2d array of string.
// vector is just a mutable array
// alternatives:
//	Java: arrayList
//	Objective-C: NSMutableArray
vector<vector<string>> table;	

// the type list is just a typedef of vector<double>
//
//		excel_time: example: 41899.486549016205	which indicating 11:40:37.835 at that day
//
list ta, xa, ya, za;			// data for accelerometer:	excel_time, x, y, z
list tg, xg, yg, zg;			// data for gyroscope:		excel_time, x, y, z
list t_geo, speed_geo;			// data for gps:			excel_time, speed
list xa_mean, ya_mean, za_mean;	// mean data for accelerometer
string output_filename;

#ifdef PYPLOT
PyPlot plt;
#endif

int main(int argc, char **argv) {

	// gflags or google related things are not the core of the algorithm
	// all gflags related things are just about parsing command line arguments (i.e. the -al flags for command 'ls -al')

	// *IMPORTANT all default flags values are defined in config.cpp
    google::SetUsageMessage("Program normalizes sensors values as recording device is always in the same orientation");
    google::ParseCommandLineFlags(&argc, &argv, true);
    if (argc <= 1) {
        cerr << "No input file given" << endl;
        return -1;
    }

	//read data from disk to 2d string array
	// data format:
	//		entry;entry;entry;entry;entry.... \n
	//		entry;entry;entry;entry;entry.... \n
	//		...
    table = read_table(argv[1]);


	// ## WE MIGHT START OUR CONVERSION FROM HERE ##
	//read all data from 2d string array to lists
    parse_data(table, ta, xa, ya, za, tg, xg, yg, zg, t_geo, speed_geo);

	//FLAGS_xxx is just gflags things.
    output_filename = FLAGS_output.length() > 0 ? FLAGS_output : "norm_" + string(argv[1]);

	// calculates corresponding quantile_mean to each elements
	//		the result is also a list
	//		p.s. I don't know what quantile_mean means. If needed, see algorithm in utils.cpp
	// 
	// to_mean(list const &t, list const &x, list &res, double radius = FLAGS_sm_radius, double percent = FLAGS_sm_range_part);
	//		FLAGS_sm_radius
	//			"Radius in excel time, used in acceleration data smoothing."
	//			= 0.5 * 1.0 / EXCEL_SECOND ~= 0.00000578703				EXCEL_SECOND = 1.0 / (24 * 60 * 60)
	//		FLAGS_sm_range_part
	//			Shows how many values will be taken for smoothing, very small and large ones will be thrown aside.
	//			= 0.5
    to_mean(ta, xa, xa_mean);
    to_mean(ta, ya, ya_mean);
    to_mean(ta, za, za_mean);



	// divide the data into blocks.
	//		return an array of index which indicates the borders of blocks.
	//		each index is the starting index of the block.
	//
	// vector<int> get_block_indices(list const &t, list const &x, list const &y, list const &z,
	//        double threshold = FLAGS_block_diff_thres, double time_thres = FLAGS_block_time_thres,
	//        bool adjacent = FLAGS_adjacent);
	//
	//			FLAGS_block_diff_thres = 100.0
	//			FLAGS_block_time_thres = 3.0 * EXCEL_SECOND
	//			FLAGS_adjacent = false			uses non-adjacent algorithm
	//											I don't know the difference :(
	//
    vector<int> block_starts = get_block_indices(ta, xa_mean, ya_mean, za_mean);

    cout << block_starts.size() << " block(s) found" << endl;

	// for each block: do core normalization algorithm
    for (int i = 0; i < block_starts.size(); ++i) {
        int start = block_starts[i];
        int finish = i < block_starts.size() - 1 ? block_starts[i + 1] : (int) ta.size();

		// calculate the transformation matrix in 2d array format. To understand requires knowledge in linear algebra.
		// normally the matrix in used to rotate the data.
		//
		// for this matrix, rotate the accelerometer data to make z axis points to ground
        vector<vector<double>> rot_matrix = get_z_rotation_matrix(start, finish, xa_mean, ya_mean, za_mean);
		// uses the matrix to rotate the data
        rotate_block(start, finish, xa_mean, ya_mean, za_mean, rot_matrix);

		// get the block of gyroscope data also
		// lower_bound is used to match the time (as the time of gyroscope data and accelerometer data are different)
		// lower_bound: find the index of gyroscope time which is the right lower than but nearest to the time (ta[block_starts[i]])
		// lower_bound might be implement by ourself as it is part of c++ <algoritgm.h>
        int start2 = (int) (lower_bound(tg.begin(), tg.end(), ta[block_starts[i]]) - tg.begin());
        int finish2 = i < block_starts.size() - 1 ? (int) (lower_bound(tg.begin(), tg.end(), ta[block_starts[i + 1]]) - tg.begin()) :
                ((int) tg.size());
		// similarly, rotate the data of gyroscope using same rotation matrix
        rotate_block(start2, finish2, xg, yg, zg, rot_matrix);

		// calculate a matrix for making x axis points forward
        vector<vector<double>> rot_matrix2 = get_plane_rotation_matrix(start, finish, ta, xa_mean, ya_mean, tg, zg,
                t_geo, speed_geo);

		//rotations....
        rotate_block(start, finish, xa_mean, ya_mean, za_mean, rot_matrix2);
        rotate_block(start2, finish2, xg, yg, zg, rot_matrix2);

		// note that the first line uses the first rot_matrix
		//  this should be able to move right below rotation of ?a_mean but for safe just leave here.
        rotate_block(start, finish, xa, ya, za, rot_matrix);
        rotate_block(start, finish, xa, ya, za, rot_matrix2);
    }
	// the ta, xa, ta, za, tg, xg, yg, zg lists are now normalized
	// ## WE MIGHT END OUR CONVERSION HERE ##

	// the remaining is just the actions of writing data to files


	// placing data back to the table(2d array of string), other sensor data unchanged
    replace_data(table, ta, xa, ya, za, tg, xg, yg, zg);
	// write to file
    write_data(output_filename, table);

	// the following is just plotting data
#ifdef PYPLOT
//    list x, y;
//    dumb_track_calculation(ta, xa_mean, ya_mean, tg, zg, x, y, 1598); // 1598 to skip big pause in 2014-09-28_SensorDatafile
//    plt.plot(x, y);

    plt.plot(ta, xa_mean);
    plt.plot(ta, ya_mean);
    plt.plot(ta, za_mean);
//    plt.plot(ta, za);
//    plt.plot(ta, zg);

    list vert;
    vert.push_back(-10);
    vert.push_back(10);
    for (int i = 1; i < block_starts.size(); ++i) {
        list hor;
        hor.push_back(ta[block_starts[i]]);
        hor.push_back(ta[block_starts[i]]);
        plt.plot(hor, vert, ", c='black', lw=3");
    }

    plt.show();
#endif

    return 0;
}
