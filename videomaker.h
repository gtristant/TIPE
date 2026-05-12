#ifndef VIDEO_H
#define VIDEO_H

#include "utils.h"

#define MAX_DEGRE 1000 //en Kelvin
#define MAX_DENSITY 1000 // en jsp

enum data_type {TEMPERATURE, DENSITY, VELOCITY};

typedef enum data_type data_type;

void save_frame(int frame_num, data_type d, double data[GRID_HEIGHT+1][GRID_WIDTH+1], double rho[GRID_HEIGHT+1][GRID_WIDTH+1],double ux[GRID_HEIGHT+1][GRID_WIDTH+1],double uy[GRID_HEIGHT+1][GRID_WIDTH+1]);
void save_frame_mv(int frame_num, data_type d, double data[GRID_HEIGHT+1][GRID_WIDTH+1], double rho[GRID_HEIGHT+1][GRID_WIDTH+1],double ux[GRID_HEIGHT+1][GRID_WIDTH+1],double uy[GRID_HEIGHT+1][GRID_WIDTH+1],int lower_limit);


#endif //VIDEO_H