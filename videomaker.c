#include "videomaker.h"

void save_frame(int frame_num, data_type d, double data[GRID_HEIGHT+1][GRID_WIDTH+1], double rho[GRID_HEIGHT+1][GRID_WIDTH+1],double uj[GRID_HEIGHT+1][GRID_WIDTH+1],double ui[GRID_HEIGHT+1][GRID_WIDTH+1]) {
    char filename[64];
    sprintf(filename, "therm_%04d.ppm", frame_num);
    FILE *f = fopen(filename, "w");
    if (!f) return;
    
    // Header PPM
    fprintf(f, "P6\n%d %d\n255\n", GRID_WIDTH+1, GRID_HEIGHT+1);

    unsigned char* img = (unsigned char*) calloc( (GRID_WIDTH+1) * (GRID_HEIGHT+1) * 3, sizeof(unsigned char));

    int r, g, b;
    if (d == TEMPERATURE){
        for(int i =0; i <= GRID_HEIGHT; i++){
            for(int j= 0; j <= GRID_WIDTH; j++){
                double v_carre = uj[i][j]*uj[i][j] + ui[i][j]*ui[i][j];
                double T = (data[i][j] - 0.5*v_carre)/CV_T_VAPEUR;
                get_temperature_color((float)T, (float)MAX_DEGRE, &r, &g, &b);
                int pj = (i* (GRID_WIDTH+1) + j) * 3;
                img[pj] = (unsigned char)r;
                img[pj + 1] = (unsigned char)g;
                img[pj + 2] = (unsigned char)b;  
            }
        }
    } else if (d == DENSITY){
        for(int i =0; i <= GRID_HEIGHT; i++){
            for(int j= 0; j <= GRID_WIDTH; j++){
                get_density_color(data[i][j], MAX_DENSITY, &r, &g, &b);
                int pj = (i*(GRID_WIDTH+1) + j) * 3;
                img[pj] = (unsigned char)r;
                img[pj + 1] = (unsigned char)g;
                img[pj + 2] = (unsigned char)b;
            }
        }
    } else {
        printf("WIP\n");
    }

    fwrite(img, sizeof(unsigned char), (GRID_WIDTH+1) * (GRID_HEIGHT+1) * 3, f);
    free(img);
    fclose(f);
}

void save_frame_mv(int frame_num, data_type d, double data[GRID_HEIGHT+1][GRID_WIDTH+1], double rho[GRID_HEIGHT+1][GRID_WIDTH+1],double uj[GRID_HEIGHT+1][GRID_WIDTH+1],double ui[GRID_HEIGHT+1][GRID_WIDTH+1],int lower_limit) {
    char filename[64];
    sprintf(filename, "fluid_%04d.ppm", frame_num);
    FILE *f = fopen(filename, "w");
    if (!f) return;
    
    // Header PPM
    fprintf(f, "P6\n%d %d\n255\n", (GRID_WIDTH+1), (GRID_HEIGHT+1));

    unsigned char* img = (unsigned char*) calloc( (GRID_WIDTH+1) * (GRID_HEIGHT+1) * 3, sizeof(unsigned char));

    int r, g, b;

    for(int i =0; i <lower_limit; i++){
        for(int j= 0; j <= GRID_WIDTH; j++){
            int pj = (i*(GRID_WIDTH+1) + j) * 3;
            img[pj] = 0;
            img[pj + 1] = 0;
            img[pj + 2] = 0;  
        }
    }

    if (d == TEMPERATURE){
        for(int i =lower_limit; i <= GRID_HEIGHT; i++){
            for(int j= 0; j <= GRID_WIDTH; j++){
                double v_carre = uj[i][j]*uj[i][j] + ui[i][j]*ui[i][j];
                double T = (data[i][j] - 0.5*v_carre)/CV_T_VAPEUR;
                get_temperature_color((float)T, (float)MAX_DEGRE, &r, &g, &b);
                int pj = (i*(GRID_WIDTH+1) + j) * 3;
                img[pj] = (unsigned char)r;
                img[pj + 1] = (unsigned char)g;
                img[pj + 2] = (unsigned char)b;  
            }
        }
    } else if (d == DENSITY){
        for(int i =lower_limit; i <= GRID_HEIGHT; i++){
            for(int j= 0; j <= GRID_WIDTH; j++){
                get_density_color(data[i][j], MAX_DENSITY, &r, &g, &b);
                int pj = (i*(GRID_WIDTH+1) + j) * 3;
                img[pj] = (unsigned char)r;
                img[pj + 1] = (unsigned char)g;
                img[pj + 2] = (unsigned char)b;
            }
        }
    } else {
        printf("WIP\n");
    }

    fwrite(img, sizeof(unsigned char), (GRID_WIDTH+1) * (GRID_HEIGHT+1) * 3, f);
    free(img);
    fclose(f);
}