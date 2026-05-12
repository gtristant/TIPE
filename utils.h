#ifndef UTLIS_H
#define UTILS_H

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

/*------------------ Grandeurs et paramètres utiles à la simulation de fluides ------------------*/

//Grandeurs caractéristiques du réseau

#define DELTA 0.01 //en m
#define DT 0.1 //s

//Paramètres de la grille

#define GRID_WIDTH 49
#define GRID_HEIGHT 49
#define Q 9
#define D 2

//Grandeurs pour le fluide

#define CV_T_VAPEUR 1500 //en J/kg/K
#define CV_T_EAU 4200 // en J/kg/K (approximé, varie légèrement selon la pression)

/*------------------ Grandeurs et paramètres utiles à la simulation du  cycle thermo ------------------*/

#define þ 0.1 //en tour/minute
#define RA 8 //en unités de la sim
#define L 3*RA //en unités de la sim


void get_temperature_color(float temp, float max_degre, int *r, int *g, int *b);
void get_density_color(double density, double max_density, int* r, int* g, int* b);
void print_data(double data[GRID_HEIGHT+1][GRID_WIDTH + 1]);
void macro_physique(double rho[GRID_HEIGHT+1][GRID_WIDTH+1],double E[GRID_HEIGHT+1][GRID_WIDTH+1],double ux[GRID_HEIGHT+1][GRID_WIDTH+1],double uy[GRID_HEIGHT+1][GRID_WIDTH+1]);
void macro_physique_mv(double rho[GRID_HEIGHT+1][GRID_WIDTH+1],double E[GRID_HEIGHT+1][GRID_WIDTH+1],double ux[GRID_HEIGHT+1][GRID_WIDTH+1],double uy[GRID_HEIGHT+1][GRID_WIDTH+1],int lower_limit);
void macro_mv_data(double rho[GRID_HEIGHT+1][GRID_WIDTH+1],double E[GRID_HEIGHT+1][GRID_WIDTH+1],double ux[GRID_HEIGHT+1][GRID_WIDTH+1],double uy[GRID_HEIGHT+1][GRID_WIDTH+1],int lower_limit,int i);


/*------------------- Fonctions utilitaires pour la simulation du cycle -------------------*/

int ex(int i);
int ey(int i);


#endif //UTILS_H