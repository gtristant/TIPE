#ifndef LBM_H
#define LBM_H

#include "videomaker.h"

/*---------------------Définition des constantes utiles à la modélisation (en unités SI) ----------------------*/



#define TAU_F 0.73 //s (vapeur d'eau)
#define TAU_H 1.0 //s

#define WF (2*DT)/(2*TAU_F + DT) //s
#define WH (2*DT)/(2*TAU_H + DT) //s

#define R 8.314 // adimensionnée (J/mol/K)
#define T0 500 // adimensionnée

#define CS 1/1.7320508 //sqrt(DELTA*DELTA/(3*DT*DT)) m/s (vitesse du son du réseau)
#define T_chaud 1000 //K
#define T_froid 100 //K
#define H 1000000 //Coefficient de convection thermique


/*--------------------- DEFINITIONS DES FONCTIONS -----------------------*/

void init_uniforme();
void init_bilateral();
double feq(double Rho,double U,double V,int k);
double heq(double Rho,double U,double V,double Etot, int k);
void collision();
void diffusion();
void rebond();
void densite_velocite();
void lbm(int nb_frames,void (*init)());

//fonction de test

void lbm_test_tt(int nb_frames,void (*init)());
void lbm_test_mv(int nb_frames,void (*init)());

#endif //LBM_H