#include "lbm.h"



/*-------------Tableau global utilitaires + structure---------------*/



int opp[Q] = {0,3,4,1,2,7,8,5,6}; //tableau des opposées pour les rebonds
int dot_x[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
int dot_y[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
double w[Q] = {4./9., 1./9., 1./9., 1./9., 1./9., 1./36., 1./36., 1./36., 1./36.};
double ax = 0.;
double ay = -9.81*(DT*DT/DELTA); //en unité du réseau

//Fonctions de répartition et grandeurs thermo

double f[GRID_HEIGHT+1][GRID_WIDTH+1][Q];
double f_new[GRID_HEIGHT+1][GRID_WIDTH+1][Q];
double h[GRID_HEIGHT+1][GRID_WIDTH+1][Q];
double h_new[GRID_HEIGHT+1][GRID_WIDTH+1][Q];
double rho[GRID_HEIGHT+1][GRID_WIDTH+1];
double ux[GRID_HEIGHT+1][GRID_WIDTH+1];
double uy[GRID_HEIGHT+1][GRID_WIDTH+1];
double E[GRID_HEIGHT+1][GRID_WIDTH+1];
double T[GRID_HEIGHT+1][GRID_WIDTH+1];

//Paramètre pour le cycle

int lower_limit = 16; //pour RA = 8
int lower_limit_new;
int higher_limit = 20;
int higher_limit_new;

////////////////////////// A - Simulation basique //////////////////////////

/*----------------- A-I - Initialisation de la simulation ------------------*/

void init_uniforme(){
    for (int i = 0; i <= GRID_HEIGHT; i++) {
        for (int j = 0; j <= GRID_WIDTH;j ++) {
            E[i][j] = 400000.0;
            T[i][j] = 400000.0 / CV_T_VAPEUR;
            rho[i][j] = 400.0;
            ux[i][j] = 0.0;
            uy[i][j] = 0.0;
            
            // Initialisation à l'équilibre (vitesse nulle)
            // f_eq = rho * w_i (simplifié car u=0 au début)
            for (int k = 0; k < Q; k++) {
                f[i][j][k]=feq(rho[i][j],ux[i][j],uy[i][j],k);
                h[i][j][k]=heq(rho[i][j],ux[i][j],uy[i][j],E[i][j],k);
            }
        }
    }
}

void init_bilateral(){
    for (int i = 0; i < GRID_HEIGHT/2; i++) {
        for (int j = 0; j <= GRID_WIDTH;j ++) {
            E[i][j] = 400000.0;
            rho[i][j] = 800.0;
            ux[i][j] = 0.0;
            uy[i][j] = 0.0;
            
            // Initialisation à l'équilibre (vitesse nulle)
            // f_eq = rho * w_i (simplifié car u=0 au début)
            for (int k = 0; k < Q; k++) {
                f[i][j][k]=feq(rho[i][j],ux[i][j],uy[i][j],k);
                h[i][j][k]=heq(rho[i][j],ux[i][j],uy[i][j],E[i][j],k);
            }
        }
    }
    for (int i = GRID_HEIGHT/2; i <= GRID_HEIGHT; i++) {
        for (int j = 0; j <= GRID_WIDTH;j ++) {
            E[i][j] = 800000.0;
            rho[i][j] = 800.0;
            ux[i][j] = 0.0;
            uy[i][j] = 0.0;
            
            // Initialisation à l'équilibre (vitesse nulle)
            // f_eq = rho * w_i (simplifié car u=0 au début)
            for (int k = 0; k < Q; k++) {
                f[i][j][k]=feq(rho[i][j],ux[i][j],uy[i][j],k);
                h[i][j][k]=heq(rho[i][j],ux[i][j],uy[i][j],E[i][j],k);
            }
        }
    }
}

/*------------------- A-II - Calcul des fonctions d'équilibres f et h ---------------------*/

double feq(double Rho,double U,double V,int k){ //U vitesse horizontale V vitesse verticale
    double cu = dot_x[k]*U + dot_y[k]*V;
    double v_carre = U*U + V*V;

    return w[k]*Rho*(1 + cu/(R*T0) + 0.5*cu*cu/(R*T0*R*T0) - 0.5*v_carre/(R*T0));
}

double heq(double Rho,double U,double V,double Etot, int k){ //U vitesse horizontale V vitesse verticale
    double cu = dot_x[k]*U + dot_y[k]*V;
    double v_carre = U*U + V*V;

    return w[k]*Rho*R*T0*(cu/(R*T0) + cu*cu/(R*T0*R*T0) - v_carre/(R*T0) + 0.5*((dot_x[k]*dot_x[k] + dot_y[k]*dot_y[k])/(R*T0) - D)) + Etot*feq(Rho,U,V,k);
}

/*------------------- A-III - Collision BGK-Operator ---------------------*/

void collision(){
    double FEQ;
    double HEQ;
    for(int i = 0;i<=GRID_HEIGHT;i++){
        for(int j=0;j<=GRID_WIDTH;j++){
            for (int k = 0; k < Q; k++) {
                FEQ=feq(rho[i][j],ux[i][j],uy[i][j],k);
                HEQ=heq(rho[i][j],ux[i][j],uy[i][j],E[i][j],k);
                double cu = dot_x[k]*ux[i][j] + dot_y[k]*uy[i][j];
                double ca = dot_x[k]*ax + dot_y[k]*ay;
                double au = ux[i][j]*ax + uy[i][j]*ay;
                double v_carre = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j]; 
                double z = cu - v_carre;
                double F = w[k]*rho[i][j]*(ca/(R*T0) + (ca*cu)/(R*T0*R*T0) - au/(R*T0));
                double phi = f[i][j][k] - FEQ  + F*DT/2;
                double q = (w[k]*rho[i][j]*E[i][j]/(R*T0) + FEQ + (1 - WF/2)*phi)*ca;
                f_new[i][j][k]=f[i][j][k]-WF*(f[i][j][k] - FEQ) + ((2-WF)/2)*DT*F; 
                h_new[i][j][k]=h[i][j][k]-WH*(h[i][j][k] - HEQ) + ((2-WF)/2)*DT*q + (WH - WF)*z*phi;
            }
        }
    }
}

/*------------------- A-IV - Diffusion ---------------------*/

void diffusion(){
    int id;
    int jd;
    for(int i = 0;i<=GRID_HEIGHT;i++){
        for(int j=0;j<=GRID_WIDTH;j++){
            for (int k = 0; k < Q; k++) {
                id = i - dot_y[k];
                jd = j - dot_x[k];

                if(id>=0 && id<=GRID_HEIGHT && jd>=0 && jd<=GRID_WIDTH){
                    f[i][j][k] = f_new[id][jd][k];
                    h[i][j][k] = h_new[id][jd][k];
                }
                
            }
        }
    }
}

/*------------------- A-V - Gestion des conditions aux limites ---------------------*/

void rebond(){

    for(int j=0;j<=GRID_WIDTH;j++){
        
        //Rebond sur le haut de la grille

        f[0][j][4] = f_new[0][j][2];
        f[0][j][7] = f_new[0][j][5];
        f[0][j][8] = f_new[0][j][6];

        h[0][j][4] = h_new[0][j][2];
        h[0][j][7] = h_new[0][j][5];
        h[0][j][8] = h_new[0][j][6];

        //Rebond sur le haut de la grille

        f[GRID_HEIGHT][j][2] = f_new[GRID_HEIGHT][j][4];
        f[GRID_HEIGHT][j][5] = f_new[GRID_HEIGHT][j][7];
        f[GRID_HEIGHT][j][6] = f_new[GRID_HEIGHT][j][8];

        h[GRID_HEIGHT][j][2] = h_new[GRID_HEIGHT][j][4];
        h[GRID_HEIGHT][j][5] = h_new[GRID_HEIGHT][j][7];
        h[GRID_HEIGHT][j][6] = h_new[GRID_HEIGHT][j][8];

    }
    
    for(int i=0;i<=GRID_HEIGHT;i++){

        //Rebond sur le côté gauche

        f[i][0][1] = f_new[i][0][3];
        f[i][0][8] = f_new[i][0][6];
        f[i][0][5] = f_new[i][0][7];

        h[i][0][1] = h_new[i][0][3];
        h[i][0][8] = h_new[i][0][6];
        h[i][0][5] = h_new[i][0][7];
        
        //Rebond sur le côté droit

        f[i][GRID_WIDTH][3] = f_new[i][GRID_WIDTH][1];
        f[i][GRID_WIDTH][6] = f_new[i][GRID_WIDTH][8];
        f[i][GRID_WIDTH][7] = f_new[i][GRID_WIDTH][5];

        h[i][GRID_WIDTH][3] = h_new[i][GRID_WIDTH][1];
        h[i][GRID_WIDTH][6] = h_new[i][GRID_WIDTH][8];
        h[i][GRID_WIDTH][7] = h_new[i][GRID_WIDTH][5];
    }

}

/*------------------- A-VI - Calcul des grandeurs thermodynamiques ---------------------*/

void densite_velocite(){
    for(int i=0;i<=GRID_HEIGHT;i++){
        for(int j=0;j<=GRID_WIDTH;j++){

            //calcul de la densité
            rho[i][j] = 0;
            for(int k=0;k<Q;k++){
                rho[i][j] += f[i][j][k];
            }

            //calcul de la vélocité

            ux[i][j] = (f[i][j][1] + f[i][j][5] + f[i][j][8] - f[i][j][3] - f[i][j][6] - f[i][j][7] + 0.5*DT*rho[i][j]*ax) / rho[i][j];
            uy[i][j] = (f[i][j][5] + f[i][j][6] + f[i][j][2] - f[i][j][7] - f[i][j][8] - f[i][j][4] + 0.5*DT*rho[i][j]*ay) / rho[i][j];

            //calcul de l'énergie totale

            E[i][j] = 0;
            for(int k=0;k<Q;k++){
                E[i][j] += h[i][j][k] / rho[i][j];
            }

            double au = ux[i][j]*ax + uy[i][j]*ay;
            E[i][j] += 0.5*DT*au;

            // calcul de la température (sucre)
            
            double v_carre = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
            T[i][j] = (E[i][j] - 0.5*v_carre)/CV_T_VAPEUR;

        }
    }
}

/*------------------- A-VIII - Fonction principale ---------------------*/

void lbm(int nb_frames,void (*init)()){
    init();
    for(int i=0;i<nb_frames;i++){
        collision();
        diffusion();
        rebond();
        densite_velocite();
        save_frame(i,TEMPERATURE,E,rho,ux,uy);
        macro_physique(rho,E,ux,uy);
    }
}






/*##########################################################################################################################*/
/*###################################  TEST EN VUE D'IMPLEMENTER UN CYCLE DE STIRLING  #####################################*/
/*##########################################################################################################################*/






////////////////////////// B - Simulation avec un transfert thermique //////////////////////////



void rebond_test_tt(){ //On modifie le rebond pour introduire un chauffage par le côté droit de la grille (j=GRID_WIDTH)
    for(int j=0;j<=GRID_WIDTH;j++){
        
        //Rebond sur le haut de la grille

        f[0][j][4] = f_new[0][j][2];
        f[0][j][7] = f_new[0][j][5];
        f[0][j][8] = f_new[0][j][6];

        h[0][j][4] = h_new[0][j][2];
        h[0][j][7] = h_new[0][j][5];
        h[0][j][8] = h_new[0][j][6];

        //Rebond sur le haut de la grille

        f[GRID_HEIGHT][j][2]=f_new[GRID_HEIGHT][j][4];
        f[GRID_HEIGHT][j][5]=f_new[GRID_HEIGHT][j][7];
        f[GRID_HEIGHT][j][6]=f_new[GRID_HEIGHT][j][8];

        h[GRID_HEIGHT][j][2]=h_new[GRID_HEIGHT][j][4];
        h[GRID_HEIGHT][j][5]=h_new[GRID_HEIGHT][j][7];
        h[GRID_HEIGHT][j][6]=h_new[GRID_HEIGHT][j][8];

    }
    
    for(int i=0;i<=GRID_HEIGHT;i++){

        //Rebond sur le côté gauche

        f[i][0][1] = f_new[i][0][3];
        f[i][0][8] = f_new[i][0][6];
        f[i][0][5] = f_new[i][0][7];

        h[i][0][1] = h_new[i][0][3] + w[1]*H*(T_chaud - T[i][0]);
        h[i][0][8] = h_new[i][0][6] + w[8]*H*(T_chaud - T[i][0]);
        h[i][0][5] = h_new[i][0][7] + w[5]*H*(T_chaud - T[i][0]);
        
        //Rebond sur le côté droit

        f[i][GRID_WIDTH][3] = f_new[i][GRID_WIDTH][1];
        f[i][GRID_WIDTH][6] = f_new[i][GRID_WIDTH][8];
        f[i][GRID_WIDTH][7] = f_new[i][GRID_WIDTH][5];

        h[i][GRID_WIDTH][3] = h_new[i][GRID_WIDTH][1] + w[1]*H*(T_chaud - T[i][GRID_WIDTH]); //ajout d'un transfert thermique sur le côté droit
        h[i][GRID_WIDTH][6] = h_new[i][GRID_WIDTH][8] + w[8]*H*(T_chaud - T[i][GRID_WIDTH]);
        h[i][GRID_WIDTH][7] = h_new[i][GRID_WIDTH][5] + w[5]*H*(T_chaud - T[i][GRID_WIDTH]);
    }
}

void lbm_test_tt(int nb_frames,void (*init)()){
    init();
    for(int i=0;i<nb_frames;i++){
        collision();
        diffusion();
        rebond_test_tt();
        densite_velocite();
        save_frame(i,TEMPERATURE,E,rho,ux,uy);
        macro_physique(rho,E,ux,uy);
    }
}



////////////////////////// C - Simulation avec des mouvements de murs //////////////////////////



/*------------------- C-III - Collision BGK-Operator ---------------------*/

void collision_mv(){
    double FEQ;
    double HEQ;
    for(int i = lower_limit;i<=GRID_HEIGHT;i++){
        for(int j = 0;j<=GRID_WIDTH;j++){
            for (int k = 0; k < Q; k++) {
                FEQ=feq(rho[i][j],ux[i][j],uy[i][j],k);
                HEQ=heq(rho[i][j],ux[i][j],uy[i][j],E[i][j],k);

                double cu = dot_x[k]*ux[i][j] + dot_y[k]*uy[i][j];
                double ca = dot_x[k]*ax + dot_y[k]*ay;
                double au = ux[i][j]*ax + uy[i][j]*ay;
                double v_carre = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j]; 

                double z = cu - v_carre;
                
                double F = w[k]*rho[i][j]*(ca/(R*T0) + (ca*cu)/(R*T0*R*T0) - au/(R*T0));
                double phi = f[i][j][k] - FEQ  + F*DT/2;
                double q = (w[k]*rho[i][j]*E[i][j]/(R*T0) + FEQ + (1 - WF/2)*phi)*ca;
                

                f_new[i][j][k]=f[i][j][k]-WF*(f[i][j][k] - FEQ) + ((2-WF)/2)*DT*F; 
                h_new[i][j][k]=h[i][j][k]-WH*(h[i][j][k] - HEQ) + ((2-WF)/2)*DT*q + (WH - WF)*z*phi;
            }
        }
    }
}


/*------------------- C-IV - Diffusion ---------------------*/


void diffusion_mv(int i){
    int id;
    int jd;
    lower_limit_new = ex(i);
    if(lower_limit == lower_limit_new){
        for(int i = lower_limit;i<=GRID_HEIGHT;i++){
            for(int j=0;j<=GRID_WIDTH;j++){
                for (int k = 0; k < Q; k++){
                    id = i - dot_y[k];
                    jd = j - dot_x[k];

                    if(id>=lower_limit && id<=GRID_HEIGHT && jd>=0 && jd<=GRID_WIDTH){
                        f[i][j][k] = f_new[id][jd][k];
                        h[i][j][k] = h_new[id][jd][k];
                    }

                }
            }
        }
    }else if( lower_limit_new == lower_limit + 1){
        for(int i = lower_limit + 2;i<=GRID_HEIGHT;i++){
            for(int j=0;j<=GRID_WIDTH;j++){
                for (int k = 0; k < Q; k++){
                    id = i - dot_y[k];
                    jd = j - dot_x[k];

                    if(id>=lower_limit + 2 && id<=GRID_HEIGHT && jd>=0 && jd<=GRID_WIDTH){
                        f[i][j][k] = f_new[id][jd][k];
                        h[i][j][k] = h_new[id][jd][k];
                    }
                }
            }
        }
        for(int j=0;j<=GRID_WIDTH;j++){// on "compacte" sur une rangée, celle déjà présente et celle qui se fait pousser → rebond sur la rangé précédente
            for (int k = 0; k < Q; k++){
                f[lower_limit + 1][j][k] = f_new[lower_limit][j][k] + f_new[lower_limit + 1][j][k];
                h[lower_limit + 1][j][k] = h_new[lower_limit][j][k] + f_new[lower_limit + 1][j][k];
            }
        }


    }else if( lower_limit_new == lower_limit - 1){
        for(int i = lower_limit;i<=GRID_HEIGHT;i++){
            for(int j=0;j<=GRID_WIDTH;j++){
                for (int k = 0; k < Q; k++){
                    id = i - dot_y[k];
                    jd = j - dot_x[k];

                    if(id>=lower_limit - 1 && id<=GRID_HEIGHT && jd>=0 && jd<=GRID_WIDTH){ //gestion des conditions aux limites dans rebond()
                        f[i][j][k] = f_new[id][jd][k];
                        h[i][j][k] = h_new[id][jd][k];
                    }
                    
                }
            }
        }

    }else{
        printf("Lower_limit : %d  |  Lower_limit_new : %d\n",lower_limit,lower_limit_new);
        assert(false);
    }
}


/*------------------- C-V - Gestion des conditions aux limites ---------------------*/


void rebond_mv(){
    if(lower_limit == lower_limit_new){
        for(int j=0;j<=GRID_WIDTH;j++){
            
            //Rebond sur le haut de la grille

            f[lower_limit][j][4] = f_new[lower_limit][j][2];
            f[lower_limit][j][7] = f_new[lower_limit][j][5];
            f[lower_limit][j][8] = f_new[lower_limit][j][6];

            h[lower_limit][j][4] = h_new[lower_limit][j][2];
            h[lower_limit][j][7] = h_new[lower_limit][j][5];
            h[lower_limit][j][8] = h_new[lower_limit][j][6];

            //Rebond sur le haut de la grille

            f[GRID_HEIGHT][j][2] = f_new[GRID_HEIGHT][j][4];
            f[GRID_HEIGHT][j][5] = f_new[GRID_HEIGHT][j][7];
            f[GRID_HEIGHT][j][6] = f_new[GRID_HEIGHT][j][8];

            h[GRID_HEIGHT][j][2] = h_new[GRID_HEIGHT][j][4];
            h[GRID_HEIGHT][j][5] = h_new[GRID_HEIGHT][j][7];
            h[GRID_HEIGHT][j][6] = h_new[GRID_HEIGHT][j][8];

        }
        
        for(int i=lower_limit;i<=GRID_HEIGHT;i++){

            //Rebond sur le côté gauche

            f[i][0][1] = f_new[i][0][3];
            f[i][0][8] = f_new[i][0][6];
            f[i][0][5] = f_new[i][0][7];

            h[i][0][1] = h_new[i][0][3];
            h[i][0][8] = h_new[i][0][6];
            h[i][0][5] = h_new[i][0][7];
            
            //Rebond sur le côté droit

            f[i][GRID_WIDTH][3] = f_new[i][GRID_WIDTH][1];
            f[i][GRID_WIDTH][6] = f_new[i][GRID_WIDTH][8];
            f[i][GRID_WIDTH][7] = f_new[i][GRID_WIDTH][5];

            h[i][GRID_WIDTH][3] = h_new[i][GRID_WIDTH][1];
            h[i][GRID_WIDTH][6] = h_new[i][GRID_WIDTH][8];
            h[i][GRID_WIDTH][7] = h_new[i][GRID_WIDTH][5];
        }
    }else if(lower_limit_new == lower_limit + 1){ //rebond une case plus bas encore
        for(int j=0;j<=GRID_WIDTH;j++){ 
        //Rebond sur le haut de la grille

            f[lower_limit+2][j][4] = f_new[lower_limit+2][j][2];
            f[lower_limit+2][j][7] = f_new[lower_limit+2][j][5];
            f[lower_limit+2][j][8] = f_new[lower_limit+2][j][6];

            h[lower_limit+2][j][4] = h_new[lower_limit+2][j][2];
            h[lower_limit+2][j][7] = h_new[lower_limit+2][j][5];
            h[lower_limit+2][j][8] = h_new[lower_limit+2][j][6];

            //Rebond sur le haut de la grille

            f[GRID_HEIGHT][j][2] = f_new[GRID_HEIGHT][j][4];
            f[GRID_HEIGHT][j][5] = f_new[GRID_HEIGHT][j][7];
            f[GRID_HEIGHT][j][6] = f_new[GRID_HEIGHT][j][8];

            h[GRID_HEIGHT][j][2] = h_new[GRID_HEIGHT][j][4];
            h[GRID_HEIGHT][j][5] = h_new[GRID_HEIGHT][j][7];
            h[GRID_HEIGHT][j][6] = h_new[GRID_HEIGHT][j][8];

        }
        
        for(int i=lower_limit + 2;i<=GRID_HEIGHT;i++){

            //Rebond sur le côté gauche

            f[i][0][1] = f_new[i][0][3];
            f[i][0][8] = f_new[i][0][6];
            f[i][0][5] = f_new[i][0][7];

            h[i][0][1] = h_new[i][0][3];
            h[i][0][8] = h_new[i][0][6];
            h[i][0][5] = h_new[i][0][7];
            
            //Rebond sur le côté droit

            f[i][GRID_WIDTH][3] = f_new[i][GRID_WIDTH][1];
            f[i][GRID_WIDTH][6] = f_new[i][GRID_WIDTH][8];
            f[i][GRID_WIDTH][7] = f_new[i][GRID_WIDTH][5];

            h[i][GRID_WIDTH][3] = h_new[i][GRID_WIDTH][1];
            h[i][GRID_WIDTH][6] = h_new[i][GRID_WIDTH][8];
            h[i][GRID_WIDTH][7] = h_new[i][GRID_WIDTH][5];
        }
    }else if(lower_limit_new == lower_limit - 1){
        
        for(int i=lower_limit;i<=GRID_HEIGHT;i++){ // gestion d'abord des rebonds sur les côtés pour ne perdre la conservatino de l'énergie (very important)

            //Rebond sur le côté gauche

            f[i][0][1] = f_new[i][0][3];
            f[i][0][8] = f_new[i][0][6];
            f[i][0][5] = f_new[i][0][7];

            h[i][0][1] = h_new[i][0][3];
            h[i][0][8] = h_new[i][0][6];
            h[i][0][5] = h_new[i][0][7];
            
            //Rebond sur le côté droit

            f[i][GRID_WIDTH][3] = f_new[i][GRID_WIDTH][1];
            f[i][GRID_WIDTH][6] = f_new[i][GRID_WIDTH][8];
            f[i][GRID_WIDTH][7] = f_new[i][GRID_WIDTH][5];

            h[i][GRID_WIDTH][3] = h_new[i][GRID_WIDTH][1];
            h[i][GRID_WIDTH][6] = h_new[i][GRID_WIDTH][8];
            h[i][GRID_WIDTH][7] = h_new[i][GRID_WIDTH][5];
        }
        //on remplit les trous avec du "vide"(0)
        
        //gestions des coins
        f[lower_limit][0][4] = 0; //coin gauche i=lower_limit
        f[lower_limit][0][7] = 0;

        h[lower_limit][0][4] = 0;
        h[lower_limit][0][7] = 0;

        f[lower_limit][GRID_WIDTH][4] = 0; //coin droit i=lower_limit
        f[lower_limit][GRID_WIDTH][8] = 0;

        h[lower_limit][GRID_WIDTH][4] = 0;
        h[lower_limit][GRID_WIDTH][8] = 0;

        f[lower_limit - 1][0][5] = 0;  //coin gauche i=lower_limit - 1

        h[lower_limit - 1][0][5] = 0;

        f[lower_limit - 1][GRID_WIDTH][6] = 0;  //coin droit i=lower_limit - 1

        h[lower_limit - 1][GRID_WIDTH][6] = 0;

        //boucles

        for(int j=1;j<GRID_WIDTH;j++){//boucle i=lower_limit (les coins non inclus dans la boucle)
            f[lower_limit][j][4] = 0;
            f[lower_limit][j][7] = 0;
            f[lower_limit][j][8] = 0;

            h[lower_limit][j][4] = 0;
            h[lower_limit][j][7] = 0;
            h[lower_limit][j][8] = 0;
        }

        for(int j=0;j<=GRID_WIDTH;j++){//boucle i=lower_limit-1 (les coins inclus dans la boucle)
            f[lower_limit - 1][j][0] = 0;
            f[lower_limit - 1][j][1] = 0;
            f[lower_limit - 1][j][3] = 0;
            f[lower_limit - 1][j][4] = 0;
            f[lower_limit - 1][j][7] = 0;
            f[lower_limit - 1][j][8] = 0;

            h[lower_limit - 1][j][0] = 0;
            h[lower_limit - 1][j][1] = 0;
            h[lower_limit - 1][j][3] = 0;
            h[lower_limit - 1][j][4] = 0;
            h[lower_limit - 1][j][7] = 0;
            h[lower_limit - 1][j][8] = 0;
        }
    }
    
    lower_limit = lower_limit_new;
}

/*------------------- C-VI - Calcul des grandeurs thermodynamiques ---------------------*/

void densite_velocite_mv(){
    for(int i=lower_limit;i<=GRID_HEIGHT;i++){
        for(int j=0;j<=GRID_WIDTH;j++){

            //calcul de la densité
            rho[i][j] = 0;
            for(int k=0;k<Q;k++){
                rho[i][j] += f[i][j][k];
            }

            //calcul de la vélocité

            ux[i][j] = (f[i][j][1] + f[i][j][5] + f[i][j][8] - f[i][j][3] - f[i][j][6] - f[i][j][7] + 0.5*DT*rho[i][j]*ax) / rho[i][j];
            uy[i][j] = (f[i][j][5] + f[i][j][6] + f[i][j][2] - f[i][j][7] - f[i][j][8] - f[i][j][4] + 0.5*DT*rho[i][j]*ay) / rho[i][j];

            //calcul de l'énergie totale

            E[i][j] = 0;
            for(int k=0;k<Q;k++){
                E[i][j] += h[i][j][k] / rho[i][j];
            }

            double au = ux[i][j]*ax + uy[i][j]*ay;
            E[i][j] += 0.5*DT*au;

            // calcul de la température (sucre)
            
            double v_carre = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
            T[i][j] = (E[i][j] - 0.5*v_carre)/CV_T_VAPEUR;

        }
    }
}

/*------------------- C-VIII - Fonction principale ---------------------*/

void lbm_test_mv(int nb_frames,void (*init)()){
    init();
    for(int i=0;i<nb_frames;i++){
        collision_mv();
        diffusion_mv(i);
        rebond_mv();
        densite_velocite_mv();
        save_frame_mv(i,DENSITY,rho,rho,ux,uy,lower_limit);
        macro_mv_data(rho,E,ux,uy,lower_limit,i);
    }
}
