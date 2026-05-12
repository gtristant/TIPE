#include "utils.h"

void get_temperature_color(float temp, float max_degre, int *r, int *g, int *b) {
    // On sature la température entre 0 et max_degre pour éviter les erreurs
    if (temp < 0) temp = 0;
    if (temp > max_degre) temp = max_degre;

    float ratio = temp / max_degre; // Valeur entre 0.0 et 1.0

    if (ratio < 0.33) {
        // Segment 1 : Bleu Sombre -> Rouge
        // On passe de (0, 0, 100) à (255, 0, 0)
        float f = ratio / 0.33f;
        *r = (int)(f * 255);
        *g = 0;
        *b = (int)((1.0f - f) * 100);
    } 
    else if (ratio < 0.66) {
        // Segment 2 : Rouge -> Jaune
        // On passe de (255, 0, 0) à (255, 255, 0)
        float f = (ratio - 0.33f) / 0.33f;
        *r = 255;
        *g = (int)(f * 255);
        *b = 0;
    } 
    else {
        // Segment 3 : Jaune -> Blanc
        // On passe de (255, 255, 0) à (255, 255, 255)
        float f = (ratio - 0.66f) / 0.34f;
        *r = 255;
        *g = 255;
        *b = (int)(f * 255);
    }
}

void get_density_color(double density, double max_density, int* r, int* g, int* b) {
    // On normalise la densité entre 0 et 255
    int val = (int)((density / max_density) * 255);
    if (val > 255) val = 255;
    if (val < 0) val = 0;

    *r = val; // Plus c'est dense, plus c'est blanc
    *g = val;
    *b = val;
}

void print_data(double data[GRID_HEIGHT+1][GRID_WIDTH + 1]){
    printf("\n\n");
    for(int i=0;i<GRID_HEIGHT;i++){
        for(int j=0;j<GRID_WIDTH;j++){
            printf("%d ",(int)floor(data[i][j]));
        }
        printf("\n");
    }
    printf("\n\n");
}

void macro_physique(double rho[GRID_HEIGHT+1][GRID_WIDTH+1],double E[GRID_HEIGHT+1][GRID_WIDTH+1],double ux[GRID_HEIGHT+1][GRID_WIDTH+1],double uy[GRID_HEIGHT+1][GRID_WIDTH+1]){
    double global_rho = 0;
    double global_T = 0;
    int v2;
    for(int i=0;i<=GRID_HEIGHT;i++){
        for(int j=0;j<=GRID_WIDTH;j++){
            global_rho += rho[i][j];
            v2 = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
            global_T += ((E[i][j] - v2/2)/CV_T_VAPEUR)*rho[i][j];
        }
    }
    global_T /= global_rho;
    global_rho /= (GRID_HEIGHT*GRID_WIDTH);
    
    printf("Température globale : %lf | Densité globale: %lf\n",global_T,global_rho);
}

void macro_physique_mv(double rho[GRID_HEIGHT+1][GRID_WIDTH+1],double E[GRID_HEIGHT+1][GRID_WIDTH+1],double ux[GRID_HEIGHT+1][GRID_WIDTH+1],double uy[GRID_HEIGHT+1][GRID_WIDTH+1],int lower_limit){
    double global_rho = 0;
    double global_T = 0;
    int v2;
    for(int i=lower_limit;i<=GRID_HEIGHT;i++){
        for(int j=0;j<=GRID_WIDTH;j++){
            global_rho += rho[i][j];
            v2 = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
            global_T += ((E[i][j] - v2/2)/CV_T_VAPEUR)*rho[i][j];
        }
    }
    global_T /= global_rho;
    global_rho /= (GRID_HEIGHT*GRID_WIDTH);
    
    printf("Température globale : %lf | Densité globale: %lf\n",global_T,global_rho);
}

int ex(int i){
    double theta = þ * (i*DT/60) * 2 * 3.141592;
    double ex = sqrt(L*L - RA*RA*sin(theta)*sin(theta)) - RA*(2 - cos(theta));
    return (int)ex;
}

int ey(int i){
    double theta = þ * (i*DT/60) * 2 * 3.141592;
    double ey = sqrt(L*L - RA*RA*cos(theta)*cos(theta)) - RA*(2 + sin(theta));
    return (int)ey;
}

void macro_mv_data(double rho[GRID_HEIGHT+1][GRID_WIDTH+1],double E[GRID_HEIGHT+1][GRID_WIDTH+1],double ux[GRID_HEIGHT+1][GRID_WIDTH+1],double uy[GRID_HEIGHT+1][GRID_WIDTH+1],int lower_limit,int i){
    if(i%100==0){
        
        double global_rho = 0;
        double global_T = 0;
        int v2;
        for(int i=lower_limit;i<=GRID_HEIGHT;i++){
            for(int j=0;j<=GRID_WIDTH;j++){
                global_rho += rho[i][j];
                v2 = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
                global_T += ((E[i][j] - v2/2)/CV_T_VAPEUR)*rho[i][j];
            }
        }
        global_T /= global_rho;
        global_rho /= (GRID_HEIGHT*GRID_WIDTH);

        printf("Frame %d : Température globale : %lf | Densité globale: %lf\n",i,global_T,global_rho);
    }

    
}
