#include "lbm.h"

/*---------------------- MAIN -------------------------*/

int main(int argc,char** argv){ // 1 argument : temps de la simulation (nb de frame)
    assert(argc == 2);
    int nb_frames = atoi(argv[1]);
    lbm_test_mv(nb_frames,&init_uniforme);
    return EXIT_SUCCESS;
}