all: clean comp run build

comp: main.c lbm.c videomaker.c utils.c
	gcc main.c lbm.c videomaker.c utils.c -o lattice -lm

run:
	./lattice $(NB_FRAME)

clean :
	rm -f lattice *.o *.ppm

trash:
	rm -f lattice *.o *.ppm *.mp4

build_ :
	ffmpeg -i fluid_%04d.ppm -c:v libx264 -vf "fps=50,format=yuv420p" sim.mp4

build: build_ clean

debug:
	gcc main.c lbm.c videomaker.c utils.c -o lattice -Wall -lm

save:
	mv *.mp4 saved_sim

name:
	mv sim.mp4 $(NAME)
