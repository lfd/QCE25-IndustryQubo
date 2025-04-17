.PHONY: all

all: extract plots pre_compile

extract:
	7z x results/Annealing1.7z -oresults/
	7z x results/Annealing2.7z -oresults/
	7z x results/LRQAOA.7z -oresults/
	7z x results/Simulated.7z -oresults/

pre_compile:
	./pre_compile_pics.sh

plots: 
	Rscript PaperPlots.r


