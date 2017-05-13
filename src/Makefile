all: ML.so MLe.so MLeq.so MLWPUE.so MLNPUE.so MSM1S.so MSM23.so 

ML.so: ML.cpp
	Rscript --vanilla -e "TMB::compile('ML.cpp')"

MLe.so: MLe.cpp
	Rscript --vanilla -e "TMB::compile('MLe.cpp')"

MLeq.so: MLeq.cpp
	Rscript --vanilla -e "TMB::compile('MLeq.cpp')"

MLWPUE.so: MLWPUE.cpp
	Rscript --vanilla -e "TMB::compile('MLWPUE.cpp')"
	
MLNPUE.so: MLNPUE.cpp
	Rscript --vanilla -e "TMB::compile('MLNPUE.cpp')"
	
MSM1S.so: MSM1S.cpp
	Rscript --vanilla -e "TMB::compile('MSM1S.cpp')"
	
MSM23.so: MSM23.cpp
	Rscript --vanilla -e "TMB::compile('MSM23.cpp')"

clean:
	rm -rf *o
	