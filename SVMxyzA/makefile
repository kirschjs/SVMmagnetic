.SUFFIXES: .o .cpp .h
CFILES = CFILES
S = $(wildcard $(CFILES)/*.cpp)
SO = $(subst $(CFILES),OFILES,$(subst cpp,o,$(S)))

COMPILER =  c++
FLAGSopt =  -Ofast
#FLAGSdbg =  -Wall -fbounds-check
FLAGS = $(FLAGSopt)
#FLAGS = $(FLAGSdbg)

LIBS=./eigen/Eigen

svm: $(SO) 
	$(COMPILER) $(FLAGS) -o svm.x $^ 

all: $(SO) 
	$(COMPILER) $(FLAGS) -o svm.x -I$(LIBS) $^ 

OFILES/%.o: $(CFILES)/%.cpp
	$(COMPILER) -o $@ -c $(FLAGS)   $< -I$(LIBS)

clean::
	-rm -f *.o *.mod *~ *.x *# OFILES/*.o

