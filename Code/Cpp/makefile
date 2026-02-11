CXX = g++ 	
SRC = src/main.cpp src/solver.cpp src/writer.cpp src/formule.cpp src/teste.cpp # fichiers sources
CXXFLAGS = -Wall -DDEBUG -std=c++17 # options de compilation
OBJDIR = build
# fichiers objets
OBJS = $(OBJDIR)/main.o $(OBJDIR)/solver.o $(OBJDIR)/writer.o $(OBJDIR)/formule.o $(OBJDIR)/teste.o
OPTIM_FLAGS = -O3
RES = resultats

#Targets
Target = main

all: $(Target) #cible par défaut

$(OBJDIR): #crée le dossier build s'il n'existe pas
	mkdir -p $(OBJDIR)
	
$(RES):
	mkdir -p $(RES)

$(Target): $(OBJS) 	#cible principale
	$(CXX) $(CXXFLAGS) -o $(Target) $(OBJS)

$(OBJDIR)/%.o: src/%.cpp $(OBJDIR) # règle de compilation des fichiers .cpp en .o
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean: #cible nettoyage
	rm -rf $(OBJDIR) $(Target) $(RES)

