#include "writer.hpp"
#include "solver.hpp"



std::map<std::string, double> parameters(const std::string& nomfichier) {

	std::map<std::string, double> parameters;
	std::string key;
	std::ifstream file(nomfichier);
	double valeur;

	// petit erreur affichage même si lu correctement
	if (!file.is_open()) {
		std::cerr << "Erreur lors de l'ouverture du fichier : " << nomfichier << std::endl;
		return parameters;
	}
	
	// erreur dans la lecture du fichier
	while (file >> key >> valeur) {
		parameters[key] = valeur;
	}

	file.close();
	return parameters;
}

