#include "writer.hpp"


std::map<std::string, double> params(const std::string& nomfichier) {

	std::map<std::string, double> parameters;
	std::string key;
	std::ifstream file(nomfichier);
	double valeur;

	if (!file.is_open()) {
		std::cerr << "Erreur d'ouverture du fichier: " << nomfichier << std::endl;
		return parameters;
	}

	while (file >> key >> valeur) {
		if (key[0] == '#') {
			std::string ignore;
			std::getline(file, ignore);
			continue;
		}

		parameters[key] = valeur;
	}

	file.close();
	return parameters;
}

