#include "writer.hpp"
#include "solver.hpp"

#include <sstream>



std::map<std::string, double> parameters(const std::string& nomfichier) {

	std::map<std::string, double> parameters;
	std::ifstream file(nomfichier);
	std::string line;

	// petit erreur affichage même si lu correctement
	if (!file.is_open()) {
		std::cerr << "Erreur lors de l'ouverture du fichier : " << nomfichier << std::endl;
		return parameters;
	}
	
	while (std::getline(file, line)) {
		if (line.empty()) {
			continue;
		}

		const std::size_t comment_pos = line.find('#');
		if (comment_pos != std::string::npos) {
			line = line.substr(0, comment_pos);
		}

		std::istringstream iss(line);
		std::string key;
		double valeur;
		if (iss >> key >> valeur) {
			parameters[key] = valeur;
		}
	}

	file.close();
	return parameters;
}

