#include"solver.hpp"
#include"writer.hpp"
#include"formule.hpp"



int main(int argc, char** argv){

    if (argc != 6) {
        std::cerr << "Usage: NADH_m ; ADP_m; Dpsi ; Ca_m ; tfinal" << std::endl;
        return 1;
    }

    // création de l'état initials a partire des données
    etat x;
    x.NADH_m = std::stod(argv[1]);
    x.ADP_m = std::stod(argv[2]);
    x.deltaPsi = std::stod(argv[3]);
    x.Ca_m = std::stod(argv[4]);


    // paramètres du modèle
    std::map<std::string,double>params = parameters("parameters.txt");

    std::cout << params["p4"]<<std::endl;  //test d'affichage
    //on génère un fichier d'ecriture de sortie

    std::ofstream sortie("resultats/resultat.txt");
    sortie << "t" << "\t"<< "NADH_m" << "\t" << "ADP_m" << "\t" << "deltaPsi" << "\t" << "Ca_m" << "\t" << "Ca_c" << "\t" << "FBP" << std::endl;
    double t=0.0;
    double dt=0.001; // pas de temps en milisecondes
    double tfinal = std::stod(argv[5]);


    // boucle sur le temps pour generé un fichier lisible pour intérpretation graphique
    while (t <= tfinal){
        sortie << t << "\t" << x.NADH_m << "\t" << x.ADP_m << "\t" << x.deltaPsi << "\t" << x.Ca_m << "\t" << Ca_c(t) << "\t" << params["FBP"] << std::endl;
        etat dx_dt = solver(x, params, t, dt);
        x = dx_dt;
        t += dt;
    }

    sortie.close();

    return 0;  
}