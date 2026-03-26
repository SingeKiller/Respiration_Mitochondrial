import numpy as np 
import matplotlib.pyplot as plt
from scipy import stats
import subprocess
from pathlib import Path

CPP_READY = True

dossier_source = "normalized_respiration_data"
courbe_1 = np.load(f"{dossier_source}/data14h16.npy")
courbe_2 = np.load(f"{dossier_source}/data15h12.npy")
courbe_3 = np.load(f"{dossier_source}/data15h37.npy")
courbe_4 = np.load(f"{dossier_source}/data16h00.npy")
courbe_5 = np.load(f"{dossier_source}/data17h03.npy")



# séparation en section pour le lin reg

t_1 = [[14.6 , 15.5], [16.25 , 16.8],[17.6 , 18.5],[19.7 , 20.5]]
t_2 = [[10.8, 12.8], [13.3 , 14.35], [16, 17],[17.25 , 17.8]]
t_3 = [[10.6 , 11.3],[12.3, 12.6], [13.1 , 13.55], [14.3 , 15.7]]
t_4 = [[10.75 , 11.75], [12.5 , 13.75], [14.25, 15], [15.4, 16.3]]
t_5 = [[10.3, 11.45], [12.2 , 12.85], [13.75 , 14.8], [15.4, 17]]

def decoupage_reg(courbe , decoupage, affichage = True):

    droite_a_b= []
    #affichage
    if affichage :
        plt.subplot(1,len(decoupage)+1,1)
        plt.plot(courbe[0,:], courbe[1,:], color= 'blue', linestyle = ':')


    for i in range (len(decoupage)):
        #calcul pente
        masque = (courbe[0, :] > decoupage[i][0]) & (courbe[0, :] < decoupage[i][1])
        reg = stats.linregress(courbe[0,:][masque], courbe[1,:][masque] )
        droite = reg.slope *courbe[0,:][masque]+ reg.intercept
        #sauvegarde des equation
        droite_a_b.append([float(reg.slope) , float(reg.intercept)])
        #affichage
        if affichage :
            plt.subplot(1,len(decoupage)+1,1)
            plt.plot(courbe[0,:][masque], droite, color = "red" )
            plt.subplot(1,len(decoupage)+1,i+2)
            plt.plot(courbe[0,:][masque], droite, color = "red" )
            plt.plot(courbe[0,:][masque], courbe[1,:][masque], color = "blue" , linestyle=":" )

            
    if affichage:
        plt.show()     
    return droite_a_b

droite = []
droite.append(decoupage_reg(courbe_1 , t_1, affichage= False))
droite.append(decoupage_reg(courbe_2 , t_2, affichage= False))
droite.append(decoupage_reg(courbe_3 , t_3, affichage= False))
droite.append(decoupage_reg(courbe_4 , t_4, affichage= False))
droite.append(decoupage_reg(courbe_5 , t_5, affichage= False))


if CPP_READY :
    cpp_dir = Path(__file__).resolve().parents[1] / "Cpp"
    programme_cpp = cpp_dir / "main.exe"
    if not programme_cpp.exists():
        programme_cpp = cpp_dir / "main"

    cpp_args = [
        "400", "1", "0.01", "0.6", "100", "177", "5", "7", "100", "177",
        "5", "120", "10000", "190", "8.5", "35", "0.002", "-0.03", "0.35",
        "2", "0.01", "1.1", "0.001", "0.016", "0.037", "10000", "15000",
        "1.8", "0.01", "0.0005", "500", "15000", "5", "0.1", "1", "25",
        "16", "330", "19", "1000", "o2_pulses_2p5_10_normalized.txt"
    ]

    subprocess.run([str(programme_cpp), *cpp_args], check=True, cwd=str(cpp_dir))


# il faut récuperer les différent fichier et prendre la reg lineaire qui ets la plus proche pour cela 
"""
    recuperer les point
    lancer la fonction decoupage_reg(courbe , bon_decoupage)
    passer dans un focntion comparaison qui arde la meilleure
        la fonction comparaison(reg_ref , reg_cpp )
            erreur slope = abs(a_ref-reg_cpp) 
            erreur b = abs(b_ref - b_reg)
            return erreur slope et erreur b
    
    dans un boucle on lance le programme et on garde l'erreur la plus petite 
"""
