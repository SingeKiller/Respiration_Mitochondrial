import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode

# --- DÉFINITION DES TAILLES ---
size_algebrica = 33
size_state = 7
size_constant = 47

# --- DÉFINITION DES INDICES (CONSTANTES) ---
indice_constant = list(range(size_constant))
gamma, NADtot, fm, Cmito, Am_tot, Ac_tot, khyd, Jhydbas, Fhold, Ftest, Fton, Ftoff, Chold, Ctest, Cton1, Cton2, Cton3, Ctoff1, Ctoff2, Ctoff3, p1, p2, p3, p21, p22, p23, p24, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, FRT, kgpdh, delta = indice_constant 

# --- DÉFINITION DES INDICES (ALGÉBRIQUES) ---
indice_algebrica = list(range(size_algebrica))
NADm, ATPM, heav_on, heav_Cton1, MM1, RATm, heav_off, heav_Cton2, Jo, Atpc, Fproto, heav_Cton3, Jgpdh, Jh_yd, heav_Ctoff1, Jpdh, Juni, pulse1, Jnaca, heav_Ctoff2, Jh_leak, Jmito, pulse2, MM2, heav_Ctoff3, Jh_res, pulse3, b13, Cproto, Jh_atp, b2, Jf1f0, Jant = indice_algebrica

# --- DÉFINITION DES INDICES (ÉTAT) ---
indice_state = list(range(size_state))
NADHm, ADPm, PSIm, Cam, c, ADPc, FBP = indice_state

def init_constant():
    constant = [0.0]*size_constant

    # Mitochondrial_constante
    constant[gamma] =  0.001 
    constant[NADtot] = 10.0 
    constant[fm] = 0.01 
    constant[Cmito] = 1.8 
    constant[Am_tot] = 15.0

    # Cytosol_constante
    constant[Ac_tot] = 2500.0
    constant[khyd] = 0.00005 
    constant[Jhydbas] = 0.00005
    constant[delta] = 3.90000/53.2000

    # clamp_constant (FIGURE 12 CONFIG)
    constant[Fhold] = 1.0 
    constant[Ftest] = 5.0 
    constant[Fton] = 90000.0 # 1.5 min * 60 * 1000 = 90000 ms
    constant[Ftoff] = 330000.0 # 5.5 min * 60 * 1000 = 330000 ms
    
    constant[Chold] = 0.1
    constant[Ctest] = 0.1
    constant[Cton1] = 120000.0
    constant[Cton2]= 180000.0
    constant[Cton3] = 240000.0
    constant[Ctoff1] = 150000.0
    constant[Ctoff2] = 210000.0
    constant[Ctoff3] = 270000.0

    # parametre_constant 
    constant[p1] = 400.0
    constant[p2] = 1.0
    constant[p3] = 0.01 
    constant[p21] = 0.01 
    constant[p22] = 1.1 
    constant[p23] = 0.001 
    constant[p24] = 0.016 
    constant[p4] =  0.6 
    constant[p5] = 0.1 
    constant[p6] = 177.0 
    constant[p7] = 5.0  
    constant[p8] = 7.0 
    constant[p9] = 0.1 
    constant[p10] = 177.0 
    constant[p11] = 5.0 
    constant[p12] = 120.0 
    constant[p13] = 10.0 
    constant[p14] = 190.0 
    constant[p15] = 8.5 
    constant[p16] = 35.0 
    constant[p17] = 0.002 
    constant[p18] = -0.03 
    constant[p19] = 0.35 
    constant[p20] = 2.0 

    # j_variable_constant
    constant[FRT] = 0.037 
    constant[kgpdh] = 0.0005 

    return constant

def init_state():
    state = [0.0]*size_state
    # mito
    state[NADHm] = 0.4 #(mM) 
    state[ADPm] = 12.0 #(mM)
    state[PSIm] = 150.0 #(mV)
    state[Cam] = 0.05 #(uM)

    # clamp
    state[c] = 0.1 #(uM)
    state[FBP] = 1.0 #(uM)
    # cytosol
    state[ADPc] = 1850.0 #(uM)

    return state

def custom_piecewise(cases):
    """Compute result of a piecewise function using Numpy Select"""
    return np.select(cases[0::2], cases[1::2])

def calcul_rate(voi, state, constant):
    algebrica = [0.0]*size_algebrica
    rate = [0.0]*size_state
    
    # --- SECURITE NUMERIQUE CRITIQUE ---
    # Utilisation de np.maximum pour éviter l'erreur "cannot be interpreted as integer"
    # Cela fonctionne aussi bien sur les scalaires que les tableaux
    NADHm_val = np.maximum(1e-12, state[NADHm])
    ADPm_val  = np.maximum(1e-12, state[ADPm])
    PSIm_val  = state[PSIm] # Voltage peut etre negatif
    Cam_val   = np.maximum(1e-12, state[Cam])
    c_val     = np.maximum(1e-12, state[c])
    ADPc_val  = np.maximum(1e-12, state[ADPc])
    FBP_val   = np.maximum(0.0, state[FBP])
    # -----------------------------------

    # calcul évolution FBP
    # Utilisation de np.greater_equal pour la compatibilité
    algebrica[heav_on] = custom_piecewise([np.greater_equal(voi - constant[Fton], 0.0), 1.0, True, 0.0])
    algebrica[heav_off] = custom_piecewise([np.greater_equal(voi - constant[Ftoff], 0.0), 1.0, True, 0.0])
    algebrica[Fproto] = constant[Fhold] + (constant[Ftest] - constant[Fhold]) * (algebrica[heav_on] - algebrica[heav_off])
    rate[FBP] = (algebrica[Fproto] - FBP_val)/100.0 

    # calcul évolution NADHm
    algebrica[NADm] = constant[NADtot] - NADHm_val
    algebrica[Jgpdh] = constant[kgpdh] * (np.sqrt(FBP_val)) 
    
    # Calculs intermédiaires sécurisés
    term_pdh_1 = constant[p1] / (constant[p2] + NADHm_val/algebrica[NADm])
    term_pdh_2 = Cam_val / (constant[p3] + Cam_val)
    algebrica[Jpdh] = term_pdh_1 * term_pdh_2 * algebrica[Jgpdh]
    
    algebrica[MM1] = (constant[p4]*NADHm_val)/(constant[p5]+NADHm_val)
    algebrica[Jo] =  algebrica[MM1]/(1.0 + np.exp((PSIm_val-constant[p6])/constant[p7]))
    rate[NADHm] = constant[gamma] * (algebrica[Jpdh] - algebrica[Jo])

    # calcul évolution Cam
    algebrica[Jnaca] = ((constant[p23]*Cam_val)/c_val) * np.exp(constant[p24]*PSIm_val) 
    algebrica[Juni] = (constant[p21]*PSIm_val - constant[p22]) * (c_val**2) 
    algebrica[Jmito] = algebrica[Jnaca] - algebrica[Juni]
    rate[Cam] = -1 * constant[fm] * algebrica[Jmito]

    # calcul évolution c (avec 3 pulse)
    algebrica[heav_Cton1] = custom_piecewise([np.greater_equal(voi - constant[Cton1], 0.0), 1.0, True, 0.0])
    algebrica[heav_Ctoff1] = custom_piecewise([np.greater_equal(voi - constant[Ctoff1], 0.0), 1.0, True, 0.0])
    algebrica[pulse1] = (constant[Ctest] - constant[Chold])*(algebrica[heav_Cton1] - algebrica[heav_Ctoff1]) 
    
    algebrica[heav_Cton2] = custom_piecewise([np.greater_equal(voi - constant[Cton2], 0.0), 1.0, True, 0.0])
    algebrica[heav_Ctoff2] = custom_piecewise([np.greater_equal(voi - constant[Ctoff2], 0.0), 1.0, True, 0.0])
    algebrica[pulse2] = (constant[Ctest]-constant[Chold])*(algebrica[heav_Cton2] - algebrica[heav_Ctoff2]) 

    algebrica[heav_Cton3] = custom_piecewise([np.greater_equal(voi - constant[Cton3], 0.0), 1.0, True, 0.0])
    algebrica[heav_Ctoff3] =  custom_piecewise([np.greater_equal(voi - constant[Ctoff3], 0.0), 1.0, True, 0.0])
    algebrica[pulse3] = (constant[Ctest]-constant[Chold])*(algebrica[heav_Cton3] - algebrica[heav_Ctoff3]) 

    algebrica[Cproto] = constant[Chold]+ algebrica[pulse1] + algebrica[pulse2] + algebrica[pulse3]
    rate[c] = (algebrica[Cproto] - c_val)/100.0 

    # calcul évolution ADPm 
    algebrica[ATPM] = constant[Am_tot] - ADPm_val
    algebrica[RATm] = algebrica[ATPM] / ADPm_val
    algebrica[b2] = (constant[p16]*constant[p13])/(constant[p13] + algebrica[ATPM]) 
    algebrica[Jf1f0] = algebrica[b2]/(1.0 + np.exp((constant[p14] - PSIm_val)/constant[p15])) 
    algebrica[Jant] = ((constant[p19]* algebrica[RATm])/(algebrica[RATm] + constant[p20])) * np.exp(0.5*constant[FRT]*PSIm_val)
    rate[ADPm] = constant[gamma] * (algebrica[Jant] - algebrica[Jf1f0])

    # calcul évolution PSIm
    algebrica[MM2] = (constant[p8]*NADHm_val)/(constant[p9]+NADHm_val) 
    algebrica[Jh_res] = algebrica[MM2]/(1.0 + np.exp((PSIm_val - constant[p10])/constant[p11])) 
    algebrica[b13] =  (constant[p12]*constant[p13])/(constant[p13] + algebrica[ATPM])
    algebrica[Jh_atp] =  algebrica[b13]/(1.0 + np.exp((constant[p14] - PSIm_val)/constant[p15])) 
    algebrica[Jh_leak] = constant[p17]*PSIm_val + constant[p18] 
    
    # Equation du potentiel (Eq 5)
    ddPsi = (algebrica[Jh_res] - algebrica[Jh_atp] - algebrica[Jant] - algebrica[Jh_leak] - algebrica[Jnaca] - 2*algebrica[Juni])
    rate[PSIm] = ddPsi / constant[Cmito]

    # calcul évolution ADPc 
    algebrica[Atpc] = constant[Ac_tot] - ADPc_val 
    algebrica[Jh_yd] = ((constant[khyd]*c_val+constant[Jhydbas])*algebrica[Atpc])
    rate[ADPc] = -constant[delta]*algebrica[Jant]+algebrica[Jh_yd]*1.00000
    return rate

def calcul_algebrica(voi, state, constant):
    # Pour le post-traitement (graphiques), state est une matrice
    algebrica = np.array([[0.0] * len(voi)] * size_algebrica)
    state = np.array(state)
    voi = np.array(voi)

    # --- SECURITE TABLEAUX ---
    NADHm_val = np.maximum(1e-12, state[NADHm])
    ADPm_val  = np.maximum(1e-12, state[ADPm])
    PSIm_val  = state[PSIm]
    Cam_val   = np.maximum(1e-12, state[Cam])
    c_val     = np.maximum(1e-12, state[c])
    ADPc_val  = np.maximum(1e-12, state[ADPc])
    FBP_val   = np.maximum(0.0, state[FBP])
    # -------------------------

    algebrica[heav_on] = custom_piecewise([np.greater_equal(voi- constant[Fton], 0.0), 1.0, True, 0.0])
    algebrica[heav_off] = custom_piecewise([np.greater_equal(voi-constant[Ftoff] , 0.0), 1.0 , True, 0.0])
    algebrica[Fproto] = constant[Fhold] + (constant[Ftest] - constant[Fhold]) * (algebrica[heav_on] - algebrica[heav_off])

    algebrica[NADm] = constant[NADtot] - NADHm_val
    algebrica[Jgpdh] =  constant[kgpdh] * np.sqrt(FBP_val)
    algebrica[Jpdh] =(constant[p1]/(constant[p2]+NADHm_val/algebrica[NADm]))*(Cam_val/(constant[p3]+Cam_val))*algebrica[Jgpdh]
    algebrica[MM1] = (constant[p4]*NADHm_val)/(constant[p5]+NADHm_val)
    algebrica[Jo] =  algebrica[MM1]/(1.0 + np.exp((PSIm_val-constant[p6])/constant[p7]))

    algebrica[Jnaca] = ((constant[p23]*Cam_val)/c_val)*np.exp(constant[p24]*PSIm_val) 
    algebrica[Juni] = (constant[p21]*PSIm_val - constant[p22])*(c_val**2) 
    algebrica[Jmito] = algebrica[Jnaca] - algebrica[Juni]

    return algebrica

def solve_model(fbp):
    constant = init_constant()
    state = init_state()
    constant[Ftest] = fbp
    # Timeline en ms (7 minutes = 420000 ms)
    voi = np.linspace(0, 420000, 5000) 
    
    r = ode(calcul_rate)
    r.set_integrator('vode', method='bdf', atol=1e-06, rtol=1e-06, max_step=50.0) 
    r.set_initial_value(state, voi[0])
    r.set_f_params(constant)

    # Matrice de resultats
    state_res = np.zeros((size_state, len(voi)))
    state_res[:, 0] = state

    print("Debut simulation...")
    for i in range(1, len(voi)):
        r.integrate(voi[i])
        if r.successful():
            state_res[:, i] = r.y
        else:
            print(f"Erreur integration a t={voi[i]}")
            break
            
    print("Simulation terminee.")
    algebrica_res = calcul_algebrica(voi, state_res, constant)
    return (voi, state_res, algebrica_res)

if __name__ == "__main__":
    
    # Valeurs de FBP testées dans l'article
    fbp_levels = [5.0, 10.0, 15.0]
    # Styles de ligne pour correspondre à la figure originale
    styles = [':', '--', '-'] 
    labels = ['5 uM', '10 uM', '15 uM']

    # Configuration de la grille de 6 graphiques (3 lignes, 2 colonnes)
    fig, axes = plt.subplots(3, 2, figsize=(12, 12))
    plt.subplots_adjust(hspace=0.4, wspace=0.3)

    for val, sty, lab in zip(fbp_levels, styles, labels):
        print(f"Simulation pour FBP = {val} uM...")
        voi, state_res, algebrica_res = solve_model(val)
        
        # Le temps en minutes (ignorer t=100)
        t_min = voi[100:] / 60000.0
        idx = slice(100, None) 

        # (A) FBP - Input
        axes[0, 0].plot(t_min, state_res[FBP, idx], 'k', linestyle=sty)
        axes[0, 0].set_ylabel('FBP (uM)')
        axes[0, 0].set_title('(A)')

        # (B) Cam - Calcium mitochondrial
        axes[0, 1].plot(t_min, state_res[Cam, idx], 'k', linestyle=sty)
        axes[0, 1].set_ylabel('Cam (uM)')
        axes[0, 1].set_title('(B)')
        axes[0, 1].set_ylim(0, 0.3)

        # (C) NADHm - NADH mitochondrial
        axes[1, 0].plot(t_min, state_res[NADHm, idx], 'k', linestyle=sty)
        axes[1, 0].set_ylabel('NADHm (mM)')
        axes[1, 0].set_title('(C)')

        # (D) PSIm - Potentiel de membrane
        axes[1, 1].plot(t_min, state_res[PSIm, idx], 'k', linestyle=sty)
        axes[1, 1].set_ylabel('Delta Psi (mV)')
        axes[1, 1].set_title('(D)')

        # (E) Jo - Consommation d'Oxygène
        axes[2, 0].plot(t_min, algebrica_res[Jo, idx], 'k', linestyle=sty)
        axes[2, 0].set_ylabel('Jo (uM/ms)')
        axes[2, 0].set_xlabel('Time (min)')
        axes[2, 0].set_title('(E)')

        # (F) ATPm - ATP mitochondrial (Calculé par conservation)
        # Am_tot = 15mM d'après le Tableau 1
        atpm = 15.0 - state_res[ADPm, idx]
        axes[2, 1].plot(t_min, atpm, 'k', linestyle=sty, label=lab)
        axes[2, 1].set_ylabel('ATPm (mM)')
        axes[2, 1].set_xlabel('Time (min)')
        axes[2, 1].set_title('(F)')

    # Ajout d'une légende unique pour le dernier panneau
    axes[2, 1].legend()
    
    print("Affichage des graphiques...")
    plt.show()