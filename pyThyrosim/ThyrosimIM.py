import numpy as np
from scipy.integrate import solve_ivp

def initialize(
    dial = [1.0, 0.88, 1.0, 0.88],
    scale_Vp = True,
    height=1.70,
    weight=70,
    sex=True,               #true = male, false = female,
    fitting_index = [],     # needed in fitting
    p_being_optimized = [], # needed in fitting
    fixed_parameters = [],  # (a, b) means fix p[a] at b 
    scale_plasma_ode = True,
    scale_slow_ode = False,
    scale_fast_ode = False,
    scale_allometric_exponent = False,
    scale_clearance_by_gender = True):


    # initial conditions
    ic       = np.zeros(19)
    ic[0:19] = [0.322114215761171 , 0.201296960359917 , 0.638967411907560, 0.00663104034826483,
                0.0112595761822961, 0.0652960640300348, 1.78829584764370 , 7.05727560072869   ,
                7.05714474742141  , 0                 ,0                 , 0                  ,
                0                 , 3.34289716182018  , 3.69277248068433 , 3.87942133769244   ,
                3.90061903207543 , 3.77875734283571 , 3.55364471589659]

    # Parameter values
    p = np.zeros(100)
    p[0] = 0.0027785399344 #S4 (fitted)
    p[1] = 8               #tau
    p[2] = 0.868           #k12
    p[3] = 0.108           #k13
    p[4] = 584             #k31free
    p[5] = 1503            #k21free
    p[6] = 0.000289        #A
    p[7] = 0.000214        #B
    p[8] = 0.000128        #C
    p[9] = -8.83*10e-6    #D
    p[10] = 0.88           #k4absorb; originally 0.881
    p[11] = 0.0189         #k02
    p[12] = 0.012101809339 #VmaxD1fast (fitted)
    p[13] = 2.85           #KmD1fast
    p[14] = 6.63*10e-4     #VmaxD1slow
    p[15] = 95             #KmD1slow
    p[16] = 0.00074619     #VmaxD2slow
    p[17] = 0.075          #KmD2slow
    p[18] = 3.3572*10e-4   #S3
    p[19] = 5.37           #k45
    p[20] = 0.0689         #k46
    p[21] = 127            #k64free
    p[22] = 2043           #k54free
    p[23] = 0.00395        #a
    p[24] = 0.00185        #b
    p[25] = 0.00061        #c
    p[26] = -0.000505      #d
    p[27] = 0.88           #k3absorb
    p[28] = 0.184972339613 #k05 (fitted)
    p[29] = 450            #Bzero (fixed so max TSH is about 1000)
    p[30] = 219.7085301388 #Azero (fitted)
    p[31] = 0              #Amax (set to 0 because 1976 weeke says hypothyroid patients should have no oscillations)
    p[32] = -3.71          #phi
    p[33] = 0.53           #kdegTSH-HYPO
    p[34] = 0.226          #VmaxTSH (originally it's 0.037 but this is probably a typo because eq4 of 2010 eigenberg it not a real hill function)
    p[35] = 23             #K50TSH
    p[36] = 0.058786935033 #k3 (fitted)
    p[37] = 0.29           #T4P-EU
    p[38] = 0.006          #T3P-EU
    p[39] = 0.037          #KdegT3B
    p[40] = 0.0034         #KLAG-HYPO
    p[41] = 5              #KLAG
    p[42] = 1.3            #k4dissolve
    p[43] = 0.12           #k4excrete; originally 0.119 (change with dial 2)
    p[44] = 1.78           #k3dissolve
    p[45] = 0.12           #k3excrete; originally 0.118 (change with dial 4)
    p[46] = 3.2            #Vp
    p[47] = 5.2            #VTSH

    #parameters for hill functions in f_circ and SRtsh
    p[48] = 3.001011022378 #K_circ (fitted)
    p[49] = 3.094711690204 #K_SR_tsh (fitted)
    p[50] = 5.674773816316 #n, hill exponent in f_circ (fitted)
    p[51] = 6.290803221796 #m, hill exponent in SR_tsh (fitted)
    p[52] = 8.498343729591 #K_f4 for f4 (fitted)
    p[53] = 14.36664496926 #l, hill exponent for f4 (fitted)

    # p[54] = 0.0 # T4 oral dose
    # p[55] = 0.0 # T3 oral dose

    # dial parameters 
    p[56] = dial[0] # controls T4 secretion rate
    p[57] = dial[1] # controls T4 excretion rate
    p[58] = dial[2] # controls T3 secretion rate
    p[59] = dial[3] # controls T3 excretion rate

    # variance parameters for T4/T3/TSH and schneider error (these are used only for parameter estimation!)
    p[60] = 5.003761571969437   # σ for T4 in Blakesley (fixed to reasonable value before fitting)
    p[61] = 0.11122955089297369 # σ for T3 Blakesley and Jonklaas (fixed to reasonable value before fitting)
    p[62] = 0.4                 # σ for TSH in Blakesley and Jonklaas (fixed to reasonable value before fitting)
    p[63] = 0.1                 # σ for FT4 in Jonklaas (fixed to reasonable value before fitting)

    # Blakesley reference BMI
    p[64] = 21.82854404275587 # (male, fitted)
    p[65] = 22.99050845201536 # w / h^2 (female, fitted)

    # Vtsh scaling factor
    p[66] = 1.0 

    # extra parameter
    # p[67] = 22.5 # w / h^2 (female)

    # Volume scaling ratio
    p[68] = 1.0 # Plasma volume ratio
    p[69] = -1.0 # Plasma volume (forgot what this is supposed to represent)
    p[70] = 1.0 # allometric exponent for plasma volume

    # slow compartment scaling ratio
    p[71] = 1.0 # fat-free constant
    p[72] = 0.0 # fat constant
    p[73] = 1.0 # scaling ratio for slow compartment

    # fast compartment scaling ratio
    p[74] = 1.0

    # allometric exponent for k05 
    p[75] = 0.75 # for male (fixed)
    p[76] = 0.75 # for female (fixed, 0.75 works well)

    # ref height for male and femalejul
    p[77] = 1.7608716659237555 # (fitted)
    p[78] = 1.6696106891941103 # (fitted)

    # clearance scale (male / female)
    p[79] = 1.0499391485135692 # male clearance (fitted)

    # infusion parameters
    p[80] = 0.0 # T4 infusion
    p[81] = 0.0 # T3 infusion

    # change fitting parameters
    if len(fitting_index) > 0:
        p[fitting_index] = p_being_optimized

    # scale plasma parameters
    ref_bmi = p[64] if sex else p[65]
    if scale_plasma_ode:
        # for now, assume male and females have the same ref Vp (ie average male/female ref Vp)
        ref_Vp = (reference_Vp(ref_bmi, True, p[77]) + reference_Vp(ref_bmi, False, p[78])) / 2
        p[68] = predict_Vp(height, weight, sex) / ref_Vp

    p[70] = 0.75 if scale_allometric_exponent else p[70]

    # scale slow compartment
    if scale_slow_ode:
        ref_weight = p[64] * p[77]*2 if sex else p[65] * p[78]**2
        ref_fat_free_mass = reference_fat_free_mass(sex, male_ref_height=p[77], female_ref_height=p[78])
        ref_fat_mass = ref_weight - ref_fat_free_mass
        slow_compartment_scale = (p[71] * fat_free_mass(sex, height) + p[72] * (weight - fat_free_mass(sex, height))) / (p[71] * ref_fat_free_mass + p[72] * ref_fat_mass)
        p[73] = slow_compartment_scale


    # scale fast compartment
    p[74] = 1.0 if scale_fast_ode else p[74]

    if scale_Vp:
        Vp, Vtsh = plasma_volume(height, weight, sex, p[66], ref_bmi, p[77], p[78])
        p[46] = Vp
        p[47] = Vtsh

    clearance_allometric_exp = p[75] if sex else p[76]
    if scale_clearance_by_gender:
        ref_weight = p[64] * p[77]**2 if sex else p[65] * p[78]**2
    else:
        ref_weight = (p[64] * p[77]**2 + p[65] * p[78]**2) / 2

    p[28] *= (weight / ref_weight)**clearance_allometric_exp
    p[28] *= p[79] if sex else p[28] # scale male k05 by prefactor

    # fix parameters declared by users
    for (a, b) in fixed_parameters:
        p[a] = b

    return ic, p

def plasma_volume(h, w, sex:bool,
    Vtsh_scale = 1.0, ref_bmi = 22.5,
    male_ref_height = 1.7, female_ref_height=1.63
    ):
    # for now, assume male and females have the same ref Vp (ie average male/female ref Vp)
    ref_Vp = (reference_Vp(ref_bmi, True, male_ref_height) + reference_Vp(ref_bmi, False, female_ref_height)) / 2
    Vp_new = predict_Vp(h, w, sex) * 3.2 / ref_Vp

    # scale Vtsh according to Vtsh_new = Vtsh_old + c(Vp_new - Vp_old) 
    Vtsh_new = 5.2 + Vtsh_scale * (Vp_new - 3.2)

    return Vp_new, Vtsh_new

def reference_Vp(ref_BMI, sex:bool, ref_height):
    # calculate weight for specified ref_BMI. Ideal weight (iw) is fitted to Feldschush's data
    if sex:
        iw = 176.3 - 220.6 * ref_height + 93.5 * ref_height**2
    else:
        iw = 145.8 - 182.7 * ref_height + 79.55 * ref_height**2
  
    w = ref_BMI * ref_height**2

    return predict_Vp(ref_height, w, sex)

def predict_Vp(h, w, sex:bool):
    # hematocrit level, set to .45 for male and .4 for females
    Hem = 0.40 + 0.05 * sex

    # calculate Ideal Weight fitted to Feldschush's data
    if sex:
        iw = 176.3 - 220.6 * h + 93.5 * h**2
    elif not sex:
        iw = 145.8 - 182.7 * h + 79.55 * h**2

    # power law fitted to Feldchush data
    a, n = 1.26975706e+03, 3.72981228e-01
    Δiw = (w - iw) / iw * 100  #deviation from ideal weight, in percentage
    Vb_per_kg = a * (100.0 + Δiw)**(n - 1)
    Vb = Vb_per_kg * w / 1000
    
    return Vb * (1 - Hem)

def blood_volume(h, w, sex:bool):
    Hem = 0.40 + 0.05 * sex #.45 for male and .4 for females (by default)
    BMI = w / h**2

    # calculate Ideal Weight fitted to Feldschush's data
    if sex:
        iw = 176.3 - 220.6 * h + 93.5 * h**2
    elif not sex:
        iw = 145.8 - 182.7 * h + 79.55 * h**2

    # power law fitted to Feldchush data
    a, n = 1.26975706e+03, 3.72981228e-01
    Δiw = (w - iw) / iw * 100  #deviation from ideal weight, in percentage
    Vb_per_kg = a * (100.0 + Δiw)**(n - 1)
    return Vb_per_kg * w / 1000

# a0...a101 are just p. Python argument passing to solveIVP sucks
def thyrosim(t,q,a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80,a81,a82,a83,a84,a85,a86,a87,a88,a89,a90,a91,a92,a93,a94,a95,a96,a97,a98,a99): #changed order for ode_int
    p = [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28,a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49,a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63,a64,a65,a66,a67,a68,a69,a70,a71,a72,a73,a74,a75,a76,a77,a78,a79,a80,a81,a82,a83,a84,a85,a86,a87,a88,a89,a90,a91,a92,a93,a94,a95,a96,a97,a98,a99]
    kdelay = 5/8
    dq = np.zeros(19)

    # scaling the mass/concentration of compartments
    plasma_volume_ratio = p[68]**p[70]
    slow_volume_ratio = p[73]**p[70]
    fast_volume_ratio = p[74]**p[70]

    # scale comparment sizes
    q1 = q[0] * 1 / p[68]
    q2 = q[1] * 1 / p[74]
    q3 = q[2] * 1 / p[73]
    q4 = q[3] * 1 / p[68]
    q5 = q[4] * 1 / p[74]
    q6 = q[5] * 1 / p[73]
    q7 = q[6] * 1 / p[68]

    # Auxillary equations
    q4F = (p[23]+ p[24] * q1 + p[25] * q1**2 + p[26] * q1**3) * q4 #FT3p
    q1F = (p[6] + p[7] * q1 + p[8] * q1**2 + p[9] * q1**3) * q1  #FT4p
    SR3 = (p[18] * p[58] * q[18])                                        #Brain delay (dial 3)
    SR4 = (p[0] * p[56] * q[18])                                         #Brain delay (dial 1)
    fCIRC = q[8]**p[50] / (q[8]**p[50] + p[48]**p[50])
    SRTSH = (p[29]+p[30]*fCIRC*np.sin(np.pi/12*t-p[32]))*(p[49]**p[51]/(p[49]**p[51] + q[8]**p[51]))
    fdegTSH = p[33] + p[34] / (p[35] + q7)
    fLAG = p[40] + 2*q[9]**11 / (p[41]**11 + q[7]**11)
    f4 = p[36]*(1 + 5*(p[52]**p[53]) / (p[52]**p[53]+q[7]**p[53]))
    NL = p[12] / (p[13] + q2)

    # ODEs
    dq[0]  = p[80] + (SR4 + p[2] * q2 + p[3] * q3 - (p[4] + p[5]) * q1F) * plasma_volume_ratio + p[10] * q[10] #T4dot (need to remove u1)
    dq[1]  = (p[5] * q1F - (p[2] + p[11] + NL) * q2) * fast_volume_ratio                                    #T4fast
    dq[2]  = (p[4] * q1F -(p[3] + p[14] / (p[15] + q3) + p[16] /(p[17] + q3)) * q3) * slow_volume_ratio  #T4slow
    dq[3]  = p[81] + (SR3 + p[19] * q5 + p[20] * q6 - (p[21] + p[22]) * q4F) * plasma_volume_ratio + p[27] * q[12] #T3pdot
    dq[4]  = (p[22] * q4F + NL * q2 - (p[19] + p[28]) * q5) * fast_volume_ratio                         #T3fast
    dq[5]  = (p[21] * q4F + p[14] * q3 / (p[15] + q3) + p[16] * q3 / (p[17] + q3) -(p[20])*q6) * slow_volume_ratio #T3slow
    dq[6]  = (SRTSH - fdegTSH * q7) * plasma_volume_ratio                                           #TSHp
    dq[7]  = f4 / p[37] * q1 + p[36] / p[38] * q4 - p[39] * q[7]          #T3B
    dq[8]  = fLAG * (q[7] - q[8])                                             #T3B LAG
    dq[9] = -p[42] * q[9]                                                   #T4PILLdot
    dq[10] =  p[42] * q[9] - (p[43] * p[57]+ p[10]) * q[10]                  #T4GUTdot: note p[44] * p[58] = p[44] * dial[2] = k4excrete
    dq[11] = -p[44] * q[11]                                                   #T3PILLdot
    dq[12] =  p[44] * q[11] - (p[45] * p[59] + p[27]) * q[12]                 #T3GUTdot: note p[46] * p[60] = p[46] * dial[4] = k3excrete

    # Delay ODEs
    dq[13] = kdelay * (q7 - q[13]) 
    dq[14] = kdelay * (q[13] - q[14])                                         #delay2
    dq[15] = kdelay * (q[14] - q[15])                                         #delay3
    dq[16] = kdelay * (q[15] - q[16])                                         #delay4
    dq[17] = kdelay * (q[16] - q[17])                                         #delay5
    dq[18] = kdelay * (q[17] - q[18])                                         #delay6

def output_equations(sol, p):
    return [777.0 * sol[0, :] / p[46], #T4
            651.0 * sol[3, :] / p[46], #T3
            5.6 * sol[6, :] / p[47]] #TSH


def set_patient_ic(ic, p, t4, t3, tsh,
                   steady_state:bool=False, set_tsh_lag:bool=False): #will need to assign variable manually, julia did global
    """
    Set initial conditions from data. Options to set other compartments to steady state,
    optionally including the TSH lag compartments.
    """
    # Set IC for observed compartments. 
    ic[0] = (p[46] * t4) / 777.0
    ic[3] = (p[46] * t3) / 651.0
    ic[6] = (p[47] * tsh) / 5.6
    
    if steady_state:
        q4F = (p[23]+ p[24] * ic[0] + p[25] * ic[0]**2 + p[26] *ic[0]**3) * ic[3] #FT3p
        q1F = (p[6] + p[7] * ic[0] + p[8] * ic[0]**2 + p[9] * ic[0]**3) * ic[0]  #FT4p
        
        B = p[5] * q1F - p[13] * (p[2] + p[11]) - p[12]
        A = -(p[2] + p[11])
        C = p[6] * p[13] * q1F
        ic[1] = (-B - np.sqrt(B**2 - 4.0 * A *C)) / (2.0 * A)
        
        B = p[4] * q1F - (p[3] + p[14] / p[15]) * p[17] - p[16]
        A = -(p[3] + p[14] / p[15])
        C = p[4] * p[17] * q1F
        ic[2] = (-B - np.sqrt(B**2 - 4.0 * A *C)) / (2.0 * A)
        
        ic[4] = (p[22] * q4F + (p[12] / (p[13] + ic[1])) * ic[1]) / (p[19] + p[28])
        ic[5] = (p[21] * q4F + p[14] * (ic[2] / (p[15] + ic[2]))+ p[16] * (ic[2] / (p[17] + ic[2]))) / p[20]
    
    if set_tsh_lag:
        # Probably not 100% correct since they're supposed to be lagged, but probably better than the default.
        ic[14:20] = ic[7]

    return ic

def find_patient_ic(ic, p, days, model = thyrosim): #will need to assign variable manually, julia did global, also sus af
    """
    Find initial conditions from approximate steady state solution. 

    This function runs a Thyrosim simulation for 30 days and sets the initial 
    contidion `ic` to the ending values for each compartment.
    """
    tspan = (0.0, 24*days)
    prob = solve_ivp(model, tspan, ic, args=(p))
    sol = prob.y
    ic = sol[-1]
    return ic

# Figure 2 of Heymsfield 2007: https://academic.oup.com/ajcn/article/86/1/82/4633194
def adipose_tissue_free_mass(sex:bool, h): #true = male, false = female, height h in meter
    h_cm = 100*h
    adp_tfm = 0.0006 * h_cm**2.21 if sex else 0.001 * h_cm**2.1 # unit kg
    return adp_tfm

# Figure 3 of Heymsfield 2007: https://academic.oup.com/ajcn/article/86/1/82/4633194
def fat_free_mass(sex:bool, h): #true = male, false = female, height h in meter
    h_cm = 100*h
    ffm = 0.0004 * h_cm**2.3 if sex else 0.0019 * h_cm**1.97 # unit kg
    return ffm

def reference_fat_free_mass(sex:bool, male_ref_height=1.7, female_ref_height=1.63):
    h = male_ref_height if sex else female_ref_height # avg male/female height
    return fat_free_mass(sex, h) # unit kg

# Table 2 of Muler 2011: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0022732
def liver(w, h): # w in kg, h in meter
    return (0.088 * w**0.54 * h**1.04)

def kidney(w, h):
    return 0.012 * w**0.72 * h**0.19
