# -*- coding: utf-8 -*-
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
from itertools import chain

def dH(A,B,C,D,E,F,H): # kJ/(mol*K)
    dH_MOLECULE = A*t + B*t**2/2 + C*t**3/3 + D*t**4/4 - E/t + F - H
    return dH_MOLECULE

def s(A,B,C,D,E,G): # J/(mol*K)
    s_MOLECULE = A*np.log(t) + B*t + C*t**2/2 + D*t**3/3 - E/(2*t**2) + G
    return s_MOLECULE
    
def G1(H,t,S):
    g = H - t*S
    return g

def dG(dH,t,dS):
    dG = dH - t*dS
    return dG

def dG2(G,G_28915):
    dG = G - G_28915
    return dG
    
# enthalpy and free energy of reaction at 298.15 K for one reaction (e.g. CO2 + H2 <-> CO + H2O)
    
def Hrxn_29815(Hf_29815_reactant1,Hf_29815_reactant2,Hf_29815_product1,Hf_29815_product2):
    Hrxn_29815 = Hf_29815_product1 + Hf_29815_product2 - Hf_29815_reactant1 - Hf_29815_reactant2
    return Hrxn_29815
    
def Srxn_29815(S_29815_reactant1,S_29815_reactant2,S_29815_product1,S_29815_product2):
    Srxn_29815 = S_29815_product1 + S_29815_product2 - S_29815_reactant1 - S_29815_reactant2
    return Srxn_29815
    
def Grxn_29815(Hrxn_29815,Srxn_29815):
    Grxn_29815 = Hrxn_29815 - 298.15*(Srxn_29815)/1000
    return Grxn_29815
    
def abundances_norm_H(n_C_wanted):
    
    #METHANE
    
    T = np.linspace(500,3000,26)
    count_1bar=0
    count_10mbar=0
    
    n_C2H2_1bar = []
    n_H2O_1bar = []
    n_CO_1bar = []
    n_CO2_1bar = []
    
    n_C2H2_10mbar = []
    n_H2O_10mbar = []
    n_CO_10mbar = []
    n_CO2_10mbar = []
    for j in range(2):
        A=[]
        if j==0:
            print("P = 1 bar")
            K1 = K1_1bar
            K2 = K2_sin_presion
            K3 = K3_1bar
            n_CH4_1bar = []
        if j==1:
            print("P = 10 mbar")
            K1 = K1_10mbar
            K2 = K2_sin_presion
            K3 = K3_10mbar
            n_CH4_10mbar = []
        A.append(-2*n_C_wanted)
        A.append(8*(K1/K2)*(n_O-n_C_wanted)**2 + 1 + 2*K1*(n_O-n_C_wanted))
        A.append(8*(K1/K2)*(n_O-n_C_wanted) + 2*K3 + K1)
        A.append(2*(K1/K2)*(1 + 8 * K3 * (n_O-n_C_wanted)) + 2*K1*K3)
        A.append(8*(K1*K3/K2))
        A.append(8*(K1*K3**2/K2))
        for i in range(26):
            print("Coeficientes (A)")
            print(A[0],'|',A[1][i],'|',A[2][i],'|',A[3][i],'|',A[4][i],'|',A[5][i])
            roots = poly.polyroots([A[0],A[1][i],A[2][i],A[3][i],A[4][i]])
            print("Raices para T=",T[i])
            print(roots)
            for w in range(len(roots)):
                two_n_C = float(-A[0]) # We know that n_CH4 should be approx. 2*n_C
                possible_n_CH4 = roots[w]
                print(possible_n_CH4)
                print(two_n_C)
                if ('{0:.1f}'.format(possible_n_CH4.real) == '{0:.1f}'.format(two_n_C) or (possible_n_CH4.real <= two_n_C and possible_n_CH4.real > 0.0)) and possible_n_CH4.imag == 0.0 :
                    print("~~~~~ MATCH 1313 ~~~~~(T=",T[i],")","possible_n_CH4.real=",'{0:.1f}'.format(possible_n_CH4.real),"two_n_C=",'{0:.1f}'.format(two_n_C))
                    if j==0:
                        n_CH4_1bar.append(possible_n_CH4.real)
                        count_1bar += 1 # This number should be 26 at the end of the loops (because we have 26 temperatures)
                    if j==1:
                        n_CH4_10mbar.append(possible_n_CH4.real)
                        count_10mbar += 1
    print("***count_1bar*** = ",count_1bar)
    print("***count_10mbar*** = ",count_10mbar)
    
    print("# METHANE")
    
    print(n_CH4_1bar,"|",len(n_CH4_1bar))
    print(n_CH4_10mbar,"|",len(n_CH4_10mbar))
    
    # ACETYLENE
    print("# ACETYLENE")

    n_C2H2_1bar = K3_1bar * np.array(n_CH4_1bar)**2.0
    n_C2H2_10mbar = K3_10mbar * np.array(n_CH4_10mbar)**2.0
    
    print(n_C2H2_1bar,"|",len(n_C2H2_1bar))
    print(n_C2H2_10mbar,"|",len(n_C2H2_10mbar))

    # WATER
    print("# WATER")

    n_H2O_1bar = (2.0 * K3_1bar * np.array(n_CH4_1bar)**2.0) + np.array(n_CH4_1bar) + (2.0 * (n_O - n_C_wanted))
    n_H2O_10mbar = (2.0 * K3_10mbar * np.array(n_CH4_10mbar)**2.0) + np.array(n_CH4_10mbar) + (2.0 * (n_O - n_C_wanted))

    print(n_H2O_1bar,"|",len(n_H2O_1bar))
    print(n_H2O_10mbar,"|",len(n_H2O_10mbar))

    # CARBON MONOXIDE
    print("# CARBON MONOXIDE")

    n_CO_1bar = K1_1bar * n_CH4_1bar * n_H2O_1bar
    n_CO_10mbar = K1_10mbar * n_CH4_10mbar * n_H2O_10mbar

    print(n_CO_1bar,"|",len(n_CO_1bar))
    print(n_CO_10mbar,"|",len(n_CO_10mbar))

    # CARBON DIOXIDE
    print("# CARBON DIOXIDE")

    n_CO2_1bar = n_CO_1bar * n_H2O_1bar /K2_sin_presion
    n_CO2_10mbar = n_CO_10mbar * n_H2O_10mbar /K2_sin_presion
    
    print(n_CO2_1bar,"|",len(n_CO2_1bar))
    print(n_CO2_10mbar,"|",len(n_CO2_10mbar))
    
    return n_CH4_1bar, n_CH4_10mbar, n_C2H2_1bar, n_C2H2_10mbar, n_H2O_1bar, n_H2O_10mbar, n_CO_1bar, n_CO_10mbar, n_CO2_1bar, n_CO2_10mbar
    
def graphics_CO(P,figN,CO_wanted): # P = "1bar" o "10mbar" o ""
    plt.figure()
    if P == "1bar":
        plt.plot(T,n_CH4_1bar,label='$CH_{4} (1bar)$',color='black',ls=':',linewidth=3)
        plt.plot(T,n_C2H2_1bar,label='$C_{2}H_{2} (1bar)$',color='y',ls=':',linewidth=3)
        plt.plot(T,n_H2O_1bar,label='$H_{2}O (1bar)$',color='m',ls=':',linewidth=3)
        plt.plot(T,n_CO_1bar,label='$CO (1bar)$',color='c',ls=':',linewidth=3)
        plt.plot(T,n_CO2_1bar,label='$CO_{2} (1bar)$',color='g',ls=':',linewidth=3)
    if P == "10mbar":
        plt.plot(T,n_CH4_10mbar,label='$CH_{4} (10mbar)$',color='black',ls=':',linewidth=1)
        plt.plot(T,n_C2H2_10mbar,label='$C_{2}H_{2} (10mbar)$',color='y',ls=':',linewidth=1)
        plt.plot(T,n_H2O_10mbar,label='$H_{2}O (10mbar)$',color='m',ls=':',linewidth=1)
        plt.plot(T,n_CO_10mbar,label='$CO (10mbar)$',color='c',ls=':',linewidth=1)
        plt.plot(T,n_CO2_10mbar,label='$CO_{2} (10mbar)$',color='g',ls=':',linewidth=1)
    if P == "":
        plt.plot(T,n_CH4_1bar,label='$CH_{4} (1bar)$',color='black',ls=':',linewidth=3)
        plt.plot(T,n_C2H2_1bar,label='$C_{2}H_{2} (1bar)$',color='y',ls=':',linewidth=3)
        plt.plot(T,n_H2O_1bar,label='$H_{2}O (1bar)$',color='m',ls=':',linewidth=3)
        plt.plot(T,n_CO_1bar,label='$CO (1bar)$',color='c',ls=':',linewidth=3)
        plt.plot(T,n_CO2_1bar,label='$CO_{2} (1bar)$',color='g',ls=':',linewidth=3)
        plt.plot(T,n_CH4_10mbar,label='$CH_{4} (10mbar)$',color='black',ls='-',linewidth=1)
        plt.plot(T,n_C2H2_10mbar,label='$C_{2}H_{2} (10mbar)$',color='y',ls='-',linewidth=1)
        plt.plot(T,n_H2O_10mbar,label='$H_{2}O (10mbar)$',color='m',ls='-',linewidth=1)
        plt.plot(T,n_CO_10mbar,label='$CO (10mbar)$',color='c',ls='-',linewidth=1)
        plt.plot(T,n_CO2_10mbar,label='$CO_{2} (10mbar)$',color='g',ls='-',linewidth=1)
    plt.xlim([500, 3000])
    plt.ylim([10**-22, 10**0])
    plt.yscale('log')
    plt.xlabel('Temperature (K)')
    plt.ylabel('$ñ_{x}$')
    plt.title('C/O = '+str(CO_wanted))
    #plt.legend( loc='best')
    plt.legend(loc=4,prop={'size':7})
    plt.savefig('wgs-nist-'+str(figN)+'.png')
    plt.show()


    
T = np.linspace(500,1700,13) # degrees K
t = T/1000
    
#       H2O
       
# 500-1700 K valid temperature range
A =   30.09200
B =   6.832514
C =   6.793435
D =  -2.534480
E =   0.082139
F =  -250.8810
G =   223.3967
H =  -241.8264

Hf_29815_H2O = -241.83 # this is Hf (kJ/(mol*K))
s_29815_H2O = 188.84 # (J/(mol*K))


dH_H2O_1 = dH(A,B,C,D,E,F,H)
h_H2O_1 = dH_H2O_1 + Hf_29815_H2O
s_H2O_1 = s(A,B,C,D,E,G)
ds_H2O_1 = s_H2O_1 - s_29815_H2O
g_H2O_1 = G1(h_H2O_1,t,s_H2O_1)
dg_H2O_1 = dG(dH_H2O_1,t,ds_H2O_1)

T = np.linspace(1800,3000,13) # degrees K
t = T/1000

# 1700-6000 K valid temperature range
A =   41.96426
B =   8.622053
C =  -1.499780
D =   0.098119
E =  -11.15764
F =  -272.1797
G =   219.7809
H =  -241.8264

dH_H2O_2 = dH(A,B,C,D,E,F,H)
h_H2O_2 = dH_H2O_2 + Hf_29815_H2O
s_H2O_2 = s(A,B,C,D,E,G)
ds_H2O_2 = s_H2O_2 - s_29815_H2O
g_H2O_2 = G1(h_H2O_2,t,s_H2O_2)
dg_H2O_2 = dG(dH_H2O_2,t,ds_H2O_2)

dH_H2O = np.array(list(dH_H2O_1)+list(dH_H2O_2))
h_H2O = np.array(list(h_H2O_1)+list(h_H2O_2))
s_H2O = np.array(list(s_H2O_1)+list(s_H2O_2))
ds_H2O = np.array(list(ds_H2O_1)+list(ds_H2O_2))
g_H2O = np.array(list(g_H2O_1)+list(g_H2O_2))
dg_H2O = np.array(list(dg_H2O_1)+list(dg_H2O_2))



print("\n                                            WATER (H2O)")
print("      Temperature |      dH       |       H       |       dS      |       S       |       G       |       *dG*      ")
print("     ---------------------------------------------------------------------------------------------------------------")
for i in np.arange(26):
    print("\t",np.linspace(500,3000,26)[i],"  |",dH_H2O[i],"|",h_H2O[i],"|",ds_H2O[i],"|",s_H2O[i],"|",g_H2O[i],"|",dg_H2O[i])

#       CH4

T = np.linspace(500,1300,9) # degrees K
t = T/1000
       
# 298-1300 K valid temperature range
A =  -0.73029
B =   108.4773
C =  -42.52157
D =   5.862788
E =   0.678565
F =  -76.84376
G =   158.7163
H =  -74.87310

Hf_29815_CH4 = -74.87 #this is Hf.
s_29815_CH4 = 186.25 #P = 1bar

dH_CH4_1 = dH(A,B,C,D,E,F,H)
h_CH4_1 = dH_CH4_1 + Hf_29815_CH4
s_CH4_1 = s(A,B,C,D,E,G)
ds_CH4_1 = s_CH4_1 - s_29815_CH4
g_CH4_1 = G1(h_CH4_1,t,s_CH4_1)
dg_CH4_1 = dG(dH_CH4_1,t,ds_CH4_1)

T = np.linspace(1400,3000,17) # degrees K
t = T/1000

# 1300-6000 K valid temperature range
A =   85.81217
B =   11.26467
C =  -2.114146
D =   0.138190
E =  -26.42221
F =  -153.5327
G =   224.4143
H =  -74.87310

dH_CH4_2 = dH(A,B,C,D,E,F,H)
h_CH4_2 = dH_CH4_2 + Hf_29815_CH4
s_CH4_2 = s(A,B,C,D,E,G)
ds_CH4_2 = s_CH4_2 - s_29815_CH4
g_CH4_2 = G1(h_CH4_2,t,s_CH4_2)
dg_CH4_2 = dG(dH_CH4_2,t,ds_CH4_2)

dH_CH4 = np.array(list(dH_CH4_1)+list(dH_CH4_2))
h_CH4 = np.array(list(h_CH4_1)+list(h_CH4_2))
s_CH4 = np.array(list(s_CH4_1)+list(s_CH4_2))
ds_CH4 = np.array(list(ds_CH4_1)+list(ds_CH4_2))
g_CH4 = np.array(list(g_CH4_1)+list(g_CH4_2))
dg_CH4 = np.array(list(dg_CH4_1)+list(dg_CH4_2))

print("\n                                            METHANE (CH4)")
print("\n      Temperature |      dH       |       H       |       dS      |       S       |       G       |       *dG*      ")
print("     ---------------------------------------------------------------------------------------------------------------")

for i in np.arange(26):
    print("\t",np.linspace(500,3000,26)[i],"  |",dH_CH4[i],"|",h_CH4[i],"|",ds_CH4[i],"|",s_CH4[i],"|",g_CH4[i],"|",dg_CH4[i])

#       CO

T = np.linspace(500,1300,9) # degrees K
t = T/1000
       
# 298-1300 K valid temperature range
A =   25.56759
B =   6.096130
C =   4.054656
D =  -2.671301
E =   0.131021
F =  -118.0089
G =   227.3665
H =  -110.5271

Hf_29815_CO = -110.53 #this is Hf.
s_29815_CO = 197.66 #P = 1bar

dH_CO_1 = dH(A,B,C,D,E,F,H)
h_CO_1 = dH_CO_1 + Hf_29815_CO
s_CO_1 = s(A,B,C,D,E,G)
ds_CO_1 = s_CO_1 - s_29815_CO
g_CO_1 = G1(h_CO_1,t,s_CO_1)
dg_CO_1 = dG(dH_CO_1,t,ds_CO_1)

T = np.linspace(1400,3000,17) # degrees K
t = T/1000

# 1300-6000 K valid temperature range
A =   35.15070
B =   1.300095
C =  -0.205921
D =   0.013550
E =  -3.282780
F =  -127.8375
G =   231.7120
H =  -110.5271

dH_CO_2 = dH(A,B,C,D,E,F,H)
h_CO_2 = dH_CO_2 + Hf_29815_CO
s_CO_2 = s(A,B,C,D,E,G)
ds_CO_2 = s_CO_2 - s_29815_CO
g_CO_2 = G1(h_CO_2,t,s_CO_2)
dg_CO_2 = dG(dH_CO_2,t,ds_CO_2)

dH_CO = np.array(list(dH_CO_1)+list(dH_CO_2))
h_CO = np.array(list(h_CO_1)+list(h_CO_2))
s_CO = np.array(list(s_CO_1)+list(s_CO_2))
ds_CO = np.array(list(ds_CO_1)+list(ds_CO_2))
g_CO = np.array(list(g_CO_1)+list(g_CO_2))
dg_CO = np.array(list(dg_CO_1)+list(dg_CO_2))



print("\n                                           CARBON MONOXIDE (CO)")
print("\n      Temperature |      dH       |       H       |       dS      |       S       |       G       |       *dG*      ")
print("     ---------------------------------------------------------------------------------------------------------------")

for i in np.arange(26):
    print("\t",np.linspace(500,3000,26)[i],"  |",dH_CO[i],"|",h_CO[i],"|",ds_CO[i],"|",s_CO[i],"|",g_CO[i],"|",dg_CO[i])

#       H2


T = np.linspace(500,1000,6) # degrees K
t = T/1000
       
# 298-1000 K valid temperature range
A =  33.066178
B = -11.363417
C =  11.432816
D = -2.772874
E = -0.158558
F = -9.980797
G =  172.707974
H =  0.0

Hf_29815_H2 = 0.0 #this is Hf.
s_29815_H2 = 130.68 #P = 1bar
g_29815_H2 = G1(Hf_29815_H2,298.15/1000,s_29815_H2)

dH_H2_1 = dH(A,B,C,D,E,F,H)
h_H2_1 = dH_H2_1 + Hf_29815_H2
s_H2_1 = s(A,B,C,D,E,G)
ds_H2_1 = s_H2_1 - s_29815_H2
g_H2_1 = G1(h_H2_1,t,s_H2_1)
dg_H2_1 = dG(dH_H2_1,t,ds_H2_1)
dg2_H2_1 = dG2(g_H2_1,g_29815_H2)

T = np.linspace(1100,2500,15) # degrees K
t = T/1000

# 1000-2500 K valid temperature range
A =  18.563083
B =  12.257357
C = -2.859786
D =  0.268238
E =  1.977990
F = -1.147438
G =  156.288133
H =  0.0

dH_H2_2 = dH(A,B,C,D,E,F,H)
h_H2_2 = dH_H2_2 + Hf_29815_H2
s_H2_2 = s(A,B,C,D,E,G)
ds_H2_2 = s_H2_2 - s_29815_H2
g_H2_2 = G1(h_H2_2,t,s_H2_2)
dg_H2_2 = dG(dH_H2_2,t,ds_H2_2)
dg2_H2_2 = dG2(g_H2_2,g_29815_H2)
    
T = np.linspace(2600,3000,5) # degrees K
t = T/1000

# 2500-6000 K valid temperature range
A =  43.413560
B = -4.293079
C =  1.272428
D = -0.096876
E = -20.533862
F = -38.515158
G =  162.081354
H =  0.0

dH_H2_3 = dH(A,B,C,D,E,F,H)
h_H2_3 = dH_H2_3 + Hf_29815_H2
s_H2_3 = s(A,B,C,D,E,G)
ds_H2_3 = s_H2_3 - s_29815_H2
g_H2_3 = G1(h_H2_3,t,s_H2_3)
dg_H2_3 = dG(dH_H2_3,t,ds_H2_3)
dg2_H2_3 = dG2(g_H2_3,g_29815_H2)

dH_H2 = np.array(list(dH_H2_1)+list(dH_H2_2)+list(dH_H2_3))
h_H2 = np.array(list(h_H2_1)+list(h_H2_2)+list(h_H2_3))
s_H2 = np.array(list(s_H2_1)+list(s_H2_2)+list(s_H2_3))
ds_H2 = np.array(list(ds_H2_1)+list(ds_H2_2)+list(ds_H2_3))
g_H2 = np.array(list(g_H2_1)+list(g_H2_2)+list(g_H2_3))
dg_H2 = np.array(list(dg_H2_1)+list(dg_H2_2)+list(dg_H2_3))
dg2_H2 = np.array(list(dg2_H2_1)+list(dg2_H2_2)+list(dg2_H2_3))



print("\n                                            HYDROGEN (H2)")
print("\n      Temperature |      dH       |       H       |       dS      |       S       |       G       |       *dG*       |       *dG*      ")
print("     ---------------------------------------------------------------------------------------------------------------")
for i in np.arange(26):
    print("\t",np.linspace(500,3000,26)[i],"  |",dH_H2[i],"|",h_H2[i],"|",ds_H2[i],"|",s_H2[i],"|",g_H2[i],"|",dg_H2[i],"|",dg2_H2[i])

#       CO2

T = np.linspace(500,1200,8) # degrees K
t = T/1000
       
# 298-1200 K valid temperature range
A =   24.99735
B =   55.18696
C =  -33.69137
D =   7.948387
E =  -0.136638
F =  -403.6075
G =   228.2431
H =  -393.5224

Hf_29815_CO2 = -393.52 #this is Hf.
s_29815_CO2 = 213.79 #P = 1bar

dH_CO2_1 = dH(A,B,C,D,E,F,H)
h_CO2_1 = dH_CO2_1 + Hf_29815_CO2
s_CO2_1 = s(A,B,C,D,E,G)
ds_CO2_1 = s_CO2_1 - s_29815_CO2
g_CO2_1 = G1(h_CO2_1,t,s_CO2_1)
dg_CO2_1 = dG(dH_CO2_1,t,ds_CO2_1)

T = np.linspace(1300,3000,18) # degrees K
t = T/1000

# 1200-6000 K valid temperature range
A =   58.16639
B =   2.720074
C =  -0.492289
D =   0.038844
E =  -6.447293
F =  -425.9186
G =   263.6125
H =  -393.5224

dH_CO2_2 = dH(A,B,C,D,E,F,H)
h_CO2_2 = dH_CO2_2 + Hf_29815_CO2
s_CO2_2 = s(A,B,C,D,E,G)
ds_CO2_2 = s_CO2_2 - s_29815_CO2
g_CO2_2 = G1(h_CO2_2,t,s_CO2_2)
dg_CO2_2 = dG(dH_CO2_2,t,ds_CO2_2)

dH_CO2 = np.array(list(dH_CO2_1)+list(dH_CO2_2))
h_CO2 = np.array(list(h_CO2_1)+list(h_CO2_2))
s_CO2 = np.array(list(s_CO2_1)+list(s_CO2_2))
ds_CO2 = np.array(list(ds_CO2_1)+list(ds_CO2_2))
g_CO2 = np.array(list(g_CO2_1)+list(g_CO2_2))
dg_CO2 = np.array(list(dg_CO2_1)+list(dg_CO2_2))


print("\n                                           CARBON DIOXIDE (CO2)")
print("\n      Temperature |      dH       |       H       |       dS      |       S       |       G       |       *dG*      ")
print("     ---------------------------------------------------------------------------------------------------------------")

for i in np.arange(26):
    print("\t",np.linspace(500,3000,26)[i],"  |",dH_CO2[i],"|",h_CO2[i],"|",ds_CO2[i],"|",s_CO2[i],"|",g_CO2[i],"|",dg_CO2[i])

#       C2H2

T = np.linspace(500,1100,7) # degrees K
t = T/1000
       
# 298-1100 K valid temperature range
A =   40.68697
B =   40.73279
C =  -16.17840
D =   3.669741
E =  -0.658411
F =   210.7067
G =   235.0052
H =   226.7314

Hf_29815_C2H2 = 226.73 #this is Hf.
s_29815_C2H2 = 200.93 #P = 1bar

dH_C2H2_1 = dH(A,B,C,D,E,F,H)
h_C2H2_1 = dH_C2H2_1 + Hf_29815_C2H2
s_C2H2_1 = s(A,B,C,D,E,G)
ds_C2H2_1 = s_C2H2_1 - s_29815_C2H2
g_C2H2_1 = G1(h_C2H2_1,t,s_C2H2_1)
dg_C2H2_1 = dG(dH_C2H2_1,t,ds_C2H2_1)

T = np.linspace(1200,3000,19) # degrees K
t = T/1000

# 1100-6000 K valid temperature range
A =   67.47244
B =   11.75110
C =  -2.021470
D =   0.136195
E =  -9.806418
F =   185.4550
G =   253.5337
H =   226.7314

dH_C2H2_2 = dH(A,B,C,D,E,F,H)
h_C2H2_2 = dH_C2H2_2 + Hf_29815_C2H2
s_C2H2_2 = s(A,B,C,D,E,G)
ds_C2H2_2 = s_C2H2_2 - s_29815_C2H2
g_C2H2_2 = G1(h_C2H2_2,t,s_C2H2_2)
dg_C2H2_2 = dG(dH_C2H2_2,t,ds_C2H2_2)

dH_C2H2 = np.array(list(dH_C2H2_1)+list(dH_C2H2_2))
h_C2H2 = np.array(list(h_C2H2_1)+list(h_C2H2_2))
s_C2H2 = np.array(list(s_C2H2_1)+list(s_C2H2_2))
ds_C2H2 = np.array(list(ds_C2H2_1)+list(ds_C2H2_2))
g_C2H2 = np.array(list(g_C2H2_1)+list(g_C2H2_2))
dg_C2H2 = np.array(list(dg_C2H2_1)+list(dg_C2H2_2))


print("\n                                           ACETYLENE (C2H2)")
print("\n      Temperature |      dH       |       H       |       dS      |       S       |       G       |       *dG*      ")
print("     ---------------------------------------------------------------------------------------------------------------")

for i in np.arange(26):
    print("\t",np.linspace(500,3000,26)[i],"  |",dH_C2H2[i],"|",h_C2H2[i],"|",ds_C2H2[i],"|",s_C2H2[i],"|",g_C2H2[i],"|",dg_C2H2[i])

# NO DA LOS VALORES DEL PAPER (molar gibbs free energy), LLORAR Y NEXT.

#reaction 1: CH4 + H2O <-> CO + 3 H2

Hrxn1_29815 = Hrxn_29815(Hf_29815_CH4,Hf_29815_H2O,Hf_29815_CO,3*Hf_29815_H2)
Srxn1_29815 = Srxn_29815(s_29815_CH4,s_29815_H2O,s_29815_CO,3*s_29815_H2)
Grxn1_29815 = Grxn_29815(Hrxn1_29815,Srxn1_29815)

print('deltaH_29815_1 = {0:.3f}'.format(Hrxn1_29815)) 
print('deltaS_29815_1 = {0:.3f}'.format(Srxn1_29815))
print('deltaG_29815_1 = {0:.3f}'.format(Grxn1_29815))

T = np.linspace(500,3000,26)

Hrxn1 = Hrxn1_29815 + dH_CO + 3*dH_H2 - dH_CH4 - dH_H2O
Grxn1 = Hrxn1 - T*(s_CO + 3*s_H2 - s_CH4 - s_H2O)/1000

#print(Grxn1)

plt.figure()
plt.plot(T,Grxn1, label='$\Delta G_{rxn1}$')
plt.plot(T,Hrxn1, label='$\Delta H_{rxn1}$')
plt.xlabel('Temperature (K)')
plt.ylabel('(kJ/mol)')
plt.legend( loc='best')
plt.savefig("wgs-nist-1.png")
#plt.show()

# Equilibrium constant calculation (K')

R = 8.314*10**(-3) # kJ/(mol K)
P_0 = 1 #bar
P_1 = 1 #bar
P_2 = 0.01 #bar
K1_1bar = (P_0/P_1)**2 * np.exp(-Grxn1/(R*T))
K1_10mbar = (P_0/P_2)**2 * np.exp(-Grxn1/(R*T))

plt.figure()
plt.plot(T,K1_1bar,label='$K1_{1bar}$')
plt.plot(T,K1_10mbar,label='$K1_{10mbar}$')
plt.xlim([500, 3000])
plt.yscale('log')
plt.xlabel('Temperature (K)')
plt.ylabel('Equilibrium constant rxn1')
plt.legend( loc='best')
plt.savefig('wgs-nist-1-K.png')
#plt.show()

#reaction 2: CO2 + H2 <-> CO + H2O

Hrxn2_29815 = Hrxn_29815(Hf_29815_CO2,Hf_29815_H2,Hf_29815_CO,Hf_29815_H2O)
Srxn2_29815 = Srxn_29815(s_29815_CO2,s_29815_H2,s_29815_CO,s_29815_H2O)
Grxn2_29815 = Grxn_29815(Hrxn2_29815,Srxn2_29815)

print('deltaH_29815_2 = {0:.3f}'.format(Hrxn2_29815)) 
print('deltaS_29815_2 = {0:.3f}'.format(Srxn2_29815))
print('deltaG_29815_2 = {0:.3f}'.format(Grxn2_29815))

T = np.linspace(500,3000,26)

Hrxn2 = Hrxn2_29815 + dH_CO + dH_H2O - dH_CO2 - dH_H2
Grxn2 = Hrxn2 - T*(s_CO + s_H2O - s_CO2 - s_H2)/1000

#print(Grxn2)

plt.figure()
plt.plot(T,Grxn2, label='$\Delta G_{rxn2}$')
plt.plot(T,Hrxn2, label='$\Delta H_{rxn2}$')
plt.xlabel('Temperature (K)')
plt.ylabel('(kJ/mol)')
plt.legend( loc='best')
plt.savefig("wgs-nist-2.png")
#plt.show()

# Equilibrium constant calculation (K')

K2_1bar = (P_0/P_1)**2 * np.exp(-Grxn2/(R*T))
K2_10mbar = (P_0/P_2)**2 * np.exp(-Grxn2/(R*T))
K2_sin_presion = np.exp(-Grxn2/(R*T))

plt.figure()
plt.plot(T,K2_1bar,label='$K´_{2}(1bar)$')
plt.plot(T,K2_10mbar,label='$K´_{2}(10mbar)$',ls='--')
plt.plot(T,K2_sin_presion,label='$K´_{2}$',ls='--')
plt.xlim([500, 3000])
plt.yscale('log')
plt.xlabel('Temperature (K)')
plt.ylabel('Equilibrium constant rxn2')
plt.legend( loc='best')
plt.savefig('wgs-nist-2-K.png')
#plt.show()

plt.figure()
plt.plot(T,1/K2_sin_presion,label='$1/K´_{2}$',color='r')
plt.xlim([500, 3000])
plt.yscale('log')
plt.xlabel('Temperature (K)')
plt.ylabel('$1/K´_{2}$')
plt.legend( loc='best')
plt.savefig('wgs-nist-fig1.png')
#plt.show()

#reaction 3: 2 CH4 <-> C2H2 + 3 H2

Hrxn3_29815 = Hrxn_29815(2*Hf_29815_CH4,0.0,Hf_29815_C2H2,3*Hf_29815_H2)
Srxn3_29815 = Srxn_29815(2*s_29815_CH4,0.0,s_29815_C2H2,3*s_29815_H2)
Grxn3_29815 = Grxn_29815(Hrxn3_29815,Srxn3_29815)

print('deltaH_29815_3 = {0:.3f}'.format(Hrxn3_29815)) 
print('deltaS_29815_3 = {0:.3f}'.format(Srxn3_29815))
print('deltaG_29815_3 = {0:.3f}'.format(Grxn3_29815))

T = np.linspace(500,3000,26)

Hrxn3 = Hrxn3_29815 + dH_C2H2 + 3*dH_H2 - 2*dH_CH4
Grxn3 = Hrxn3 - T*(s_C2H2 + 3*s_H2 - 2*s_CH4)/1000

#print(Grxn3)

plt.figure()
plt.plot(T,Grxn3, label='$\Delta G_{rxn3}$')
plt.plot(T,Hrxn3, label='$\Delta H_{rxn3}$')
plt.xlabel('Temperature (K)')
plt.ylabel('(kJ/mol)')
plt.legend( loc='best')
plt.savefig("wgs-nist-3.png")
#plt.show()

# Equilibrium constant calculation (K')

K3_1bar = (P_0/P_1)**2 * np.exp(-Grxn3/(R*T))
K3_10mbar = (P_0/P_2)**2 * np.exp(-Grxn3/(R*T))

plt.figure()
plt.plot(T,K3_1bar,label='$K3_{1bar}$')
plt.plot(T,K3_10mbar,label='$K3_{10mbar}$')
plt.xlim([500, 3000])
plt.yscale('log')
plt.xlabel('Temperature (K)')
plt.ylabel('Equilibrium constant rxn3')
plt.legend( loc='best')
plt.savefig('wgs-nist-3-K.png')
#plt.show()

plt.figure()
plt.plot(T,K1_10mbar,label='$K´_{1}(10mbar)$',color='orange',ls='-')
plt.plot(T,K1_1bar,label='$K´_{1}(1bar)$',color='orange',ls='-',linewidth=3)
plt.plot(T,K2_sin_presion,label='$K´_{2}$',color='r',ls='-.')
plt.plot(T,K3_10mbar,label='$K´_{3}(10mbar)$',color='k',ls='--')
plt.plot(T,K3_1bar,label='$K´_{3}(1bar)$',color='k',ls='--',linewidth=3)
plt.xlim([500, 3000])
plt.yscale('log')
plt.xlabel('Temperature (K)')
plt.ylabel('Normalised Equilibrium Constants')
plt.legend( loc='best')
plt.savefig('wgs-nist-fig2.png')
#plt.show()

# LAS 3 REACCIONES FUNCIONAN!

ratio_CO_solar = 0.46 # Por definición

# Si quieres generar abundancias con C/O = X, por ejemplo, simplemente aumenta el carbono por un factor de (X/[C/O_solar]), de manera que la nueva razón de carbono a oxígeno será X.

# ABUNDANCIAS ELEMENTALES SOLARES

n_H = 1.0
n_He = 0.0968277856261
n_Ti = 9.54992586021*10**(-8)
n_V = 1.09647819614*10**(-8)
n_O = 0.000602559586074
n_C = 0.000275422870334
n_N = 8.12830516164*10**(-5)
n_S = 1.62181009736*10**(-5)
n_Na = 2.23872113857*10**(-6)
n_K = 1.44543977075*10**(-7)
n_Fe = 3.2359365693*10**(-5)

rate_CO_solar = n_C/n_O # This number is C/O = 0.46

# ABUNDANCES WITH SOLAR C/O=0.46

n_CH4_1bar, n_CH4_10mbar, n_C2H2_1bar, n_C2H2_10mbar, n_H2O_1bar, n_H2O_10mbar, n_CO_1bar, n_CO_10mbar, n_CO2_1bar, n_CO2_10mbar = abundances_norm_H(n_C)

# IT WORKS! NOW GRAPHICS!

graphics_CO("1bar", "fig3-a", '{0:.3f}'.format(rate_CO_solar))

# YASS!!

# NOW WE MODIFY C/O AND GENERATES MORE GRAPHICS

n_O = 5 * 10**(-4)

# C/O = 0.1    

CO_wanted = 0.1
n_C_wanted = (CO_wanted/rate_CO_solar) * n_C

n_CH4_1bar, n_CH4_10mbar, n_C2H2_1bar, n_C2H2_10mbar, n_H2O_1bar, n_H2O_10mbar, n_CO_1bar, n_CO_10mbar, n_CO2_1bar, n_CO2_10mbar = abundances_norm_H(n_C_wanted)
graphics_CO("", "fig4-a",'{0:.3f}'.format(CO_wanted))

# C/O = 1    


CO_wanted = 1.0
n_C_wanted = (CO_wanted/rate_CO_solar) * n_C

n_CH4_1bar, n_CH4_10mbar, n_C2H2_1bar, n_C2H2_10mbar, n_H2O_1bar, n_H2O_10mbar, n_CO_1bar, n_CO_10mbar, n_CO2_1bar, n_CO2_10mbar = abundances_norm_H(n_C_wanted)
graphics_CO("", "fig4-b", '{0:.3f}'.format(CO_wanted))


# C/O = 10          

CO_wanted = 10.0
n_C_wanted = (CO_wanted/rate_CO_solar) * n_C

n_CH4_1bar, n_CH4_10mbar, n_C2H2_1bar, n_C2H2_10mbar, n_H2O_1bar, n_H2O_10mbar, n_CO_1bar, n_CO_10mbar, n_CO2_1bar, n_CO2_10mbar = abundances_norm_H(n_C_wanted)
graphics_CO("", "fig4-c", '{0:.3f}'.format(CO_wanted))


# C/O = 100          

CO_wanted = 100.0
n_C_wanted = (CO_wanted/rate_CO_solar) * n_C  

n_CH4_1bar, n_CH4_10mbar, n_C2H2_1bar, n_C2H2_10mbar, n_H2O_1bar, n_H2O_10mbar, n_CO_1bar, n_CO_10mbar, n_CO2_1bar, n_CO2_10mbar = abundances_norm_H(n_C_wanted)
graphics_CO("", "fig4-d", '{0:.3f}'.format(CO_wanted))
 

# Algo no está funcinando bien con el H2O, CO y CO2. Con 1 bar desde 1100K  con 10mbar con 1500K







