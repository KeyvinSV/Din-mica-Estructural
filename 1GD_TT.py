import Funciones as fun
import matplotlib.pyplot as plt
import numpy as np

# INPUT DE PROPIEDADES DINÁMICAS
T=0.5                        #Periodo    [s]
β1=5                         #Amortiguamiento critico [%]

# LECTURA DE REGISTRO SÍSMICO
ug=np.genfromtxt("./Sismo.txt")  #Registro sísmico [cm/s2]
Δt=0.002                         #Intervalo de tiempo [s]
N=len(ug)
t=[i*Δt for i in range(N)]

# CÁLCULO DE LAS RESPUESTAS SÍSMICAS
#x,xp,xpp=fun.DCentral(-ug,Tn=T,ζ=β1,Δt=Δt)
x,xp,xpp=fun.Newmark(-ug,Tn=T,ζ=β1,Δt=Δt,β=1/4,γ=1/2)

# GUARDADO DE RESULTADOS EN ARCHIVO TXT
fp=open("1GDL.txt","w")
fp.write("t [s]\tU''g [cm/s2] \tX'' [cm/s2] \tX' [cm/s]\tX [cm]")
for i in range(N-1):
    fp.write("\n%f\t%f\t%f\t%f\t%f" % (i*Δt,ug[i],xpp[i],xp[i],x[i]))
fp.close()










