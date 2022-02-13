import numpy as np

def Newmark(p,  u0 = 0 , up0 = 0, Tn =0.1, ζ =5, Δt =0.01, β = 1/6 , γ = 1/2 ):
        """
        γ: parametro de presición, generalmente 1/2
        β: razón de la variacion de la aceleración, generalmente entre 1/4 y 1/6 """

        # Parámetros Dinámicos
        m=1.0
        k=m*(2*np.pi/Tn)**2
        Wn=2*np.pi/Tn
        c=2*ζ/100*m*Wn

        # Calcular a1, a2, a3 y kt
        a1=m/(β*Δt**2)+γ*c/(β*Δt)
        a2=m/(β*Δt)+(γ/β-1)*c
        a3=(1/(2*β)-1)*m+Δt*(γ/(2*β)-1)*c
        kt=k+a1

        # Condiciones Iniciales
        u=np.zeros(len(p))
        u[0]=u0
        up=np.zeros(len(p))
        up[0]=up0
        upp=np.zeros(len(p))
        upp[0]=(p[0]-c*up0-k*u0)/m

        # Proceso Iterativo
        for i in range(len(p)-1):
            pti_1=p[i+1]+a1*u[i]+a2*up[i]+a3*upp[i]
            u[i+1]=pti_1/kt
            up[i+1]=(γ/(β*Δt))*(u[i+1]-u[i])+(1-γ/β)*up[i]+Δt*(1-γ/(2*β))*upp[i]
            upp[i+1]=(u[i+1]-u[i])/(β*Δt**2)-up[i]/(β*Δt)-(1/(2*β)-1)*upp[i]

        return u,up,upp


def DCentral(p,  u0 = 0 , up0 = 0, Tn =0.1, ζ =5, Δt =0.01):
        """
        El método es convergente siempre que el intervalo de tiempo se menor al
        valor crítico Δt<Δtcr=T/4   """

        # Parámetros Dinámicos
        m=1.0
        k=m*(2*np.pi/Tn)**2
        Wn=2*np.pi/Tn
        c=2*ζ/100*m*Wn

        # Calcular a0, a1, a2 y a3
        a0=1/Δt**2;
        a1=1/(2*Δt);
        a2=2*a0;
        a3=1/a2

        # Condiciones Iniciales
        u=np.zeros(len(p));
        u[0]=u0
        up=np.zeros(len(p))
        up[0]=up0
        upp=np.zeros(len(p))
        upp[0]=(p[0]-c*up0-k*u0)/m

        Xa=u[0]-Δt*up[0]+a3*upp[0]
        X=u[0]

        # Proceso Iterativo
        for i in range(len(p)-1):
            Xb=1/(m*a0+c*a1)*(p[i]-(k-a2*m)*X-(a0*m-a1*c)*Xa)
            if i!=0:
                upp[i]=a0*(Xb-2*X+Xa)
                up[i]=a1*(Xb-Xa)
                u[i]=X
            Xa=X; X=Xb

        return u,up,upp