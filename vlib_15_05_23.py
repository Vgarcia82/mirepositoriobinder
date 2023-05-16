from numpy import sin, cos
import numpy as np



#state[0]=angulo del brazo
#state[1]=velocidad del brazo
#state[2]=angulo del pendulo
#state[3]=velocidad del pendulo

def derivs(state,t,M2,M1,L1,L2,J1,J2,G,l1,l2,tau1,tau2):# calculo coeficientes y 

    Gr=9.81
    theta1 = M2*L1*L1+J1
    theta2 = M2*l2*l2
    theta3 = L1*l2*M2
    theta4 = M2*l2*l2+J2
    theta5 = l2*M2*Gr
    tehta6 = 0
    theta7 = 0
    theta8 = 0
    theta9 = 0
    # Matriz masas
    M11 = theta1+theta2*sin(state[2])*sin(state[2])   
    M12 = theta3*cos(state[2]) 
    M22 = theta4 
    MatM = np.matrix([[M11, M12], [M12 ,M22]])
    #-------------------------------Coriolis
    
    C11 = 0.5*theta2*state[3]*sin(2*state[2])
    C12 = -theta3*state[3]*sin(state[2])+0.5*theta2*state[1]*sin(2*state[2])
    C21 = -0.5*theta2*state[1]*sin(2*state[2])
    
    MatC = np.matrix([[C11, C12], [C21, 0]])
    #------------------------------------------
    

    g22 = -theta5*sin(state[2]) #potencial 
    #-----------------------------------------------------------------
    MatA = np.matrix([[tau1-C11*state[1]-C12*state[3]], [tau2-C21*state[1]-g22]])
   
   
    #---------------------------------------------
    K=np.matmul(np.linalg.inv(MatM),MatA)
    #---------------------------------------------------
    
    dqdt = np.zeros_like(state)
    dqdt[0] = state[1]
   
    dqdt[1] = K[0,0]#
    dqdt[2] = state[3]
   
    dqdt[3] = K[1,0]
  
    return dqdt


def fqq(state,M2,M1,L1,L2,J1,J2,G,l1,l2,kw,ke,ktheta,kdelta):# calculo lypanov funcion control y tau

    Gr=9.81
    theta1 = M2*L1*L1+J1
    theta2 = M2*l2*l2
    theta3 = L1*l2*M2
    theta4 = M2*l2*l2+J2
    theta5 = l2*M2*Gr
    tehta6 = 0
    theta7 = 0
    theta8 = 0
    theta9 = 0
    #----------
    M11 = theta1+theta2*sin(state[2])*sin(state[2])    
    M12 = theta3*cos(state[2]) 
    M22 = theta4  
    MatM = np.matrix([[M11, M12], [M12 ,M22]])
    detM=np.linalg.det(MatM)
    #----------------- def F página 81
   
    f1=-theta4*theta2*sin(2*state[2])*state[1]*state[3]
    f2=-0.5*theta3*theta2*cos(state[2])*sin(2*state[2])*state[1]*state[1]
    f3=theta4*theta3*sin(state[2])*state[3]*state[3]
    f4=-theta3*theta5*cos(state[2])*sin(state[2])
    #print("state",state,"T",t)
    F=f1+f2+f3+f4
    #--------------------------------------------- lupyanov            
    Energia=theta5*(cos(state[2])-1)+0.5*(M11*state[1]*state[1]+2*M12*state[1]*state[3]+M22*state[3]*state[3])
    print("Energia= ",Energia)
    
    V=ke*0.5*Energia*Energia+kw*0.5*state[1]*state[1]+ktheta*0.5*state[0]*state[0]
    print("V",V)
    print("V(0)",2*ke*theta5*theta5)
    #------------------------------------------------------------------------
    tau1=-kw*F-detM*(kdelta*state[1]+ktheta*state[0])
    tau2=detM*ke*Energia+kw*theta4
    tauf=tau1/tau2
    #------------------------ley de control
    r=2*theta5*(theta1+theta2)
    if kw/ke>r:
        print("control low= Ok")
      
    else:
         print("control low= Falsch")     
    print("kw/ke=",kw/ke," r=",r)
    return tauf,Energia
def angulogrados(input):
    for vy2 in range(len(input)):
   
        if input[vy2]>360 :
            vdp=np.trunc(input[vy2]/360)

            if vdp!=0:
                input[vy2]=input[vy2]-vdp*360
        elif input[vy2]<-360:        
            vdp=np.trunc(input[vy2]/360)

            if vdp!=0:
                input[vy2]=360+(input[vy2]-vdp*360)
        elif input[vy2]>-360 and input[vy2]<0 :
                input[vy2]=360+(input[vy2])
  
    return input