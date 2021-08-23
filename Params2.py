'''Parameters of system - Inputs'''
Nmet = 4 # number of metabolites in system
constext = True # do constant exterior?
Rc = 1.0e-5 #cm
Kab = 0.5e3 #microM
Kbc = 15.0e3 #microM
Kbx = 10.0e3 #microM
kcatab = 300.0  #reactions/s
kcatbc = 55.0 #reactions/s
kcatbx = 55.0 #reactions/s
Nab = 1500.0 #Pdu a-b active sites
Nbc = 2500.0 # Pdu b-c active sites
Nbx = 2500.0 # Pdu b-x active sites
DA = 1.0e-5 #cm^2/s
DB = 1.0e-5 #cm^2/s
DC = 1.0e-5 #cm^2/s
DX = 1.0e-5 #cm^2/s
kcA = 1.e-4 #cm/s
kcB = kcA #cm/s
kcC = kcA #cm/s
kcX = kcA #cm/s
Acyt = 30e3 #microM cytosolic concentration
Bcyt = 0 #microM cytosolic concentration
Ccyt = 0 #product cytocolic concentration
Xcyt = 0 #product cytocolic concentration
tstop = 10#30000 #s final time point for integration 
numgrid = 100.0 #number of grid points the radius is split into
MCPmil = 1.e12

Amultiplier = 1.
Bmultiplier = 1.

'''Parameter combinations - do not touch'''
import numpy as np
VMCP = np.divide(4.,3.)*np.pi*np.power(Rc,3)
Vratio = 20000#np.divide(1-VMCP*MCPmil,VMCP*MCPmil) #ratio of MCP volume to free bath volume per MCP 
Na = 6.022e23 #avogadro's number - constant
MCPMolar = MCPmil*1000/Na
Vab = Amultiplier*np.divide(kcatab*Nab*1.0e9,(np.divide(4,3.0)*np.pi*(Rc**3)*Na))
Vbc = Bmultiplier*np.divide(kcatbc*Nbc*1.0e9,(np.divide(4,3.0)*np.pi*(Rc**3)*Na))
Vbx = Bmultiplier*np.divide(kcatbx*Nbx*1.0e9,(np.divide(4,3.0)*np.pi*(Rc**3)*Na))
kappa = Kab/Kbc
gamma = Vbc/Vab
beta = Vbx/Vab
eta = Kbx/Kbc
xiA = (Kbc*DA)/(Vab*(Rc**2))
xiB = (Kbc*DB)/(Vab*(Rc**2))
xiC = (Kbc*DC)/(Vab*(Rc**2))
xiX = (Kbc*DX)/(Vab*(Rc**2))
tau = Kbc/Vab
chiA = (Kbc*kcA)/(Vab*Rc)
chiB = (Kbc*kcB)/(Vab*Rc)
chiC = (Kbc*kcC)/(Vab*Rc)
chiX = (Kbc*kcX)/(Vab*Rc)
acyt = Acyt/Kab
bcyt = Bcyt/Kbc
ccyt = Ccyt/Kbc
xcyt = Xcyt/Kbc
Aext = 15.0e4 / Kab
Bext = 12.5e4 / Kbc
Cext = 12.5e4 / Kbc
Xext = 13.0e4 / Kbc
dm = 1/float(numgrid)

tfinal = int(np.divide(tstop,tau))


# nM = 3*np.int(np.reciprocal(dm))+3
# initials = np.zeros(nM)

#p = np.array([xiP, xiA, xiC, kappa, gamma, chiP, chiA, chiC, pcyt, acyt, Vratio, dm])
p = np.array([xiA,xiB,xiC,xiX,chiA,chiB,chiC,chiX,Aext,Bext,Cext,Xext,kappa,gamma,beta,eta,Vratio,Nmet,dm])