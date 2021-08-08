import numpy as np
# p = np.array([xiA,xiB,xiC,xiX,chiA,chiB,chiC,chiX,Aext,Bext,Cext,Dext,kappa,gamma,beta,eta,Vratio,Nmet,dm])
#               0   1   2   3   4    5    6    7    8    9    10   11   -7    -6    -5   -4  -3     -2   -1
def deriv2(allvar,time, params,constext=False):
#     constext = False   #######
    kappa = params[-7]
    gamma = params[-6]
    beta = params[-5]
    eta = params[-4]
    Vratio = params[-3]
    Nmet = np.int(params[-2])
    dm = params[-1]
    nM =int(np.floor(np.divide(1,dm)))
    d = np.zeros(Nmet*(nM+1))

    if constext:
        variables = allvar
        ext = params[8:8+Nmet]
        allder = np.zeros(Nmet*(nM+1))
    else:
        variables = allvar[:-Nmet]
        ext = allvar[-Nmet:]
        allder = np.zeros(Nmet*(nM+2))
        
    for i in range(1,nM):
        d[Nmet*i + 0] = np.divide(9*(((i+0.5)**(4/3.0))*(variables[Nmet*(i+1)+0]-variables[Nmet*i+0]) -((i-0.5)**(4/3.0))*(variables[Nmet*i+0] - variables[Nmet*(i-1)+0]))*params[0],dm**(2/3.0)) - np.divide(variables[Nmet*i],kappa*(1+variables[Nmet*i]))
        d[Nmet*i + 1] = np.divide(9*(((i+0.5)**(4/3.0))*(variables[Nmet*(i+1)+1]-variables[Nmet*i+1]) -((i-0.5)**(4/3.0))*(variables[Nmet*i+1] - variables[Nmet*(i-1)+1]))*params[1],dm**(2/3.0)) + np.divide(variables[Nmet*i],(1+variables[Nmet*i])) - np.divide(gamma*variables[Nmet*i+1],1+variables[Nmet*i+1])
        d[Nmet*i + 2] = np.divide(9*(((i+0.5)**(4/3.0))*(variables[Nmet*(i+1)+2]-variables[Nmet*i+2]) -((i-0.5)**(4/3.0))*(variables[Nmet*i+2] - variables[Nmet*(i-1)+2]))*params[2],dm**(2/3.0)) + np.divide(gamma*variables[Nmet*i+1],(1+variables[Nmet*i+1]))
        if Nmet==4:
            d[Nmet*i + 1] -= np.divide(beta*variables[Nmet*i+1],eta+variables[Nmet*i+1])
            d[Nmet*i + 3] = np.divide(9*(((i+0.5)**(4/3.0))*(variables[Nmet*(i+1)+3]-variables[Nmet*i+3]) -((i-0.5)**(4/3.0))*(variables[Nmet*i+3] - variables[Nmet*(i-1)+3]))*params[3],dm**(2/3.0)) + np.divide(beta*variables[Nmet*i+1], eta+variables[Nmet*i+1])

    d[0] = np.divide(18*(0.5**(4/3.0))*(variables[Nmet+0] - variables[0])*params[0],dm**(2/3.0)) - np.divide(variables[0],kappa*(1+variables[0]))
    d[1] = np.divide(18*(0.5**(4/3.0))*(variables[Nmet+1] - variables[1])*params[1], dm**(2/3.0)) + np.divide(variables[0],1+variables[0]) - np.divide(gamma*variables[1],1+variables[1])
    d[2] = np.divide(18*(0.5**(4/3.0))*(variables[Nmet+2] - variables[2])*params[2], dm**(2/3.0)) + np.divide(gamma*variables[1],1+variables[1])
    if Nmet == 4:
        d[1] -= np.divide(beta*variables[1],eta+variables[1])
        d[3] = np.divide(18*(0.5**(4/3.0))*(variables[Nmet+3] - variables[3])*params[3], dm**(2/3.0)) + np.divide(beta*variables[1],eta+variables[1])

    d[-Nmet+0] = np.divide(3*((nM+0.5)**(2/3.0))*(ext[0] - variables[-Nmet+0])*params[4],dm**(1/3.0)) - np.divide(9*((nM-0.5)**(4/3.0))*(variables[-Nmet+0] - variables[-2*Nmet+0])*params[0],dm**(2/3.0))
    d[-Nmet+1] = np.divide(3*((nM+0.5)**(2/3.0))*(ext[1] - variables[-Nmet+1])*params[5],dm**(1/3.0)) - np.divide(9*((nM-0.5)**(4/3.0))*(variables[-Nmet+1] - variables[-2*Nmet+1])*params[1],dm**(2/3.0))
    d[-Nmet+2] = np.divide(3*((nM+0.5)**(2/3.0))*(ext[2] - variables[-Nmet+2])*params[6],dm**(1/3.0)) - np.divide(9*((nM-0.5)**(4/3.0))*(variables[-Nmet+2] - variables[-2*Nmet+2])*params[2],dm**(2/3.0))
    if Nmet == 4:
        d[-Nmet+3] = np.divide(3*((nM+0.5)**(2/3.0))*(ext[3] - variables[-Nmet+3])*params[7],dm**(1/3.0)) - np.divide(9*((nM-0.5)**(4/3.0))*(variables[-Nmet+3] - variables[-2*Nmet+3])*params[3],dm**(2/3.0))
    
    if constext:
        allder = d
    else:
        dext = np.zeros(Nmet)

        dext[0] = -np.divide(3*((nM+0.5)**(2/3.0))*(ext[0] - variables[-Nmet+0])*params[4],dm**(1/3.0))*np.divide(dm,Vratio)
        dext[1] = -np.divide(3*((nM+0.5)**(2/3.0))*(ext[1] - variables[-Nmet+1])*params[5],dm**(1/3.0))*np.divide(dm,Vratio)
        dext[2] = -np.divide(3*((nM+0.5)**(2/3.0))*(ext[2] - variables[-Nmet+2])*params[6],dm**(1/3.0))*np.divide(dm,Vratio)
        if Nmet == 4:
            dext[3] = -np.divide(3*((nM+0.5)**(2/3.0))*(ext[3] - variables[-Nmet+3])*params[7],dm**(1/3.0))*np.divide(dm,Vratio)

        allder[:-Nmet] = d
        allder[-Nmet:] = dext
    return allder
