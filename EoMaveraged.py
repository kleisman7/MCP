import numpy as np

def derivAve(allvar,time, params):
    nM =int(np.floor(np.divide(1,params[-1])))
    d = np.zeros(6)

    variables = allvar[:-3]
    ext = allvar[-3:]
    allder = np.zeros(9)
    
    d[0] = np.multiply(9*((nM-0.5)**(4/3.0))*(variables[-3] - variables[0])*params[0],params[-1]**(1/3.0))- np.divide(variables[0],params[3]*(1+variables[0]))
    d[1] = np.multiply(9*((nM-0.5)**(4/3.0))*(variables[-2] - variables[1])*params[1],params[-1]**(1/3.0))+ np.divide(variables[0],1+variables[0]) - np.divide(params[4]*variables[1],1+variables[1])
    d[2] = np.multiply(9*((nM-0.5)**(4/3.0))*(variables[-1] - variables[2])*params[2],params[-1]**(1/3.0))+ np.divide(params[4]*variables[1],2*(1+variables[1]))

    d[-3] = np.divide(3*((nM+0.5)**(2/3.0))*(ext[0] - variables[-3])*params[5],params[-1]**(1/3.0))- np.divide(9*((nM-0.5)**(4/3.0))*(variables[-3] - variables[-6])*params[0],params[-1]**(2/3.0))
    d[-2] = np.divide(3*((nM+0.5)**(2/3.0))*(ext[1] - variables[-2])*params[6],params[-1]**(1/3.0))- np.divide(9*((nM-0.5)**(4/3.0))*(variables[-2] - variables[-5])*params[1],params[-1]**(2/3.0))
    d[-1] = np.divide(3*((nM+0.5)**(2/3.0))*(ext[2] - variables[-1])*params[7],params[-1]**(1/3.0))- np.divide(9*((nM-0.5)**(4/3.0))*(variables[-1] - variables[-4])*params[2],params[-1]**(2/3.0))


    dext = np.zeros(3)
    
    dext[0] = (-np.divide(3*((nM+0.5)**(2/3.0))*(ext[0] - variables[-3])*params[5],params[-1]**(1/3.0)))*np.divide(params[-1],params[-2])
    dext[1] = (-np.divide(3*((nM+0.5)**(2/3.0))*(ext[1] - variables[-2])*params[6],params[-1]**(1/3.0)))*np.divide(params[-1],params[-2])
    dext[2] = (-np.divide(3*((nM+0.5)**(2/3.0))*(ext[2] - variables[-1])*params[7],params[-1]**(1/3.0)))*np.divide(params[-1],params[-2])


    allder[:-3] = d
    allder[-3:] = dext
    return allder
