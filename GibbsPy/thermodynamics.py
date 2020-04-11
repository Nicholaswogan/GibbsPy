# script for gibbs energy in gas phase
import re
import numpy as np
import sys
import codecs
import os

def gibbs_helper(spec,T):
    R = 8.3145

    f = open(os.path.join(os.path.dirname(__file__), 'NEWNASA.TXT'),'r')
    lines = f.readlines()


    j = 0
    jj = 0
    for i in range(0,len(lines)):
        if lines[i][0:len(spec)]==spec:
            j+=1
            lines_spec = lines[i:i+30]
            for ii in range(0,20):
                if lines_spec[2+3*ii]=='\n':
                    break

                LB = float(lines_spec[2+3*ii].split()[0])
                UB = float(lines_spec[2+3*ii].split()[1])
                if LB<T<UB:
                    jj+=1
                    line_1 = lines_spec[2+3*ii+1].replace('D','E').replace('\n','')
                    line_2 = lines_spec[2+3*ii+2].replace('D','E').replace('\n','')
                    C = [line_1[k:k+16] for k in range(0, len(line_1), 16)]
                    C = C+[line_2[k:k+16] for k in range(0, len(line_2), 16)]
                    C = [float(k) for k in C]
                    break

    if j==0:
        print("Couldn't find it")
    if j>1:
        print("Not specific enough of an entry!")
        sys.exit()
    if jj==0:
        print('Too hot or too cold!')
        sys.exit()

    a1,a2,a3,a4,a5,a6,a7,a8,a9,a10 = C

    DH = R*T*(-a1*T**-2+a2*np.log(T)/T+a3+a4*T/2+a5*T**2/3+a6*T**3/4+a7*T**4/5+a9/T)
    DS = R*(-a1*T**-2/2-a2/T+a3*np.log(T)+a4*T+a5*T**2/2+a6*T**3/3+a7*T**4/4+a10)

    DG = DH-T*DS
    return DG

def gibbs_helper_speedy(spec,T):
    R = 8.3145

    f = open(os.path.join(os.path.dirname(__file__), 'NEWNASA_tiny.TXT'),'r')
    lines = f.readlines()


    j = 0
    jj = 0
    for i in range(0,len(lines)):
        if lines[i][0:len(spec)]==spec:
            j+=1
            lines_spec = lines[i:]
            for ii in range(0,20):
                if lines_spec[2+3*ii]=='\n':
                    break

                LB = float(lines_spec[2+3*ii].split()[0])
                UB = float(lines_spec[2+3*ii].split()[1])
                if LB<T<UB:
                    jj+=1
                    line_1 = lines_spec[2+3*ii+1].replace('D','E').replace('\n','')
                    line_2 = lines_spec[2+3*ii+2].replace('D','E').replace('\n','')
                    C = [line_1[k:k+16] for k in range(0, len(line_1), 16)]
                    C = C+[line_2[k:k+16] for k in range(0, len(line_2), 16)]
                    C = [float(k) for k in C]
                    break

    if j==0:
        print("Couldn't find it")
    if j>1:
        print("Not specific enough of an entry!")
        sys.exit()
    if jj==0:
        print('Too hot or too cold!')
        sys.exit()

    a1,a2,a3,a4,a5,a6,a7,a8,a9,a10 = C

    DH = R*T*(-a1*T**-2+a2*np.log(T)/T+a3+a4*T/2+a5*T**2/3+a6*T**3/4+a7*T**4/5+a9/T)
    DS = R*(-a1*T**-2/2-a2/T+a3*np.log(T)+a4*T+a5*T**2/2+a6*T**3/3+a7*T**4/4+a10)

    DG = DH-T*DS
    return DG


def makeAs(aa):
    items =np.array([])

    for a in aa:
        match = re.match(r"([a-z]+)([0-9]+)", a, re.I)
        if match:
            items = np.append(items,match.groups())

        else:
            items = np.append(items,[a,1])

    if len(aa)==1:
        items = np.resize(items,[2,len(aa)])
    else:
        items = np.resize(items,[2,len(aa)]).T
    return items

def convert2standard(aa,T):
    case = ['Br','I','N','Cl','H','O','F']
    for jj in range(0,len(aa)):
        for cas in case:
            if aa[jj]==cas:
                aa[jj]=aa[jj]+'2'
        if aa[jj]=='C':
            aa[jj]=aa[jj]+'(gr)'
        if aa[jj]=='S':
            if T<368.3007:
                aa[jj]=aa[jj]+'(a)'
            elif T<388.3607:
                aa[jj]=aa[jj]+'(b)'
            else:
                aa[jj]=aa[jj]+'(L)'
    return aa

def gibbsG(spec,T):

    G = gibbs_helper(spec,T)

    tmp = re.findall(r"([A-Z][a-z]?)(\d*)", spec.strip(' '))
    sp = []
    for i in range(0,len(tmp)):
        if tmp[i][1]=='':
            sp.append((tmp[i][0],1))
        else:
            sp.append((tmp[i][0],int(tmp[i][1])))

    a = []
    aa = []
    for aaa in sp:
        a.append(aaa[1])
        aa.append(aaa[0])

    aa = convert2standard(aa,T)
    speciala = makeAs(aa)

    elementg = []
    for i in range(0,len(sp)):
        spec1 = aa[i]+'  '
        ig = gibbs_helper(spec1,T)
        ig = ig/float(speciala[1,i])
        elementg.append(ig)

    G = G-np.sum(np.array(a)*np.array(elementg))
    return G


def gibbs(spec,T):
    G = gibbsG(spec,298.15)+gibbs_helper(spec,T)-gibbs_helper(spec,298.15)
    return G


def gibbsG_speedy(spec,T):

    G = gibbs_helper_speedy(spec,T)

    tmp = re.findall(r"([A-Z][a-z]?)(\d*)", spec.strip(' '))
    sp = []
    for i in range(0,len(tmp)):
        if tmp[i][1]=='':
            sp.append((tmp[i][0],1))
        else:
            sp.append((tmp[i][0],int(tmp[i][1])))

    a = []
    aa = []
    for aaa in sp:
        a.append(aaa[1])
        aa.append(aaa[0])

    aa = convert2standard(aa,T)
    speciala = makeAs(aa)

    elementg = []
    for i in range(0,len(sp)):
        spec1 = aa[i]+'  '
        ig = gibbs_helper(spec1,T)
        ig = ig/float(speciala[1,i])
        elementg.append(ig)

    G = G-np.sum(np.array(a)*np.array(elementg))
    return G


def dielectric(T,P):
    t = T-273.15
    p = P-1
    a = [-22.5713, -.032066, -.00028568, .0011832, .000027895, -.00000001476, 2300.64, -.13476]
    D0 = 4.476150
    di = np.exp((-10**6*D0+2*a[0]*p+2*a[1]*p*t+2*a[2]*p*t**2+2*a[3]*p**2+2*a[4]*p**2*t+2*a[5]*p**2*t**2+2*a[6]*t+2*a[7]*t**2)/(-(10**6)))
    return di

def gibbsAQ(spec,T,P):


    R = 8.3145

    f = open(os.path.join(os.path.dirname(__file__), 'sprons96_edited2.dat'),'r')
    lines1 = f.readlines()

    lines = []
    pp = 0
    for line in lines1:

        if line.replace(' ','').strip('\n')=='aqueousspecies':
            pp = 1
        if pp==1:
            lines.append(line)


    j = 0
    for i in range(0,len(lines)):
        if lines[i][20:].replace(' ','').strip('\n')==spec.replace(' ','') or \
        lines[i][20:].replace(' ','').strip('\n')==spec.replace(' ','')+'(0)':
            j+=1
            lines_spec = lines[i+3:i+6]
            temp = []
            for ii in range(0,len(lines_spec)):
                temp = temp+lines_spec[ii].strip('\n').split()
                temp = [float(k) for k in temp]

    if j==0:
        print("Couldn't find it")
        sys.exit()
    if j>1:
        print("Not specific enough of an entry!")
        sys.exit()


    coef = temp
    Tr = 298.15
    Pr = 1
    Psi = 2600
    Theta = 228
    Y = -5.81*10**(-5)
    Gr = coef[0]
    Hr = coef[1]
    Sr = coef[2]
    a1 = coef[3]/10
    a2 = coef[4]*100
    a3 = coef[5]
    a4 = coef[6]*10000
    c1 = coef[7]
    c2 = coef[8]*10000
    w = coef[9]*100000
    q = coef[10]
    diE = dielectric(T,P)
                        # These formulas are all parts of the formulas used in Walther's Essentials of GeoChemistry.
    G1 = Gr
    G2 = -1*Sr*(T-Tr)
    G3 = -1*c1*(T*np.log(T/Tr)-T+Tr)
    G4 = a1*(P-Pr)

    h1 = np.log( (Psi+P) / (Psi + Pr))

    G5 = a2*h1

    h2 = (1/(T-Theta))-(1/(Tr-Theta))
    h3 = (Theta-T)/Theta
    h4 = np.log(Tr*(T-Theta)/(T*(Tr-Theta)))

    G6 = -1*c2*(h2*h3-T/(Theta*Theta)*h4)
    G7 = (1/(T-Theta))*(a3*(P-Pr)+a4*h1)
    G8 = w*Y*(T-Tr)

    G = 4.184*(G1+G2+G3+G4+G5+G6+G7+G8)

    return G

def henrys_coef(key,T,P):
    R = 8.314

    H = np.exp((gibbs(key+'  ',T)-gibbsAQ(key,T,P))/(R*T))
    return H


def fugcoef(temperature,pressure,names,n):
    R = 8.314472
    if type(n) != np.ndarray:
        n=np.array(n)

    # search for species and extract data
    f = codecs.open(os.path.join(os.path.dirname(__file__),'fugacityCoefficientVariables.txt'), "r", "utf-16")
    lines=f.readlines()[1:]


    tCrit = np.array([])
    pCrit = np.array([])
    ace = np.array([])
    for name in names:
        for line in lines:
            line = line.strip('\n').split()
            if line[0]== name:
                tCrit = np.append(tCrit,float(line[1]))
                pCrit = np.append(pCrit,float(line[2])*10)
                ace = np.append(ace,float(line[3]))


    ai = .42747*(R**2)*(tCrit**2)/pCrit
    bi = .08664*R*tCrit/pCrit
    alphai = (1+(0.48508+1.55171*ace-0.15613*(ace**2))*(1-np.sqrt(temperature/tCrit)))**2

    aaTotal = 0
    aai = np.zeros([len(ai),len(ai)])
    for q in range(0,len(ai)):
        for p in range(0,len(ai)):
            aai[q,p] = (1-0)*np.sqrt(ai[q]*alphai[q]*ai[p]*alphai[p])
            aaTotal += n[q]*n[p]/np.sum(n)**2*aai[q,p]


    bTotal = (n/np.sum(n))@bi                           # Calculate the value of bTotal, the coefficient representing 'b, w/o indicies'
    A = aaTotal*pressure/((R*temperature)**2)          # Calculate the value of the coefficient representing 'A'
    BTotal = bTotal*pressure/(R*temperature)         # Calculate the value of the coefficient representing 'B, w/o indicies'
    Bi = bi*pressure/(R*temperature)                  # Calculate the value of the coefficient representing 'B, w/ indicies'
    cubic = [1,-1,(A-BTotal-BTotal**2),-A*BTotal]       # Coefficient terms to solve for the compressability factor, z.
    zi = np.roots(cubic);                                  # Use the customized method 'solveCubic' to calculate the value of z,
                                                            # the compressability factor.


    z = float('inf')

    for q in range(0,3):
        if np.isreal(zi[q]) and zi[q]>0 and abs(zi[q]-1)<abs(z-1):
            z = np.real(zi[q])

    if z==float('inf'):
        print('The value of z is either non-real or less than zero.')
        sys.exit()



    PhiA = Bi/BTotal*(z-1)                               # These are all just smaller parts of the whole thing
    PhiB = np.log(z-BTotal);
    PhiC = A/BTotal*(Bi/BTotal-2/aaTotal*(aai@(n/sum(n))))
    PhiD = np.log(1+BTotal/z)
    lnPhi = PhiA-PhiB+PhiC*PhiD

    Phi = np.exp(lnPhi);
    return(Phi)
