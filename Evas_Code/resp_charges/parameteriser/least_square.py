import numpy as np
import math
from  get_all_dihedral import find_dihedrals
#from  energyDecomposition import singlepoint


def returnleastsquarefit(N,d,p,qme,mme,ang,type,nmax,phase,ini,w):
    '''
N int : number of data points ;
d int : number of dihedral angle to fit
p int : number of unique dihedral
qme array (N) : QM energy for each data points (from gaussian)
mme array (N): MMenergy for each data points (from Sire)
ang array (d)(N) dihedral angle value for each diheral being fit ( MDtraj)
type array (p)() : list of unique dihedral types
nmax int: maximum multiplicity being fit
ini array [2][2*p*nmax] : initial guess : parmchk gaff2 first term energy barrier
    '''
    if phase :
        npar = 2*p*nmax
    else : npar = p*nmax
    k=np.zeros(npar)
    C=np.zeros((npar,npar))
    ini_alt=np.zeros(2*p*nmax)
    for m in range(1,p*nmax):
        ini_alt[m] = ini[0][m]*math.cos(ini[0][p*nmax+m])
        ini_alt[p*nmax+m] = ini[0][m]*math.sin(ini[0][p*nmax+m])
    for tti in range (0, p):
        for n in range(1, nmax):

            k[tti*nmax+n] += ini[1][(tti)*nmax+n]*ini_alt[(tti)*nmax+n]
            if phase :
                print(tti, n , k ,p*nmax+(tti)*nmax+n )
                k[p*nmax+(tti)*nmax+n] += ini[1][(tti)*nmax+n]*ini_alt[p*nmax+(tti)*nmax+n]
            for i  in range (0, N) :
                #for ti in range ( 0 , len(type[tti])):
                for ti in range ( 0 , len(type[tti])):

                    k[tti*nmax+n] += w[i] * (qme[i] - mme[i]) * math.cos(n*ang[type[tti][ti]][i])
                    if phase :
                        k[p*nmax+(tti)*nmax+n] += w[i] * (qme[i] - mme[i]) * math.sin(n*ang[type[tti][ti]][i])
    for tti1 in range(p):
        for n1 in range(  nmax):
            #C[(tti1)*nmax+n1][(tti2-1)*nmax+n2] += ini[1][(tti1)*nmax+n1]
            C[(tti1)*nmax+n1][(tti1)*nmax+n1] += ini[1][(tti1)*nmax+n1]
            if phase :
                C[p*nmax+(tti1-1)*nmax+n1][p*nmax+(tti1-1)*nmax+n1] += ini[1][(tti1-1)*nmax+n1]
            for tti2 in range(tti1,p):
                for n2 in range(n1,nmax):
                    for i in range(1, N):
                        for ti1 in range( 1, len(type[tti1])):
                            for ti2 in range( 1,  len(type[tti2])):
                                C[(tti1-1)*nmax+n1][(tti2-1)*nmax+n2] += w[i] * math.cos(n1*ang[type[tti1][ti1]][i]) * math.cos(n2*ang[type[tti2][ti2]][i])
                            if phase :

                                C[p*nmax+(tti1-1)*nmax+n1][(tti2-1)*nmax+n2] += w[i] * math.cos(n1*ang[type[tti1][ti1]][i]) * math.sin(n2*ang[type[tti2][ti2]][i])
                                C[p*nmax+(tti1-1)*nmax+n1][p*nmax+(tti2-1)*nmax+n2] += w[i] * math.sin(n1*ang[type[tti1][ti1]][i]) * math.sin(n2*ang[type[tti2][ti2]][i])
                                C[(tti2-1)*nmax+n2][(tti1-1)*nmax+n1] = C[(tti1-1)*nmax+n1][(tti2-1)*nmax+n2]
            if phase :
                C[(tti2-1)*nmax+n2][p*nmax+(tti1-1)*nmax+n1] =   C[p*nmax+(tti1-1)*nmax+n1][(tti2-1)*nmax+n2]
                C[p*nmax+(tti2-1)*nmax+n2][p*nmax+(tti1-1)*nmax+n1] = C[p*nmax+(tti1-1)*nmax+n1][p*nmax+(tti2-1)*nmax+n2]
    C_inv = C.T
    delta=np.zeros(2*p*nmax)
    for m in range (1 , npar) :
        for l in range(1, npar):
            delta[m] += C_inv[m][l] * k[l]
    a_alt=np.zeros(2*p*nmax)
    for m in range(1, 2*p*nmax):
        a_alt[m] = delta[m]

        if ini[m] == 0 :
            a_alt[m] += ini_alt[m]
    a=np.zeros(2*p*nmax)
    for m in range(1, p*nmax):
        print(a_alt[m]**2)
        a[m] = math.sqrt(a_alt[m]**2 + a_alt[p*nmax+m]**2)
        a[p*nmax+m] = math.atan(a_alt[p*nmax+m]/ a_alt[m])
    return a
all_dihedrals,   all_dihedrals_type, dihedrals_heavy, dihedrals_heavy_name, torsion_names, dihedrals_heavy_index=find_dihedrals('newcharges.mol2' )

print(all_dihedrals)

def make_traj(prefix, dihedrals_heavy_index):
    list_pdb=open(prefix +'-listfiles')
    lines=list_pdb.readlines()
    if len(lines) > 0 : first=lines.pop(0).replace(' \n','')
    else :
        sys.exit('BOUH!' + prefix +   str(lines) )
        return False
    t = md.load(first, top='input.prmtop')

    for i in range(len(lines)):
        filename=lines.pop(0).replace(' \n','')
        print(filename)
        tnext= md.load(filename.replace(' \n',''), top='input.prmtop')
        t= md.join([t , tnext], check_topology=True, discard_overlapping_frames=False)
    t.save_mdcrd(prefix+ '-traj.mdcrd' )
    t.save_netcdf(prefix+ '-traj.nc' )
    return t.n_frames


N=5
d=3
p=2
qme=[1.01,1.02,1.03,1.00,2]
mme=[1.11,1.22,1.53,1.00,2]
ang = [[5,5,5,6,1],[10,20,30,45,12], [10,50,60,12,27]]  # size d*N
type=[[1],[1,2]]

nmax=2
phase=True
ini=[]
ini.append([2,2,2,2,2,2,2,0,2,0])
ini.append([2,2,2,2,2,2,2,0,2,0])
for h in range(2*p*nmax) : ini.append(0)
w=np.ones(N)
returnleastsquarefit(N,d,p,qme,mme,ang,type,nmax,phase,ini,w)
