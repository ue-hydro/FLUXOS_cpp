import numpy as np, glob, os, sys
d=sys.argv[1]
for f in sorted(glob.glob(d+"/*.txt"), key=lambda p:int(os.path.basename(p)[:-4])):
    a=np.loadtxt(f,delimiter=",",skiprows=1,ndmin=2); h=a[:,5]*100; r=a[:,0]
    zb=a[:,4]-a[:,5]
    print(os.path.basename(f),"S(r<=12) %7.0f  N(r>=18) %7.0f  chan %7.0f"%(h[r<=12].sum(),h[r>=18].sum(),h[(r>12)&(r<18)].sum()),
          " zb r13/r17: %.2f %.2f"%(zb[r==13].mean() if (r==13).any() else -1, zb[r==17].mean() if (r==17).any() else -1))
