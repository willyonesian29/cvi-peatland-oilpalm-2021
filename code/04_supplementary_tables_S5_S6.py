"""Tables S5 (leave-one-indicator-out) and S6 (pillar decomposition by dominant plantation type).
Input: analysis_input.csv (from 00d_prepare_analysis_input.py)."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.preprocessing import minmax_scale
SENS=["air_minum","air_mck","jamban","tps","konflik","cemar","bencana","ragam_eko"]
CAP=["faskes_km","nakes","posyandu","posyandu_aktif","bpjs","sd","smp","sma","kursus","listrik_kk","sinyal","inet_desa","inet_warga","mitigasi","pusat_dagang","pasar","minimarket","logistik","bank","kredit","koperasi","angkot","jalan_mobil","jalan_mutu"]
BIN=["bumdes_bin","l_tani_bin"]
df=pd.read_csv("analysis_input.csv")
s=df[df.peat_tp_ha.notna()].drop_duplicates("id_desa").reset_index(drop=True)
print('missing:',int(s[SENS+CAP].isna().sum().sum()))
for c in SENS+CAP: s[c]=s[c].fillna(s[c].median())
s["bumdes_bin"]=(s.bumdes_tani=="Ada").astype(int); s["l_tani_bin"]=(s.l_tani=="Ada").astype(int)
s["x"]=(s.peat_tp_ha/s.VillageArea_Ha).clip(upper=1)
E0=minmax_scale(s.x); dom=np.where(s.peat_sp_ha>s.peat_ip_ha,"S","I")
def cv(sens,cap): return (E0+minmax_scale(1-s[sens].mean(axis=1))+1-minmax_scale(s[cap].mean(axis=1)))/3
base=cv(SENS,CAP+BIN)
import mapclassify
np.random.seed(42); TOPN=int((base>mapclassify.NaturalBreaks(base,k=5).bins[2]).sum())   # High + Very High
print("priority group:",TOPN)
rb=pd.Series(base).rank(ascending=False,method="first"); tb=set(rb[rb<=TOPN].index)
rows=[]
for pil,lst,full in [("Sensitivity",SENS,None),("Adaptive Capacity",CAP+BIN,None)]:
    for c in lst:
        v=cv([x for x in SENS if x!=c],CAP+BIN) if pil=="Sensitivity" else cv(SENS,[x for x in CAP+BIN if x!=c])
        r=pd.Series(v).rank(ascending=False,method="first"); tv=set(r[r<=TOPN].index)
        u,p=mannwhitneyu(v[dom=="S"],v[dom=="I"])
        rows.append([c,pil,round(spearmanr(base,v)[0],3),round(len(tb&tv)/TOPN*100,1),int(np.median(np.abs(rb-r))),
                     "Yes" if np.median(v[dom=="S"])<np.median(v[dom=="I"]) else "No",f"{p:.1e}"])
t=pd.DataFrame(rows,columns=["Removed field","Pillar","Spearman rho","Priority group retained (%)","Median |Δ rank|","Smallholder < industrial","Mann-Whitney p"]).sort_values("Spearman rho")
t.to_csv("tableS5_loo.csv",index=False); print(t.to_string(index=False))
S=minmax_scale(1-s[SENS].mean(axis=1)); A=minmax_scale(s[CAP+BIN].mean(axis=1))
out=[]
for nm,arr in [("Exposure (normalized)",E0),("Sensitivity",S),("1 − Adaptive Capacity",1-A),("CVI",base)]:
    u,p=mannwhitneyu(arr[dom=="S"],arr[dom=="I"])
    out.append([nm,f"{np.median(arr[dom=='S']):.3f}",f"{np.median(arr[dom=='I']):.3f}",f"{u:,.0f}",f"{p:.1e}"])
t6=pd.DataFrame(out,columns=["Pillar","Median, smallholder-dominated (n = %d)"%(dom=="S").sum(),"Median, industrial-dominated (n = %d)"%(dom=="I").sum(),"Mann-Whitney U","p"])
t6.to_csv("tableS6.csv",index=False); print(t6.to_string(index=False))
