import csv, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
Q={'SER-01':(149.16,204.13,470.79),'SER-02':(67.0,90.99,207.81),'SER-03':(19.96,26.26,55.95),'SER-04':(428.99,579.24,1306.45)}
NM={'SER-01':'SER-01 Agioi Anargyroi','SER-02':'SER-02 Agia Varvara','SER-03':'SER-03 Agios Ioannis','SER-04':'SER-04 Krousovitis'}
col={'C0':'#e76f51','C1':'#264653','Cp':'#2a9d8f'}; lab={'C0':'C0 no channel','C1':'C1 assumed channel','Cp':'C+ larger channel'}
def qN(q,b,N):
    q=np.array(q,float); b=np.array(b,float)
    if b[0]>=N: return f"<{q[0]:.0f}"
    i=np.argmax(b>=N)
    if b[i]<N: return f">{q[-1]:.0f}"
    return f"{q[i-1]+(N-b[i-1])/(b[i]-b[i-1])*(q[i]-q[i-1]):.0f}"
rows=[]; fig,ax=plt.subplots(1,4,figsize=(16,4.2))
for k,c in enumerate(Q):
    S=list(csv.DictReader(open(f'{c}/D3/{c}_D3_steps.csv')))
    for ch in ['C0','C1','Cp']:
        s=[r for r in S if r['channel']==ch]; q=[float(r['Q_outlet']) for r in s]; b=[int(r['bld010']) for r in s]; b5=[int(r['bld050']) for r in s]
        ax[k].plot(q,b,'-o',ms=3,color=col[ch],label=lab[ch]); ax[k].plot(q,b5,'--',color=col[ch],lw=0.9)
        rows.append(dict(case=c,channel=ch,Q_T50=Q[c][0],Q_T100=Q[c][1],Q_T1000=Q[c][2],
            Q_10bld=qN(q,b,10),Q_50bld=qN(q,b,50),Q_100bld=qN(q,b,100),Q_10bld_05m=qN(q,b5,10),
            bld_at_T50=int(round(np.interp(Q[c][0],q,b))),bld_at_T100=int(round(np.interp(Q[c][1],q,b)))))
    for t,v in zip(['T50','T100','T1000'],Q[c]): ax[k].axvline(v,color='0.6',lw=0.7,ls=':'); ax[k].text(v,ax[k].get_ylim()[1]*0.97,t,fontsize=7,ha='right',rotation=90,va='top',color='0.4')
    ax[k].set_title(NM[c],fontsize=9); ax[k].set_xlabel('outlet discharge (m³/s), steps of the T1000 hydrograph'); ax[k].grid(alpha=0.3)
ax[0].set_ylabel('buildings reached (solid ≥0.1 m, dashed ≥0.5 m)'); ax[0].legend(fontsize=7)
plt.tight_layout(); plt.savefig('fig_D3_buildings_vs_Q.png',dpi=120)
with open('D3_thresholds_summary.csv','w',newline='') as fo: w=csv.DictWriter(fo,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
for r in rows: print(r)
