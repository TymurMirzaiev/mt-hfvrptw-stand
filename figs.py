import json
from collections import defaultdict
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mt.core import nondominated
plt.rcParams.update({"font.family":"serif","font.serif":["Liberation Serif"],"font.size":10})
S=json.load(open('results/summary.json'))
R=[json.loads(l) for l in open('results/main_own.jsonl')]
lab={'MO':'МО-ЛП','WS':'ЗС','NS':'NSGA-II'}
mk={'MO':('o','black','none'),'WS':('s','black','black'),'NS':('^','grey','none')}
# fig 1
fig,ax=plt.subplots(figsize=(5,3.4))
for m in ('MO','WS','NS'):
    pts=nondominated([tuple(p) for r in R if r['spec']==['own','RC',100,0] and r['cfg']==m for p in r['points']])
    ax.scatter([p[0] for p in pts],[p[1] for p in pts],marker=mk[m][0],edgecolors=mk[m][1],facecolors=mk[m][2],s=22,label=lab[m],linewidths=0.8)
ax.set_xlabel('Вартість, ум. од.'); ax.set_ylabel('Викиди, ум. од.'); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig('figs/fig1.png',dpi=300); fig.savefig('figs/fig1.eps')
# fig 2
fig,ax=plt.subplots(figsize=(5,3.2))
for m,ls in (('MO','-'),('WS','--'),('NS',':')):
    y=[S['main_own']['by_size'][str(n)]['hv'][m] for n in (25,50,100)]
    ax.plot([25,50,100],y,ls,marker=mk[m][0],color='black',mfc='white',label=lab[m])
ax.set_xticks([25,50,100]); ax.set_xlabel('Кількість клієнтів'); ax.set_ylabel('Нормований гіперобсяг'); ax.legend(frameon=False)
fig.tight_layout(); fig.savefig('figs/fig2.png',dpi=300); fig.savefig('figs/fig2.eps')
# fig 3
F=S['fact']['cfg']; order=['POM','-OM','P-M','--M','PO-','-O-','P--','---']
names=[o.replace('P','П').replace('O','О').replace('M','М') for o in order]
fig,ax=plt.subplots(1,2,figsize=(6.3,2.8))
ax[0].bar(range(8),[100*F[o]['feasible'] for o in order],color='white',edgecolor='black',hatch='//')
ax[0].set_ylabel('Допустимих у повній моделі, %'); ax[0].set_xticks(range(8)); ax[0].set_xticklabels(names,rotation=60,fontsize=8)
b=[100*F[o]['shares'][2] for o in order]; v=[100*F[o]['shares'][1] for o in order]; t=[100*F[o]['shares'][0] for o in order]
ax[1].bar(range(8),b,color='white',edgecolor='black',label='велосипед')
ax[1].bar(range(8),v,bottom=b,color='lightgrey',edgecolor='black',label='фургон')
ax[1].bar(range(8),t,bottom=[x+y for x,y in zip(b,v)],color='black',edgecolor='black',label='вантажівка')
ax[1].set_ylabel('Склад парку, %'); ax[1].set_xticks(range(8)); ax[1].set_xticklabels(names,rotation=60,fontsize=8); ax[1].legend(frameon=False,fontsize=7,loc='lower right')
fig.tight_layout(); fig.savefig('figs/fig3.png',dpi=300); fig.savefig('figs/fig3.eps')
print('ok')
