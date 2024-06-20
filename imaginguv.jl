using Unitful, UnitfulAngles, Plots
x=75.0*u"mas"
y=50.0*u"mas"
R=10.0 *u"mas"

x2=-75*u"mas"
y2=-75*u"mas"
R2=20.0 *u"mas"

x3=-75.0*u"mas"
y3=60.0*u"mas"
R3=50.0 *u"mas"

x4=75.0*u"mas"
y4=-60.0*u"mas"
R4=30.0 *u"mas"

m=Star(x,y)
n=Disk(x2,y2,R2)
o=Gauss(x3,y3,R3)
p=Ring(x4,y4,R3,R4)
q=Twomodels(m,n,0.02,0.98)
r=Twomodels(o,q,0.8,0.2)
s=Twomodels(p,r,0.2,0.8)
models=[m,n,o,p,q,r,s]
N=256
pixsize=1 *u"mas"
for i in 1:length(models)
    global t=imaginguv(models[i],pixsize,N)
    global valpiximg=sum(t)
    heatmap(t)
    type = typeof(models[i])
    name = nameof(type)
    println("Valeur de l'image de fourier de $name est $valpiximg")
    println()
    savefig("L'image de fourier de $name $i .png")
end 