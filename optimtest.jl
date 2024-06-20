using OIFITS
using Unitful
using UnitfulAngles
using BeamEstimation
using Zygote
using Optim
using Optimisers

# Charger le jeu de données OIFITS
filename = "/Users/marie/Downloads/R_CAR_all.fits"
oi = OIDataSet(filename)

# Construction du uvplane à partir de BeamEstimation
uv = [BeamEstimation.uvplane(v2.ucoord .* 1u"m", v2.vcoord .* 1u"m", v2.instr.eff_wave .* 1u"m") for v2 in oi.vis2]
uv2 = stripeunits.(uv[2])  # Utilisation de uv[2] comme exemple

x=75.0*u"mas"
y=50.0*u"mas"
R=10.0 *u"mas"

x2=-75*u"mas"
y2=-75*u"mas"
R2=20.0 *u"mas"

x3=-75.0*u"mas"
y3=60.0*u"mas"
R3=50.0*u"mas"

x4=75.0*u"mas"
y4=-60.0*u"mas"
R4=30.0 *u"mas"

m=Star(x,y)
n=Disk(x2,y2,R2)
o=Gauss(x3,y3,R3)
p=Ring(x4,y4,R3,R4)

m0=Star(0.0,0.0)
n0=Disk(0.0,0.0,0.0)
o0=Gauss(0.0,0.0,0.0)
p0=Ring(0.0,0.0,0.0,0.0)
initialmodels=[m0,n0,o0,p0]

# Données d'interférométrie présumées (utilisées comme exemple)
datastar = interferometry_fourier(m,uv2)
datadisk = interferometry_fourier(n,uv2)
datagauss= interferometry_fourier(o,uv2)
dataring = interferometry_fourier(p,uv2)
datas=[datastar,datadisk,datagauss,dataring]
rng=-127.0:14.0:128.0
for i in 1:4
    global found=findmodeloptim(initialmodels[i],rng,datas[i],uv2)
    type = typeof(initialmodels[i])
    name = nameof(type)
    v,de=destructure(found)
    if (found isa Star)
        println("l'objet est  $name($(v[1]),$(v[2])) avec les données ")
        println()
    elseif (found isa Ring)
        println("l'objet est  $name($(v[1]),$(v[2]),$(v[3]),$(v[4])) avec les données ")
        println()
    else
        println("l'objet est  $name($(v[1]),$(v[2]),$(v[3])) avec les données ")
        println()
    end 
end