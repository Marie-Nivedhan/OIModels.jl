using OIFITS
using NFFT
using Unitful
using UnitfulAngles
using BeamEstimation
using Optimisers
using Plots
using Zygote
using Optim
using LinearAlgebra
plotlyjs()

## Test Interferometry_Fourier
# Accessing to the fit file in order to have u,v,x and y
filename = "/Users/marie/Downloads/R_CAR_all.fits"
oi = OIDataSet(filename);

#using BeamEstimation, building the uvplane 
uv = [BeamEstimation.uvplane(v2.ucoord .* 1u"m"  ,v2.vcoord .* 1u"m" , v2.instr.eff_wave.* 1u"m" ) for v2 in oi.vis2] 
all_uv = hcat(uv...);

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
models=[m,n,o,p]

m0=Star(0.0,0.0)
n0=Disk(0.0,0.0,0.0)
o0=Gauss(0.0,0.0,0.0)
p0=Ring(0.0,0.0,0.0,0.0)
initialmodels=[m0,n0,o0,p0]
iterator=5000000
for i in 1:4
    truemodel = models[i]
    type0 = typeof(truemodel)
    name0 = nameof(type0)
    u,de=destructure(truemodel)

    donnees=interferometry_fourier(truemodel,uv[2])
    donnees_bruitees = donnees.+ randn(size(donnees)) * 0.5

    firstmodel = initialmodels[i]
    type1 = typeof(firstmodel)
    name1 = nameof(type1)

    global found=findmodelmcmc(firstmodel,uv[2],donnees,iterator)
    v,de=destructure(found)

    global foundbruit=findmodelmcmc(firstmodel,uv[2],donnees_bruitees,iterator[i])
    w,de=destructure(foundbruit)

    if (found isa Star)
        println("l'objet est  $name1($(v[1]),$(v[2])) avec les données de $name0($(u[1]),$(u[2]))")
        println("l'objet est  $name1($(w[1]),$(w[2])) avec les données bruitées de $name0($(u[1]),$(u[2]))")
        println()
    elseif (found isa Ring)
        println("l'objet est  $name1($(v[1]),$(v[2]),$(v[3]),$(v[4])) avec les données de $name0($(u[1]),$(u[2]),$(u[3]),$(u[4]))")
        println("l'objet est  $name1($(w[1]),$(w[2]),$(w[3]),$(w[4])) avec les données bruitées de $name0($(u[1]),$(u[2]),$(u[3]),$(u[4]))")
        println()
    else
        println("l'objet est  $name1($(v[1]),$(v[2]),$(v[3])) avec les données de $name0($(u[1]),$(u[2]),$(u[3]) )")
        println("l'objet est  $name1($(w[1]),$(w[2]),$(w[3])) avec les données bruitées de $name0($(u[1]),$(u[2]),$(u[3]))")
        println()
    end 
end
