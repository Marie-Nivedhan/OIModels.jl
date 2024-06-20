using OIFITS
using NFFT
using Unitful
using UnitfulAngles
using BeamEstimation
using Plots
plotlyjs()

## Test Interferometry_Fourier
# Accessing to the fit file in order to have u,v,x and y
filename = "/Users/marie/Downloads/R_CAR_all.fits"
oi = OIDataSet(filename);

#using BeamEstimation, building the uvplane 
uv = [BeamEstimation.uvplane(v2.ucoord .* 1u"m"  ,v2.vcoord .* 1u"m" , v2.instr.eff_wave.* 1u"m" ) for v2 in oi.vis2] 
all_uv = hcat(uv...);


diameter = max(1.8,maximum(oi.array[1].diameter))*1u"m"

λ =  [v2.instr.eff_wave  .* 1u"m" .|> u"μm" for v2 in oi.vis2] 
λmax = maximum(maximum.(λ))
λmin = minimum(minimum.(λ))
#Estimating the beam size in mas
rx,ry,θ = beamellipse(all_uv)

#Fixing the center of the beam as the position in rad of the Star 
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
    global t=imagingdata(models[i],uv[2],pixsize,N)
    global valpiximg=sum(t)
    heatmap(t)
    type = typeof(models[i])
    name = nameof(type)
    println("Valeur du dirty de $name est $valpiximg")
    println()
    savefig("Dirty image de $name $i .png")
end 