using Unitful, UnitfulAngles
using Plots

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
    global t=imagingxy(models[i],pixsize,N)
    global valpiximg=sum(t)
    heatmap(t)
    type = typeof(models[i])
    name = nameof(type)
    println("Valeur de l'image réelle de $name est $valpiximg")
    println()
    savefig("L'image réelle de $name $i .png")
end











###FAUT MOFIFIER LES VERIFICATION
#Verification for Star

#= for i in 1:N
    for j in 1:N
        if s[i,j]!=0
            println((j,i))
        end
    end
end  =#

#Verification for Disk
#= 
diamx=Vector{Int64}()
diamy=Vector{Int64}()
t=imagingxy(Star(n.position),pixsize,N)

for i in 1:N
    for j in 1:N
        if t[i,j]!=0
            global x0=j
            global y0=i
        end
    end
end

for k in 1:N
    if s[y0,k]!=0
        push!(diamx,k)
    end
    if s[k,x0]!=0
        push!(diamy,k)
    end
end
if (x0-diamx[1])*pixsize*u"mas^-1"==R*u"mas^-1"
    println("diamètre selon x est bon")
end
if (y0-diamy[1])*pixsize*u"mas^-1"==R*u"mas^-1"
    println("diamètre selon y est bon")
end 
 =#
#Verification for Gauss
#= diamx=Vector{Int64}()
diamy=Vector{Int64}()
t=imagingxy(Star(o.position),pixsize,N)

for i in 1:N
    for j in 1:N
        if t[i,j]!=0
            global x0=j
            global y0=i
        end
    end
end

for k in 1:N
    if !isapprox(s[y0,k], 0, atol=0.0005)
        push!(diamx,k)
    end
    if !isapprox(s[y0,k], 0, atol=0.0005)
        push!(diamy,k)
    end
end
if round(Int64,length(diamx)*pixsize*u"mas^-1")==R*2*u"mas^-1"
    println("diamètre selon x est bon")
end
if round(Int64,length(diamy)*pixsize*u"mas^-1")==R*2*u"mas^-1"
    println("diamètre selon y est bon")
end
 =#

#Verification for Ring

#= diamx=Vector{Int64}()
diamy=Vector{Int64}()
t=imagingxy(Star(p.position),pixsize,N)

for i in 1:N
    for j in 1:N
        if t[i,j]!=0
            global x0=j
            global y0=i
        end
    end
end

for k in 1:N
    if s[y0,k]!=0
        push!(diamx,k)
    end
    if s[k,x0]!=0
        push!(diamy,k)
    end
end

if (x0-diamx[1])*pixsize*u"mas^-1"==R*u"mas^-1"
    println("diamètre selon x est bon")
end
if (y0-diamy[1])*pixsize*u"mas^-1"==R*u"mas^-1"
    println("diamètre selon y est bon")
end 




println(valpiximg)
heatmap(t) =#