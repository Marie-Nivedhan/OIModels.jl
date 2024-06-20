module OIModels

using Unitful,
    UnitfulAngles,
    NFFT,
    SpecialFunctions,
    Functors,
    Optimization,
    Optim,
    ArrayTools,
    Zygote,
    Turing,
    StatsPlots,
    Distributions,
    ForwardDiff,
    LinearAlgebra

export Model,
           Star,
           Disk,
           Gauss,
           Ring,
           Twomodels,
           interferometry_fourier,
           interferometry_image,
           imagingdata,
           imagingxy,
           imaginguv,
           findimgmcmc,
           findimg,
           findobj,
           findmodelmcmc,
           findmodeloptim,
           stripeunits
           

abstract type Model
end

interferometry_image(M::Model,x,y; atol=0.5) = error("interferometry_image not implemented")
interferometry_fourier(M::Model,x,y) = error("interferometry_fourier not implemented")
interferometry_fourier(M::Model,uv) = interferometry_fourier(M,uv[1,..],uv[2,..])


###########################################

###Definition of subtype Star, defined by his position coordinates ('x','y')
struct Star{T} <: Model
	position::T
end
@functor Star
Star(x,y) = Star([x,y])

##Functions for model Star

#Definition of interferometry_fourier(M::Star,u,v)
"""
    This fuction calculates the fourier transform of a star 'M' for a spectrum data coordinates ('u','v'). 
    It will be used to compare the experiments and the model.

    Parameters
    ----------
    M : Star
        object that represents the star
    u : Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}}
        horizontal fourier space coordinates
    v : Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}}
        vertical fourier space coordinates
    
    Returns
    -------
    result : Vector{ComplexF64}
    vector that contains the results of the normalized fourier transform
"""
function interferometry_fourier(M::Star,u,v)
    uu=stripeunits(u)
    vv=stripeunits(v)
    x,y = stripeunits(M.position)
    F=exp.(-im*2*π*(uu*x.+vv*y))
    return F

end

#Definition of interferometry_image(M::Star,x,y; atol=0.5)
"""
    This fuction compares the position of a star 'M' with the a spatial data coordinates of his image ('x','y').
    This comparison has a default tolerance atol which equals to 0.5 that represents a half pixel.  
    It will be used to compare the experiments and the model.

    Parameters
    ----------
    M : Star
        object that represents the star
    x : Float64
        horizontal space coordinates of the image
    y : Float64
        vertical space coordinates of the image
    atol : Float64
        tolerance defaultly set as 0.5
    
    Returns
    -------
    result : Float64
    This float is either 1 or 0 is the conversion of the booleen that compares the experiment to the model
"""

function interferometry_image(M::Star,x,y; atol=0.5)

    Mx,My = M.position
    f=Float64.(isapprox.(x, Mx; atol=atol) .&& isapprox.(y, My; atol=atol))
    return f

end


############################################

###Definition of subtype Disk, defined by the space coordinates of his center('x','y') and his radius 'R'
struct Disk{T,S} <: Model
	position::T 
    Radius::S
end
@functor Disk
Disk(x,y,R) = Disk([x,y],[R])
##Functions for model Disk

#Definition interferometry_fourier(M::Disk,u,v)
"""
    Interferometry_Fourier(M::Disk,u,v) 
    This fuction calculates the fourier transform of a uniform Disk 'M' for a spectrum data coordinates ('u','v'). 
    It will be used to compare the experiments and the model.

    Parameters
    ----------
    M : Disk
        object that represents the uniform disk
    u : Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}}
        horizontal fourier space coordinates
    v : Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}}
        vertical fourier space coordinates
    
    Returns
    -------
    result : Vector{ComplexF64}
    vector that contains the results of the normalized fourier transform
"""
function interferometry_fourier(M::Disk,u,v)

    x,y = stripeunits.(M.position)
    R=stripeunits.(M.Radius[1])
    uu=stripeunits.(u)
    vv=stripeunits.(v)
    ρ=sqrt.(uu.^2 .+ vv.^2)
    F= zeros(Complex{Float64}, size(uu)...)
    for idx in eachindex(uu)
        if ρ[idx]==0
            F[idx]=1 #bessel cardinal =1/2 when x->0
        else
            F[idx]=2 * besselj1(2*π*R*ρ[idx]) / 
                    (2*π*R*ρ[idx]) *
                    exp(-im*2*π*(uu[idx]*x+vv[idx]*y))# la phase 
        end
    end
    F0=1 
    return F./F0

end

#Definition of a interferometry_image(M::Disk,x,y; atol=0.5)
"""
    This fuction compares the position of a Disk 'M' with the a spatial data coordinates of his image ('x','y').
    This comparison has a default tolerance atol which equals to 0.5 that represents a half pixel.  
    It will be used to compare the experiments and the model.

    Parameters
    ----------
    M : Disk
        object that represents a uniform disk
    x : Float64
        horizontal space coordinates of the image
    y : Float64
        vertical space coordinates of the image
    atol : Float64
        tolerance defaultly set as 0.5
    
    Returns
    -------
    result : Float64
    This float is either 1 or 0 is the conversion of the booleen that compares the experiment to the model
"""
function interferometry_image(M::Disk,x,y;atol=0.5) 
    
    Mx,My=M.position
    MR=M.Radius[1]
    f=Float64.(sqrt.((y .- My).^2+(x .- Mx).^2).<MR)
    f0=π*(MR*u"mas^-1"|>NoUnits)^2
    return f./f0

end


##############################################

###Definition of subtype Gauss, defined by his space coordinates of his center ('x','y') and his full width at half maximum 'R'
struct Gauss{T,S} <: Model
	position::T 
    Radius::S 
end
@functor Gauss
Gauss(x, y, R) = Gauss([x,y],[R])
##Functions for model Gauss

#Definition interferometry_fourier(M::Gauss,u,v)
""" 
    This fuction calculates the fourier transform of a gaussian form 'M' for a spectrum data coordinates ('u','v'). 
    It will be used to compare the experiments and the model.

    Parameters
    ----------
    M : Gauss
        object that represents the gaussian form
    u : Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}}
        horizontal fourier space coordinates
    v : Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}}
        vertical fourier space coordinates
    
    Returns
    -------
    result : Vector{ComplexF64}
    vector that contains the results of the normalized fourier transform
"""
function interferometry_fourier(M::Gauss,u,v)

    x,y = stripeunits.(M.position)
    R=stripeunits.(M.Radius[1])
    uu=stripeunits.(u)
    vv=stripeunits.(v)
    ρ=(uu.^2 .+ vv.^2) 
    F=exp.(-((π*R)^2 .*ρ)./ 4*log(2)).* exp.(-im*2*π*(uu.*x.+vv.*y))
    F0=1
    return F./F0 

end

#Definition of a interferometry_image(M::Gauss,x,y; atol=0.5)
"""
    This fuction compares the position of a Gauss 'M' with the a spatial data coordinates of his image ('x','y').
    This comparison has a default tolerance atol which equals to 0.5 that represents a half pixel.  
    It will be used to compare the experiments and the model.

    Parameters
    ----------
    M : Disk
        object that represents a gaussian form object
    x : Float64
        horizontal space coordinates of the image
    y : Float64
        vertical space coordinates of the image
    atol : Float64
        tolerance defaultly set as 0.5
    
    Returns
    -------
    result : Float64
    This float is either 1 or 0 is the conversion of the booleen that compares the experiment to the model
"""
function interferometry_image(M::Gauss,x,y;atol=0.5)
    
    Mx,My = M.position
    R=M.Radius[1]
    f=exp.(-4*log(2).*((x.-Mx).^2 .+(y.-My).^2) ./R^2)
    f0=(π*(R*u"mas^-1"|>NoUnits)^2)/(4*log(2))
    return f./f0

end


###############################################

###Definition of subtype Ring, defined by his space coordinates of his center ('x','y') and his internal and external radiuses 'R1' 'R2'
struct Ring{T} <: Model
	position::T 
    Radius::T 
end
@functor Ring
Ring(x,y,R1,R2) = Ring([x,y],[R1,R2])
##Functions for model Ring

#Definition interferometry_fourier(M::Ring,u,v)
""" 
    This fuction calculates the fourier transform of a ring form 'M' for a spectrum data coordinates ('u','v'). 
    It will be used to compare the experiments and the model.

    Parameters
    ----------
    M : Ring
        object that represents the ring form
    u : Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}}
        horizontal fourier space coordinates
    v : Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}}
        vertical fourier space coordinates
    
    Returns
    -------
    result : Vector{ComplexF64}
    vector that contains the results of the normalized fourier transform
"""
function interferometry_fourier(M::Ring,u,v)

    x,y = stripeunits.(M.position)
    if M.Radius[1]<M.Radius[2]
        R1 = M.Radius[1]
        R2 = M.Radius[2]
    else
        R1 = M.Radius[2]
        R2 = M.Radius[1]
    end      
    R1=stripeunits.(R1)
    R2=stripeunits.(R2)
    uu=stripeunits.(u)
    vv=stripeunits.(v)
    ρ=stripeunits.(sqrt.(uu.^2 .+ vv.^2))
    F= zeros(Complex{Float64}, size(uu)...)
    for idx in eachindex(uu)
        if ρ[idx]==0
            F[idx]=π * (R2^2 - R1^2)
        else
            F[idx]=(R2*besselj1(ρ[idx]*2*π*R2)-R1*besselj1(ρ[idx]*2*π*R1))*exp(-im*2*π*(uu[idx]*x+vv[idx]*y))/ρ[idx]
        end
    end
    F0=π * (R2^2 - R1^2)
    return F./F0

end

#Definition of a interferometry_image(M::Ring,x,y; atol=0.5)
"""
    This fuction compares the position of a Ring 'M' with the a spatial data coordinates of his image ('x','y').
    This comparison has a default tolerance atol which equals to 0.5 that represents a half pixel.  
    It will be used to compare the experiments and the model.

    Parameters
    ----------
    M : Ring
        object that represents a gaussian form object
    x : Float64
        horizontal space coordinates of the image
    y : Float64
        vertical space coordinates of the image
    atol : Float64
        tolerance defaultly set as 0.5
    
    Returns
    -------
    result : Float64
    This float is either 1 or 0 is the conversion of the booleen that compares the experiment to the model
"""
function interferometry_image(M::Ring,x,y;atol=0.5)
    
    Mx,My = M.position
    if M.Radius[1]<M.Radius[2]
        R1 = M.Radius[1]
        R2 = M.Radius[2]
    else
        R1 = M.Radius[2]
        R2 = M.Radius[1]
    end      
    f = Float64.(sqrt.((y .- My).^2+(x .- Mx).^2).<R2 .&& sqrt.((y .- My).^2+(x .- Mx).^2).>R1 )
    f0 = π * (R2^2 - R1^2)*u"mas^-2"|>NoUnits
    return f./f0

end


################################################

###Definition of subtype Twomodels, defined by two models 'M' and their flux ratios 'fluxm' and 'fluxn'
struct Twomodels{} <: Model
	models::Vector{Model}
    flux::Vector{Float64}
end
@functor Twomodels
Twomodels(m,n,fluxm,fluxn) = Twomodels([m,n],[fluxm,fluxn])
##Functions for model Twomodels

#Definition of interferometry_fourier(M::Twomodels,u,v)
""" 
    This fuction calculates the fourier transform of two objects at a time represented by 'M' for a spectrum data coordinates ('u','v'). 
    It will be used to compare the experiments and the model.

    Parameters
    ----------
    M : Twomodels
        object that represents two objects
    u : Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}}
        horizontal fourier space coordinates
    v : Vector{Quantity{Float64, NoDims, Unitful.FreeUnits{(rad^-1,), NoDims, nothing}}}
        vertical fourier space coordinates
    
    Returns
    -------
    result : Vector{ComplexF64}
    vector that contains the results of the normalized fourier transform
"""
function interferometry_fourier(M::Twomodels,u,v)

    m,n = M.models
    fluxm,fluxn =M.flux
    Fm=interferometry_fourier(m,u,v)
    Fn=interferometry_fourier(n,u,v)
    F=fluxm.*Fm+Fn.*fluxn
    return F 

end

#Definition of a interferometry_image(M::Twomodels,x,y; atol=0.5)
"""
    This fuction compares the position of a Twomodels 'M' with the a spatial data coordinates of his image ('x','y').
    This comparison has a default tolerance atol which equals to 0.5 that represents a half pixel.  
    It will be used to compare the experiments and the model.

    Parameters
    ----------
    M : Twomodels
        object that represents two objects
    x : Float64
        horizontal space coordinates of the image
    y : Float64
        vertical space coordinates of the image
    atol : Float64
        tolerance defaultly set as 0.5
    
    Returns
    -------
    result : Float64
    This float is either 1 or 0 is the conversion of the booleen that compares the experiment to the model
"""
function interferometry_image(M::Twomodels,x,y;atol=0.5)
    
    m,n = M.models
    fluxm,fluxn =M.flux
    fm=interferometry_image(m,x,y;atol)
    fn=interferometry_image(n,x,y;atol)
    f=fluxm.*fm+fn.*fluxn
    f0=sum(f)
    return f./f0

end


#################################################

include("OIModelsfunction.jl")

end # module OIModels


