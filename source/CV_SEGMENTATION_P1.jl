### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 437aabde-e7f5-4ed4-8f07-b91ecf5221ba
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ 10028c35-cf65-43d5-b8f8-4a160416fa5b
begin
	using Images, Plots
	using TestImages, DSP
	using StatsBase
	using PlutoUI
end

# ╔═╡ 59ef7fc4-4f90-4154-8dd8-0260f25ace12
begin
	img = testimage("livingroom") 
end;	

# ╔═╡ 90343e63-6739-421d-99e0-bf3bc0289fa5
md"# Segmentation des images"

# ╔═╡ 38bf4f37-d915-4375-90eb-3de73fc9257a
md"""
__1. Definition :__

La segmetation est l'operation qui consiste a partitionner les pixels de l'image en plusieurs ensembles en se basant sur un critere de similarite ou de disparite. Dans le premier cas, il s'agit de la segmentation oriente region et dans le deuxieme la segmentation oriente contour  
"""

# ╔═╡ 3d15739d-f368-4717-8b60-3b2353bf1bfe
md"""
## Segmentation oriente contour

Elle est realise en procedant de la maniere ci-apres:

- Calcul de l'image gradient en utilisant des noyaux derivatifs appropries;
- Estimation des maximas d'intensite dans la direction normal du gardient.
- Seuillage de l'image des maximas.
- post-traitement de l'image seuillee pour generer les frontieres des regions de l'image.
"""

# ╔═╡ 93f35d91-6dbc-4a12-a176-d397f335fba6
md" ### Image du gradient" 

# ╔═╡ 57bd4f49-44fa-48a5-b0c1-44b42d42af58
md"""
#### Methodes spatiales

##### Differences finies et noyaux de convolution

__1. Schema de differences finis en avant :__
```math
\begin{align*}
\dfrac{\partial I}{\partial x} [i,j] &\approx I[i,j+1] - I[i,j]\\
\dfrac{\partial I}{\partial y} [i,j] &\approx I[i+1,j] - I[i,j]
\end{align*}
```
$([0 0 0;-1 1 0; 0 0 0] |> display)
"""

# ╔═╡ e9c02c13-630b-4991-89cc-1ef91b1ef82b
let
	Dfx = [0 0 0;-1.0 1.0 0; 0 0 0] |> centered
	img_dfx = imfilter(view(img,:,:),reflect(Dfx), "replicate")
	img_dfy = imfilter(view(img,:,:),reflect(Dfx'), "replicate")
	[10img_dfx img 10img_dfy]
end

# ╔═╡ 5473a579-e4b0-4dbc-b277-9096bb992d4c
md"""
__2. Schema de differences centrees :__
```math
\begin{align*}
\dfrac{\partial I}{\partial x} [i,j] &\approx I[i,j+1] - I[i,j-1]\\
\dfrac{\partial I}{\partial y} [i,j] &\approx I[i+1,j] - I[i-1,j]
\end{align*}
```
$([0 0 0;-1 0 1; 0 0 0] |> display)
"""

# ╔═╡ 966a5176-c131-4dc0-85e8-faa9bc79e327
let
	kern = [0 0 0; -1 0 1; 0 0 0] |> centered
	img_dfx = imfilter(view(img,:,:),kern,"replicate")
	img_dfy = imfilter(view(img,:,:), kern',"replicate")
	[10img_dfx img 10img_dfy]
end

# ╔═╡ 7bc3bd10-06df-49df-bcf0-9c72d4da6ef2
md"""
__3. Effet du bruit sur l'image du gradient__

_3.1 Bruit impulsionnel :_

density = $(@bind den Slider(0.025:0.025:0.2,default=0.1,show_value=true))
"""

# ╔═╡ 4ff2504e-07b4-4cfc-a4ad-8454989a78d1
let
	nsamples = ceil(Int,den * *(size(img)...))
	sp_indices = sample(CartesianIndices(img)[:],nsamples)
	salt_idx = sample(sp_indices, nsamples÷2)
	pepper_idx = Vector{eltype(salt_idx)}()
	for idx in sp_indices
		idx ∉ salt_idx && push!(pepper_idx, idx)
	end
	img_sp = copy(img)
	img_sp[salt_idx] .= 1N0f8
	img_sp[pepper_idx] .= 0N0f8
	kern = [0 0 0;-1.0 0 1.0;0 0 0] |> centered
	img_dfx = imfilter(view(img_sp,:,:),kern,"replicate")
	img_dfy = imfilter(view(img_sp,:,:), kern',"replicate")
	[10img_dfx img_sp 10img_dfy]
end

# ╔═╡ 6b6e5883-0af1-4c69-95a5-0d521ee5bf05
md"""
_3.2 Bruit gaussien :_

σ = $(@bind σ Slider([0.01, 0.025, 0.05, 0.075, 0.1],default=0.05,show_value=true))
"""

# ╔═╡ 4a8226e4-1237-434e-be0f-a51ce55b5687
let
	using Random
	img_gn = view(img,:,:) .+ σ * randn(MersenneTwister(1234), size(img))
	kern = [0 0 0;-1.0 0.0 1.0; 0 0 0] |> centered
	img_dfx = imfilter(view(img_gn,:,:),kern,"replicate")
	img_dfy = imfilter(view(img_gn,:,:), kern',"replicate")
	[10img_dfx img_gn 10img_dfy]
end

# ╔═╡ b550f65e-f97f-4aba-a48f-02694628c304
md"""
__4. Schema de differences ameliorees:__

+ Calcul de differences dans une direction
+ Lissage dans la direction orthogonale
"""

# ╔═╡ 3ece7e97-28e7-460d-84f4-9d1905585271
md"""
_4.1 Lissage par filtre a boite 1D (Filtre de Prewitt):_
$(Kernel.prewitt()[1] |> display)
"""

# ╔═╡ 43577282-bd3b-4d1b-92e9-dea1dd5a853c
let
	img_gn = view(img,:,:) .+ σ * randn(MersenneTwister(1234), size(img))
	img_dx = imfilter(view(img_gn,:,:),Kernel.prewitt()[1],"replicate")
	img_dy = imfilter(view(img_gn,:,:),Kernel.prewitt()[2],"replicate")
	[6img_dx img_gn 6img_dy]
end

# ╔═╡ 1a4a993c-fa45-47de-a400-973f53f2ca5c
md"""
_4.2 Lissage par gaussienne (Filtre de Sobel) :_
$(Kernel.sobel()[1] |> display)
"""

# ╔═╡ d82645ff-31e7-43a0-84ee-4dd612e40972
let 
	img_gn = view(img,:,:) .+ σ * randn(MersenneTwister(1234), size(img))
	img_dx = imfilter(view(img_gn,:,:),Kernel.sobel()[1],"replicate")
	img_dy = imfilter(view(img_gn,:,:),Kernel.sobel()[2],"replicate")
	[8img_dx img_gn 8img_dy]
end

# ╔═╡ ab64a060-8f14-44bd-8604-580575e1f31b
md"""
__5. Convolution et derivation__

Considerons un filtre de lisaage gaussien ``h_g(x,y)`` de variance ``\sigma^2``. On l'applique a l'image ``I(x,y)``  
```math
I_{Blur} = h_g(x,y) * I(x,y) 
```
Calculons la derivee de l'image lissee
```math
\mathbf{\nabla}I_{Blur}(x.y) = \mathbf{\nabla}h_g(x.y) * I(x,y)
```
ce qui mene a:
```math
\begin{align*}
\mathbf{\nabla}I_{Blur}(x.y) &= -\frac{1}{2\pi\sigma^4}\exp(-\frac{x^2+y^2}{2\sigma^2})\begin{bmatrix}x\\y\\\end{bmatrix} * I(x,y)
\end{align*}
```

"""

# ╔═╡ 3220f95e-74a9-4578-97ba-0e1123731ef3
md"""
#### Methodes spectrales

##### Filtre passe-haut

- exemple : Butterworth

```math
h(u,v) = \frac{1}{\sqrt{1-(D_0/D)^{2n}}}\;\text{avec } D=\sqrt{u^2+v^2} 
```
FilterType = $(@bind ftype TextField(default="disque"))
D₀ = $(@bind D₀ Slider(0.05:0.05:0.45, default=0.25, show_value=true))
n = $(@bind n Slider(1:5, default=2, show_value=true))

"""

# ╔═╡ 5d735d83-0556-470b-920d-1a98c7bf5c7f
let
	fft = DSP.fft
	fftshift = DSP.fftshift
	img_gn = (reinterpret(N0f8,img) .|> float) .+ σ * randn(MersenneTwister(1234),Float32,size(img))
	IMG = fftshift(fft(view(img_gn,:,:)))
	U = repeat(LinRange(-0.5,0.5,size(img,1)),1,size(img,1))
	V = U'
	D = ones(size(img))
	D[U.^2+V.^2 .< D₀^2] .= 0.0
	BW = 1 ./ sqrt.(1 .+ (D₀^2 ./(U.^2 + V.^2)).^n)
	f = ftype=="butterworth" ? view(BW,:,:) : view(D,:,:)
	img_gradient = DSP.ifft(DSP.ifftshift(IMG .* f)) .|> real .|> abs .|> Gray
	
	[Gray.(img_gn) 5img_gradient]
	
	#[img_gn img_gradient]
	#[img 10(img_dx .* cos.(α) + img_dy .* sin.(α))]
	#[10abs.(img_dx) img_gn 10abs.(img_dy)]
end

# ╔═╡ 55bf3ae5-d514-4af3-a269-0262c2d74f0a
md"### Cartographie des maxima"

# ╔═╡ bebb8072-33c7-468a-a9e9-828c028d1294
md"""
#### Maxima dans la direction du gradient

```math
\mathbf{\nabla}_n I(x,y) = \mathbf{\nabla} I(x,y) \cdot \mathbf{n}
```
"""

# ╔═╡ 877552e7-6ac5-4bd2-9b8d-bd2a4b6a87b6
md"""
#### Simplification du calcul de la cartographie des maxima : 

__1. norme 1 :__

```math
\left|\nabla I(x,y)\right| = \left|\dfrac{\partial}{\partial x}I(x,y)\right| + \left|\dfrac{\partial}{\partial y}I(x,y)\right|
``` 
"""

# ╔═╡ 8f4b9bcf-09c4-4750-aec5-337d6c8a065a
md"""
__2. Norme 2 :__

```math
|\mathbf{\nabla}I(x,y)| = \sqrt{\left(\frac{\partial I}{\partial x} (x,y)\right)^2 + \left(\frac{\partial I}{\partial y} (x,y)\right)^2 }
```
"""

# ╔═╡ 709d40c3-d185-4e06-9fc8-e61dc6a6a680
md"""
__3. Calcul du gradient de la cartographie des maxima :__

```math
\mathbf{\nabla}^2 I(x,y) = \mathbf{\nabla} \cdot \mathbf\nabla I(x,y)
```
"""

# ╔═╡ a4dde3e2-b13f-4a9c-8fbd-a83565cff09d
let
	x = repeat(-2:2,1,length(-2:2))
	y = copy(x')
	myLoG = centered((x.^2+y.^2)/0.75^4 .- 2/0.75^2).*Kernel.gaussian(0.75)
	img = testimage("livingroom")
	img += σ * randn(MersenneTwister(1234), size(img)) 
	[img 10imfilter(img,myLoG,"replicate")]
end

# ╔═╡ d31b3819-e2a3-4fec-bd73-1f6334990228
md"""
__4. Supression des non-maximas :__

_4.1 Procedure :_
- Calculez de la direction θ(x,y) pointee par le gradient pour chaque pixel de l'image :
- Dans cette direction _estimez_ la norme du gradient de par et d'autre du pixel situee en (x,y).
- Si ``|\nabla I(x,y)|`` est localement maximal dans la direction du gradient on garde ce maxima sinon on le supprime 

__4.2 Identifiez les incovenients de cette methode et proposez d'autre alternative plus facile a mettre en oeuvre__.

__4.3 Cette methode est operationnel pour la cartographie des maxima de l'image de gradient. Comment doit-on proceder pour l'image du Laplacien ?__ 
"""

# ╔═╡ 800e166e-5079-451b-b257-a6949063964d
md"### Seuillage de la cartographie des maximas: "

# ╔═╡ 64743737-01c9-40c1-bf10-1367b3e12aa9
md"""
##### Seuillage simple :

Seuil = $(@bind seuil Slider(1:50,default=25,show_value=true))
"""

# ╔═╡ 41623cc9-e642-487b-8d25-e0d04539d102
md"""
##### Seuillage double:

- Dans ce type de seuillage on specifie deux seuils l'un elevee et l'autre relativement bas et on construit deux nouvelles images de gradient: 

```math
\begin{align*}
|\mathbf{\nabla}I(x,y)|_H &= |\mathbf{\nabla}I(x,y)| > T_H \\
|\mathbf{\nabla}I(x,y)|_L &= |\mathbf{\nabla}I(x,y)| > T_L
\end{align*}
```
- On met a jours ``|\mathbf{\nabla}I(x,y)|_L`` en lui soustrayant les valeurs des pixels de ``|\mathbf{\nabla}I(x,y)|_H``.

- On visite le pixel p(x,y) de ``|\mathbf{\nabla}I(x,y)|_H`` et on extrait son voisnage V₈.

- On ne s'interesse qu'aux pixels ayant une valeur nulle dans ce voisinage. Si leurs homologues dans ``|\mathbf{\nabla}I(x,y)|_L`` sont a 1 on les mets aussi a 1 dans ``\mathbf{\nabla}I(x,y)|_H``.

- Cette operation est repetee pour chaque ``p \in \mathbf{\nabla}I(x,y)|_H``. Une fois terminee la matrice  ``\mathbf{\nabla}I(x,y)|_H`` contient les vrai pixel du contour de l'image.
"""

# ╔═╡ c9594400-ffa2-4f20-a969-0b1f77e06609
md""" 
### Post-traitement de l'image de contour

Cette phase consiste a appliquer divers operations:pour :

- Chainner les contours discontinus;
- Amincir les contours epais
- supprimer les points et les segments orphelins.
- Coder les contours en entites geometriques exploitables.
"""

# ╔═╡ defa321d-19a9-4855-aa41-3cd55c05a846


# ╔═╡ 711ad271-a520-425b-8f5e-fa4bfc1a4256
md"""
##### Filtrage optimal 
__1. Filtre de Canny__
"""

# ╔═╡ 1b5c36e7-55eb-431f-b4b7-1170a0b324c8
function dgauss(σ)
	xlim = ceil(Int,5σ)
	x= -xlim:xlim
	y=copy(x)
	kern = Matrix{Float64}(undef,length(x), length(y))
	for (i,vx) in enumerate(x)
		for (j,vy) in enumerate(y)
			kern[i,j] = -vx/(2π*σ^4)*exp(-(vx^2+vy^2)/(2σ^2))
		end
	end
	return centered(kern)
end

# ╔═╡ 3c64d218-fb3c-4571-86b7-6623ed75b987
let
	img_gn = view(img,:,:) .+ σ * randn(MersenneTwister(1234), size(img))
	img_dx = imfilter(view(img_gn,:,:),dgauss(0.75),"replicate")
	img_dy = imfilter(view(img_gn,:,:),dgauss(0.75)',"replicate")
	#α = atan.(img_dy,img_dx)
	#[img 10(img_dx .* cos.(α) + img_dy .* sin.(α))]
	[10abs.(img_dx) img_gn 10abs.(img_dy)]
end


# ╔═╡ 08ca2ac8-f5e2-44a5-b714-06a5029c7b96
let
	img_gn = view(img,:,:) .+ σ * randn(MersenneTwister(1234), size(img))
	img_dx = imfilter(view(img_gn,:,:),dgauss(0.75),"replicate")
	img_dy = imfilter(view(img_gn,:,:),dgauss(0.75)',"replicate")
	α = atan.(img_dy,img_dx)
	[img 10(img_dx .* cos.(α) + img_dy .* sin.(α))]
	#[10abs.(img_dx) img_gn 10abs.(img_dy)]
end

# ╔═╡ 43fd6e4c-a8c7-4fc8-bab8-42a24352bd4e
let
	img_gn = view(img,:,:) .+ σ * randn(MersenneTwister(1234), size(img))
	img_dx = imfilter(view(img_gn,:,:),dgauss(0.75),"replicate")
	img_dy = imfilter(view(img_gn,:,:),dgauss(0.75)',"replicate")
	#α = atan.(img_dy,img_dx)
	#[img 10(img_dx .* cos.(α) + img_dy .* sin.(α))]
	[img_gn 10(abs.(img_dx) + abs.(img_dy))]
end

# ╔═╡ abe395d0-04d9-45f8-a974-8b540c00ba66
let
	img_gn = view(img,:,:) .+ σ * randn(MersenneTwister(1234), size(img))
	img_dx = imfilter(view(img_gn,:,:),dgauss(0.75),"replicate")
	img_dy = imfilter(view(img_gn,:,:),dgauss(0.75)',"replicate")
	[img_gn 10sqrt.(img_dx.^2 + img_dy.^2)]
end

# ╔═╡ 97b62445-ce53-4b56-a0f2-243cac4bf04f
let
	img_gn = view(img,:,:) .+ σ * randn(MersenneTwister(1234), size(img))
	img_dx = imfilter(view(img_gn,:,:),dgauss(0.75),"replicate")
	img_dy = imfilter(view(img_gn,:,:),dgauss(0.75)',"replicate")
	#α = atan.(img_dy,img_dx)
	#[img 10(img_dx .* cos.(α) + img_dy .* sin.(α))]
	# [img_gn 10(abs.(img_dx) + abs.(img_dy))]
	img_grad = view(img_dx,:,:) .^2 + view(img_dy,:,:).^2 .|> sqrt
	p1 = plot(imhist(img_grad,256,0,1)[2], label="")
	vline!(p1,[seuil])
	p2 = heatmap(img_grad .> seuil/255, framestyle=:none, yflip=true, c=:grays)
	#Gray.(img_grad .> 25/255)
	plot(p1,p2, legend=:none, aspect_ratio=:none, colorbar=false)
end

# ╔═╡ Cell order:
# ╟─437aabde-e7f5-4ed4-8f07-b91ecf5221ba
# ╟─10028c35-cf65-43d5-b8f8-4a160416fa5b
# ╟─59ef7fc4-4f90-4154-8dd8-0260f25ace12
# ╟─90343e63-6739-421d-99e0-bf3bc0289fa5
# ╟─38bf4f37-d915-4375-90eb-3de73fc9257a
# ╟─3d15739d-f368-4717-8b60-3b2353bf1bfe
# ╟─93f35d91-6dbc-4a12-a176-d397f335fba6
# ╟─57bd4f49-44fa-48a5-b0c1-44b42d42af58
# ╟─e9c02c13-630b-4991-89cc-1ef91b1ef82b
# ╟─5473a579-e4b0-4dbc-b277-9096bb992d4c
# ╟─966a5176-c131-4dc0-85e8-faa9bc79e327
# ╟─7bc3bd10-06df-49df-bcf0-9c72d4da6ef2
# ╟─4ff2504e-07b4-4cfc-a4ad-8454989a78d1
# ╟─6b6e5883-0af1-4c69-95a5-0d521ee5bf05
# ╟─4a8226e4-1237-434e-be0f-a51ce55b5687
# ╟─b550f65e-f97f-4aba-a48f-02694628c304
# ╟─3ece7e97-28e7-460d-84f4-9d1905585271
# ╟─43577282-bd3b-4d1b-92e9-dea1dd5a853c
# ╟─1a4a993c-fa45-47de-a400-973f53f2ca5c
# ╟─d82645ff-31e7-43a0-84ee-4dd612e40972
# ╟─ab64a060-8f14-44bd-8604-580575e1f31b
# ╟─3c64d218-fb3c-4571-86b7-6623ed75b987
# ╟─3220f95e-74a9-4578-97ba-0e1123731ef3
# ╟─5d735d83-0556-470b-920d-1a98c7bf5c7f
# ╟─55bf3ae5-d514-4af3-a269-0262c2d74f0a
# ╟─bebb8072-33c7-468a-a9e9-828c028d1294
# ╟─08ca2ac8-f5e2-44a5-b714-06a5029c7b96
# ╟─877552e7-6ac5-4bd2-9b8d-bd2a4b6a87b6
# ╟─43fd6e4c-a8c7-4fc8-bab8-42a24352bd4e
# ╟─8f4b9bcf-09c4-4750-aec5-337d6c8a065a
# ╟─abe395d0-04d9-45f8-a974-8b540c00ba66
# ╟─709d40c3-d185-4e06-9fc8-e61dc6a6a680
# ╟─a4dde3e2-b13f-4a9c-8fbd-a83565cff09d
# ╟─d31b3819-e2a3-4fec-bd73-1f6334990228
# ╟─800e166e-5079-451b-b257-a6949063964d
# ╟─64743737-01c9-40c1-bf10-1367b3e12aa9
# ╟─97b62445-ce53-4b56-a0f2-243cac4bf04f
# ╟─41623cc9-e642-487b-8d25-e0d04539d102
# ╟─c9594400-ffa2-4f20-a969-0b1f77e06609
# ╠═defa321d-19a9-4855-aa41-3cd55c05a846
# ╟─711ad271-a520-425b-8f5e-fa4bfc1a4256
# ╟─1b5c36e7-55eb-431f-b4b7-1170a0b324c8
