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

# ╔═╡ 962e4274-38fd-11ed-3890-bfb3c8d0f09d
begin
	import Pkg
	Pkg.activate(pwd())
end

# ╔═╡ 2fb557e2-258f-4d47-81a3-495ef846b080
begin
	using Images, TestImages, TiffImages
	using Plots, LaTeXStrings, Statistics
	using PlutoUI
end

# ╔═╡ b4a6ae81-bd2c-4245-b22b-47c5c317d417
html"""
<div style="margin:10%;background-color:#0303030e;">
<h1 style="text-align:center;color:darkblue;"> IMAGE PROCESSING</h1>
<h2 style="text-align:center;color:darkorange;"> Part one</h2>
</div>
"""

# ╔═╡ 1cbf1b09-c594-4415-8301-7641944cc599
md"""
## Introduction

"""

# ╔═╡ 3308305f-a842-4dd4-b03a-500bd7fa4f8b
md"""
Le prétraitement des images est l'ensemble des opérations qui visent à préparer l'image aux phases d'analyse de niveau supérieur. Nous citons à titre d'exemple sans être exhaustif les opérations ci-après :

+ La correction radiométrique et chromatique de l'image.
+ La réduction du bruit
+ L'amélioration de la netteté de l'image
+ Le redressement géométrique de l'image
+ ...
Pour exécuter ces opérations, nous serons amenés, dans certains cas, à manipuler les pixels de l'image indépendamment de leur position géométrique et de leur entourage. On parle alors de transformation ponctuelle. Dans d'autres cas, nous changerons l'intensité de chaque pixel en prenant en compte les propriétés de son voisinage et on parle de transformation par voisinage. Et parfois on ne manipulera que grille spatiale de l'image dans le cas du redressement de l'image.
"""

# ╔═╡ 85234a14-9236-41bc-8acb-5967b4630281
md"""
## Operations ponctuelles

"""

# ╔═╡ d2f1705d-f46d-400a-9c12-10da2ae5ca2d
md"""
Soit ``r`` la valeur du pixel de l'image d'origine et ``s`` la valeur du pixel dans l'image modifiée. Dans le contexte des transformations ponctuels ``s`` ne dépend que de ``r`` a travers une opération de mappage linéaire ou non linéaire qu'on appelle transformation ponctuelle ``\mathcal{T}`` tel que :
```math
s = \mathcal{T}(r)
```

Les transformations lineaires couremment utilisees ont la forme suivante:
```math
\mathcal{T}(r) = c \cdot r + b
```
Avec ``c`` le parametre controlant le gain ou le __contraste__ de l'image et ``b`` le bias eou encore la __brillance__ de l'image. 
Le choix des valeurs de ``c`` et de ``b`` va dépendre de la distribution des intensités des pixels de l'image de départ en plus du résultat qu'on cherche à mettre en place.

Pour les transformations non lineaires, elles ont la forme generale suivante:

```math
s = c \cdot r^\gamma \text{ avec } \gamma \in \mathbb{R}
```
"""

# ╔═╡ a7d3fe04-71ab-4237-abd1-5283683c29e0
md"γ = $(@bind γ Slider([2.0^i for i=-3:3],default=1.0, show_value=true))"

# ╔═╡ 544cd5cb-368d-47e5-b59f-77c529a34ca0
let
	r = collect(Float64,0:1/255:1.0)
	plot(r, r.^γ,
			label=L"\gamma="*"$(round(γ, digits=3))",
			legend=:outerright,
			ylabel=L"s",
			xlabel=L"r",
			title =L"\gamma-" * "transformations",)
end

# ╔═╡ f66ed72c-a63c-41ab-aaaf-756c051fa153
md"### Definitions"

# ╔═╡ 749872e6-e689-4cf2-8a48-95db37e21e38
md"""
#### Histogramme
L'histogramme de l'image est un tableau de dimension ``2^n`` valeurs (avec ``n`` le nombre de bits/canal de couleur) contenant les fréquences d'occurrence des tons de couleurs/niveau de gris de l'image étudiée. 
"""

# ╔═╡ 1de1e52b-98a9-4ca1-9397-b7421d646f4e
begin
	img = testimage("barbara_color.png");
	canals = [reshape(channelview(img)[i,:,:], size(img)) for i=1:3];
	p1 = heatmap(img);
	colors = [:reds, :greens, :blues]
	plot(p1,[heatmap(canals[i], color= colors[i]) for i=1:3]...,yflip=true, framestyle=:none, layout=4, aspect_ratio=:equal,colorbar=:none)
end

# ╔═╡ d72c5e3d-ecec-4e94-912d-43b7ccbffa0a
let
	colors = [:red,:green,:blue];
	canals = [red,green,blue];
	titles = ["Reds", "Greens", "Blues"]
	plot(heatmap(img,framestyle=:none),[histogram(Int.(c[1].(img)[:]*255),bins=0:255,alpha=0.35,color=c[2], title=c[3],titlefontsize=8,titlefontcolor=c[2]) for c in zip(canals,colors,titles)]...,layout=4,grid=:none, legend=:none)
end

# ╔═╡ 1d68bc2d-42a2-49d2-9278-d42a5984f7ad
md"""
### Histogramme normalisé
```math
H_{norm} = H/(M\times N) \text{ avec } \left\{\begin{matrix}M : & nombre de lignes \\
N : & Nombre de colonnes \\
\end{matrix}\right.
```
"""


# ╔═╡ 8d45b162-7c5d-422e-9234-730635b4dead
let
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/bark_he_512.tiff") |> TiffImages.load
	plot(heatmap(img,framestyle=:none),histogram(Int.(img*255)[:],color=:gray, title="Grayscale",titlefontsize=8,titlefontcolor=:gray,bins=0:255,normalize=:pdf),layout=2,grid=:none,aspect_ratio=:none,legend=:none)
end

# ╔═╡ 25eaf2af-082c-405f-9301-a2a680508a51
md" ### Histogramme cumulé"

# ╔═╡ 43d4153f-66aa-46b3-ab1b-fc9aeaad957d
md"""
L'histogramme cumulé est calculé par la relation suivante :
```math
H_C[k] = \sum_{i=1}^k H_{norm }[i]
```
"""


# ╔═╡ 4cdbdf8f-a32c-48ca-971d-ba7ef066bcd5
let
	img = (download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/livingroom.tif") |> TiffImages.load)
	edges, H = imhist(Int.(255*img[:]),0:254)
	Hnorm = H / reduce(*,size(img))
	Hc = cumsum(Hnorm)
	p = bar(0:255,Hc, bar_width=0.8, color=:gray, alpha=0.35)
	title!("Histogramme cumule", titlefontsize=8, titlefontcolor=:gray, legend=:none)
	plot(heatmap(img,framestyle=:none), p, aspect_ratio=:none)
end

# ╔═╡ f9ee9d95-82ed-4a8c-a85a-12dc32544c05
md"""
### Brillance moyenne
C'est la moyenne arithmetique des intensites de l'image en niveau de gris. En generale, pour une image RGB on utilise la formule suivante pour obtenire l'image de distribution de luminiscence ``L`` (image en niveau de gris):
```math
L = 0.299 \times R + 0.587 \times G + 0.114 \times B
```
"""

# ╔═╡ e0198a08-1c7c-4176-b71c-6cea89621f8d
let 
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/lighthouse.png") |> load;
	mosaicview(img, Gray.(img), nrow=1)
end

# ╔═╡ 5e353c28-8882-4ff2-a2b7-1861e88b5176
md"""
La brillance moyenne de cette image est donnee par:
```math
I_{moy} = \frac{1}{M\times N}\sum_{m,n} I[m,n]
```
+ Exemple :
"""

# ╔═╡ 94e6a458-72aa-4181-8e6a-5b70449b75af
let
	img = Gray.(testimage("lighthouse.png"));
	Imoy = floor(Int,Float64(mean(img[:]))*255)
	p1 = heatmap(img, framestyle=:none);
	p2 = histogram(Int.(255*img[:]),bins=0:254, normalize=:pdf, label=L"H_{norm}", alpha=0.35)
	p2 = vline!([Imoy], lw=3.0, label=L"I_{moy}="*"$Imoy")
	plot(p1, p2,aspect_ratio=:none, grid=:none)
end

# ╔═╡ 7518e615-8a8e-48df-ae1c-8382468d3df1
md"""
### Contraste et variance

#### Contraste:

```math
C = \frac{I_{max} - I_{min}}{I_{max}+I_{min}}
```

#### Variance 

```math
\sigma^2 = \frac{1}{M \times N} \sum_{m,n} \left(I[m,n]-I_{moy}\right)^2
```
"""

# ╔═╡ 101f6ffe-d855-416c-9566-0e1217f6c8c8
md"### Exemple de transformation lineaire ponctuelle"

# ╔═╡ ee278d89-f93a-4b32-b0dc-3d4a550aab43
md"""
#### Inversion d'histogramme
```math
s = 255 - r
```
"""

# ╔═╡ 2c714091-b91a-42fb-ae8c-10950ab3df84
let 
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/house.tif") |> TiffImages.load;
	img = Gray.(channelview(img)[1,:,:]);
	_, H = imhist(Int.(img*255),0:254)
	H /=reduce(*,size(img))
	_, H_inv = imhist(255 .- Int.(img*255),0:254)
	plot(bar(0:255, H, label="original"),
		bar(0:255, H_inv, label="inversed"),
		heatmap(img, framestyle=:none),
		heatmap(1N0f8 .- img, framestyle=:none), layout=4)
end

# ╔═╡ d1e46ce1-0175-49a5-86e3-6daba62a3d58
md"""
#### Compression d'histogramme

```math
\begin{align*}
r_{norm} &= \frac{r - r_{min}}{r_{max} - r_{min}}\\
s &= \frac{r_{norm}}{MAX - MIN} + MIN\\
\end{align*}
```
"""

# ╔═╡ 86191fe0-6f8b-430c-a37a-3688c025af44
let
	img = Gray.(testimage("lighthouse.png"))
	p1 = histogram(Int.(255img[:]), bins=0:255, label="original");
	MAX, MIN  = 165/255, 55/255;
	img_comp = Float64.(img .- minimum(img))/Float64(maximum(img) -	minimum(img)) * (MAX-MIN) .+ MIN;
	p2 = histogram(round.(Int,255img_comp[:]), bins=0:255, label="compressed")
	p3 = heatmap(img, framestyle = :none)
	p4 = heatmap(Gray.(img_comp), framestyle = :none)
	plot(p1,p2,p3,p4,layout=4)
end

# ╔═╡ 8e21d673-ccad-4241-8cf1-057ffd4d398d
md"### Exemple de transformation non lineaire"

# ╔═╡ ee4901f6-7729-4530-9b4a-0f4253bd702e
md"""
#### γ-Transformation

γ = $(@bind γ1 Slider([2.0^i for i=-3:3],default=1.0, show_value=true))
	
"""

# ╔═╡ d1913389-a3da-47fc-a145-b39b4446326f
let
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/jetplane.tif") |> TiffImages.load .|> Gray
	p1 = histogram(Int.(255img[:]), bins=0:255, label="original");
	img_nl_trans = img.^(γ1)
	p2 = histogram(floor.(Int,Float32.(255*img_nl_trans[:])), bins=0:255, label=L"\gamma\, transformation")
	p3 = heatmap(img, framestyle = :none)
	p4 = heatmap(img_nl_trans, framestyle = :none)
	plot(p1,p2,p3,p4,layout=4, legend=:outertop, grid=false)
end

# ╔═╡ 627ab2d7-2546-4507-bb5a-ae0706185803
md"""
#### Transformation statistique

##### Égalisation d'histogramme

```math
s = \mathcal{T}(r) = (L-1)\int_0^r p_r(\omega)d\omega 
```
"""

# ╔═╡ e8cb0843-c626-4be5-852f-e16e2390cc46
let
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/moonsurface.tiff") |> TiffImages.load
	_, H = build_histogram(img,256,0,1)
	Hc = cumsum(H/reduce(*,size(img)))
	im_eq = Gray{N0f8}.(Hc[Int.(255*img)])
	_, Heq = build_histogram(im_eq,256,0,1)
	#for i in 1:size(img,1)
	#	for j in 1:size(img,2)
	#		im_eq[i,j] = N0f8(Hc[Int(255*img[i,j])])
	#	end
	#end
	l= @layout [a b ; c]
	plot(bar(0:256,H.parent, alpha=0.35,label="original"),bar(0:256,Heq.parent, alpha=0.35, label="Equalized"),heatmap(mosaicview(img, im_eq, nrow=1), framestyle=:none),layout=l, aspect_ratio=:none, legend=:outertop)
end

# ╔═╡ 2175d9e5-7744-406a-a8da-933a5fce11c5
md"""
##### Exercise
__Imagine and code an algorithm that perform adaptative histogram equalization on grayscale images.__
"""

# ╔═╡ ab6af097-1692-42ca-abca-e7af9c965617
md"#### Egalisation pour image couleur"

# ╔═╡ 9acd7b71-a53c-485c-8b5c-9ee6b5b93783
md"##### Can you describe what we must do in the case of RGB image ?"

# ╔═╡ 6a54929d-f965-4c28-a926-8fc85c95fc57
let
	hsv = testimage("monarch_color.png") .|> HSV
	hsv_eq = Vector{HSV}()
	V= Vector{Float32}(undef,length(hsv))
	for i =1:length(V)
		V[i] = (@view hsv[:])[i].v
	end
	V = reshape(V,size(hsv))
	Veq = adjust_histogram(V,Equalization(nbins = 256, minval = 0, maxval = 1))
	for item in  zip(hsv, Veq)
		push!(hsv_eq,HSV(item[1].h, item[1].s, item[2]))
	end
	[hsv reshape(hsv_eq,size(hsv))]
end

# ╔═╡ Cell order:
# ╟─962e4274-38fd-11ed-3890-bfb3c8d0f09d
# ╟─2fb557e2-258f-4d47-81a3-495ef846b080
# ╟─b4a6ae81-bd2c-4245-b22b-47c5c317d417
# ╟─1cbf1b09-c594-4415-8301-7641944cc599
# ╟─3308305f-a842-4dd4-b03a-500bd7fa4f8b
# ╟─85234a14-9236-41bc-8acb-5967b4630281
# ╟─d2f1705d-f46d-400a-9c12-10da2ae5ca2d
# ╟─a7d3fe04-71ab-4237-abd1-5283683c29e0
# ╟─544cd5cb-368d-47e5-b59f-77c529a34ca0
# ╟─f66ed72c-a63c-41ab-aaaf-756c051fa153
# ╟─749872e6-e689-4cf2-8a48-95db37e21e38
# ╟─1de1e52b-98a9-4ca1-9397-b7421d646f4e
# ╟─d72c5e3d-ecec-4e94-912d-43b7ccbffa0a
# ╟─1d68bc2d-42a2-49d2-9278-d42a5984f7ad
# ╟─8d45b162-7c5d-422e-9234-730635b4dead
# ╟─25eaf2af-082c-405f-9301-a2a680508a51
# ╟─43d4153f-66aa-46b3-ab1b-fc9aeaad957d
# ╟─4cdbdf8f-a32c-48ca-971d-ba7ef066bcd5
# ╟─f9ee9d95-82ed-4a8c-a85a-12dc32544c05
# ╟─e0198a08-1c7c-4176-b71c-6cea89621f8d
# ╟─5e353c28-8882-4ff2-a2b7-1861e88b5176
# ╟─94e6a458-72aa-4181-8e6a-5b70449b75af
# ╟─7518e615-8a8e-48df-ae1c-8382468d3df1
# ╟─101f6ffe-d855-416c-9566-0e1217f6c8c8
# ╟─ee278d89-f93a-4b32-b0dc-3d4a550aab43
# ╟─2c714091-b91a-42fb-ae8c-10950ab3df84
# ╟─d1e46ce1-0175-49a5-86e3-6daba62a3d58
# ╟─86191fe0-6f8b-430c-a37a-3688c025af44
# ╟─8e21d673-ccad-4241-8cf1-057ffd4d398d
# ╟─ee4901f6-7729-4530-9b4a-0f4253bd702e
# ╟─d1913389-a3da-47fc-a145-b39b4446326f
# ╟─627ab2d7-2546-4507-bb5a-ae0706185803
# ╟─e8cb0843-c626-4be5-852f-e16e2390cc46
# ╟─2175d9e5-7744-406a-a8da-933a5fce11c5
# ╟─ab6af097-1692-42ca-abca-e7af9c965617
# ╟─9acd7b71-a53c-485c-8b5c-9ee6b5b93783
# ╟─6a54929d-f965-4c28-a926-8fc85c95fc57
