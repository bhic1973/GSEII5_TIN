### A Pluto.jl notebook ###
# v0.19.40

#> [frontmatter]
#> "Filière" = "GTR S4"
#> title = "Image Processing Part 2"
#> date = "2024-03-29"
#> tags = ["Pr. Hicham Belkebir"]
#> 
#>     [[frontmatter.author]]
#>     name = "Hicham Belkebir"

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
	Pkg.activate(".")
end

# ╔═╡ 2fb557e2-258f-4d47-81a3-495ef846b080
begin
	using Images, TestImages, TiffImages
	using Plots, LaTeXStrings, Statistics
	using PlutoUI
end

# ╔═╡ 0e2fef44-65e4-45bd-a126-1c15c598121c
let
	using FFTW
	img = Gray.(testimage("lighthouse.png") )
	IMG_x = Matrix{ComplexF32}(undef,size(img))
	IMG = copy(IMG_x)
	px = plan_fft(Float32.(img[1,1:end]);flags=FFTW.ESTIMATE, timelimit=Inf)
	for i in 1:size(img,1)
		IMG_x[i,:] = px*Float32.(img[i,:])
	end
	py = plan_fft(IMG_x[1:end,1];flags=FFTW.ESTIMATE, timelimit=Inf)
	for j in 1:size(img,2)
		IMG[:,j] = py*IMG_x[:,j]
	end
	IMG_abs_dB = 20log10.(fftshift(abs.(IMG)/maximum(abs.(IMG))))
	p1 = heatmap(IMG_abs_dB, color=:jet)
	p2 =heatmap(img)
	p3 = heatmap(fftshift(angle.(IMG)),color=:grays)
	plot(p2,p1,p3,aspect_ratio=:equal, colorbar=:none, framestyle=:none, yflip=true,layout=@layout [a ;b c])
end

# ╔═╡ de935179-fd70-4a68-8f7f-e23b2b9187aa
# ╠═╡ disabled = true
#=╠═╡
begin
	using OffsetArrays
	function im_edge_extent(im,K)
		M,N = size(im)
		im_ext = OffsetMatrix(Matrix{eltype(im)}(undef,M+2K+1,N+2K+1),-K:M+K,-K:N+K)
		im_ext[1:M,1:N] = @view im[:,:]
		im_ext[-K:0,1:N] = repeat(im[begin,:]',K+1,1)
		im_ext[M+1:M+K,1:N] = repeat(im[end,:]',K,1)
		im_ext[:,-K:0] = repeat(im_ext[:,1],1,K+1)
		im_ext[:,N+1:N+K] = repeat(im_ext[:,N],1,K)
		return Float64.(im_ext)
	end

	function blf(im, K, σs, σr)
		im_filt = copy(im)
		im_ext = im_edge_extent(im,K)
		blf_space_x = ones(Int,length(-K:K)) * collect(-K:K)';
		blf_space_y = blf_space_x'
		blf_space = exp.(-(blf_space_x.^2 + blf_space_y.^2)/(2σs^2));
		M,N = size(im)
		#blf_space = (Kernel.gaussian(σs)).parent
		for i=1:M
			for j=1:N
				V = @view im_ext[i-K:i+K, j-K:j+K]
				blf_range = exp.(-(im_ext[i,j] .- V ).^2/2/σr^2) 
				blf_w = blf_range .* blf_space
				blf_w /= sum(blf_w);
				im_filt[i,j] = sum(V .* blf_w) > 1.0 ?  Gray(1N0f8) : Gray{N0f8}(sum(V .* blf_w));
			end
		end
		return im_filt
	end
end
  ╠═╡ =#

# ╔═╡ b4a6ae81-bd2c-4245-b22b-47c5c317d417
html"""
<div style="margin:10%;background-color:#0303030e;">
<h1 style="text-align:center;color:darkblue;"> IMAGE PROCESSING</h1>
<h2 style="text-align:center;color:darkorange;"> Part two</h2>
</div>
"""

# ╔═╡ 9edd311a-b2a2-46f8-89d2-e21a71edcaef
md"""
__Regardons cette video explicative :__
$(@bind play CheckBox(default=false))
"""

# ╔═╡ f1b20ca8-2e71-49a4-98a0-da4d7ac9ed1a
let
	using PyCall
	ipd = pyimport("IPython.display")
	if play
		ipd.Video("cv_pics/Conv_demo.mp4", embed=true)
	end
end

# ╔═╡ 1cbf1b09-c594-4415-8301-7641944cc599
md"""
## Filtrage lineaire

"""

# ╔═╡ 63be3bce-7a1d-41b7-918d-00cccbc575df
md"### Produit de convolution"

# ╔═╡ 3308305f-a842-4dd4-b03a-500bd7fa4f8b
md"""
C'est une transformation qu'on applique aux pixels de l'image en se basant cette fois-ci sur l'information du voisinage spatiale du pixel ``I[i,j]`` considéré:
```math
I_{dest}[i,j] = \mathcal{T}(V(I_{orig}[i,j]),I_{orig}[i,j])
```
L'opérateur ``\mathcal{T}`` doit répondre à deux caractéristiques pour que son action puisse prendre place :
+ Il doit être linéaire
+ Il doit être spatialement invariant

Une fois ceci est vérifié, la transformation ``\mathcal{T}`` peut être effectuée en utilisant le produit de convolution :
```math
I_{dest}[i,j] = h*I = \sum_{m,n}h[m,n]I[i-m,j-n]
```
Avec ``h`` est le noyau (kernel) représentant la transformation ``\mathcal{T}``. ``h`` est appelé aussi *__fonction d'étalement du point__* (PSF).
"""

# ╔═╡ 9eb8287e-b81b-48c7-8b97-72cdcdc86f20
md"""
### Implementation 1D du produit de convolution 2D

Considerons ``K`` un kernel de convolution et ``\mathbf{K}`` sa representation matricielle. Procedons a la decomposition en valeur singuliere de ``\mathbf{K}``. Il en resulte :
```math
\mathbf{K} = \sum_i \sigma_i \mathbf{u}_i \mathbf{v}_i^T  
```
On peux pratiquer une separation du kernel 2D en deux kernels 1D a actions orthogonales si ``\sigma_0 \ne 0``.

```math
h[m,n] = h_x[m] * h_y[n] \text{ avec } h_x = \sqrt{\sigma_0} \mathbf{u_0}\text{ et } h_y = \sqrt{\sigma_0} \mathbf{v}_0 
```
On commence alors par appliquer ``h_x[n]`` au niveau de chaque ligne de l'image d'origine :

```math
I_x[i,j] = \sum_{n=-K}^K h_x[n] I[i, j-n]\;\forall\; 1\leq i \leq M
```
Ensuite, on applique ``h_y[m]`` dans la direction orthogonale sur l'image ``I_x`` :

```math
I_{xy}[i,j] = \sum_{m=-K}^K h_y[m] I_x[i-m,j]\; \forall\; 1\leq j \leq N 
```

Il devient alors envisageable d'implementer ces operations sous forme de produit entre matrice creuse ``\mathbf{H}`` et vecteur concu a partir des lignes ou colonnes de l'image :

```math
\begin{align*}
I_x &= \mathbf{H}\,I[i,1:N]\;\forall\; i\\
I_{xy} &=\mathbf{H}\,I_x[1:M,j]\;\forall\;j\\
\end{align*}
```
"""

# ╔═╡ 98203216-6988-4e62-b1d5-eaf355cb6c1a
md"""
### Probleme des bords de l'image

En produit de convolution le traitement des pixels de l'interieur de l'image differe de ceux qui se trouvent sur les bandes de bord a cause de la nature finie de l'image. Pour palier a ce probleme on peut proceder de quatre facons differentes:
+ Considerer que l'image s'etend au dela de ses bords naturelle en attribuant aux pixels de cette zone une valeur constante ( zero ou autre):
```math
I[i,j]=0\;\forall\; i < 1 \textbf{ or } i>M  \textbf{ or } j < 1 \textbf{ or } j > N
```
+ Etendre la grille de l'image au dela de ses limites naturelles en repliquant ses bords et ses coins autant de fois qu necessaire:
```math
\begin{align*}
I[i,j] &= I[i,1] &\forall\, j \in[-K, 0]\;\forall\, 1\leq i\leq M \\
I[i,j] &= I[i,N] &\forall\, j \in[N+1, N+K]\;\forall\, 1\leq i\leq M \\
I[i,j] &= I[1,j] &\forall\, i \in[-K, 0]\;\forall\, 1\leq j\leq N \\
I[i,j] &= I[M,j] &\forall\, i \in[M+1, M+K]\;\forall\, 1\leq j\leq N \\
\end{align*}
```
+ Etendre l'image en considerons qu'elle est la realisation d'un procede periodique
```math
\begin{align*}
I[i,j] &= I[i,N+j] &\forall\, j \in[-K, 0]\;\forall\, 1\leq i\leq M \textbf{ and } j>N \\
I[i,j] &= I[i,] &\forall\, j \in[N+1, N+K]\;\forall\, 1\leq i\leq M \\
I[i,j] &= I[1,j] &\forall\, i \in[-K, 0]\;\forall\, 1\leq j\leq N \\
I[i,j] &= I[M,j] &\forall\, i \in[M+1, M+K]\;\forall\, 1\leq j\leq N \\
\end{align*}
```
+ Appliquer une symetrie miroir autour des bords de l'image :
```math
\begin{align*}
I[i,k] &= I[i,-k] & \textbf{ for } j=1 \textbf{ and } \forall k \in [-K, 0]\; \forall i\\ 
I[i,k] &= I[i,N-k] & \textbf{ for } j=N \textbf{ and } \forall k \in [1, K]\; \forall i\\ 
I[i,k] &= I[-k,j] & \textbf{ for } i=1 \textbf{ and } \forall k \in [-K, 0]\; \forall j\\ 
I[i,k] &= I[M-k,j] & \textbf{ for } i=M \textbf{ and } \forall k \in [1, K]\; \forall j\\ 
\end{align*}
```
"""

# ╔═╡ bf3c567d-3cc5-472c-a4c0-bb98bebfa1d6
md"""
### Exemples de filtrage linéaire

#### Moyenne mobile ou filtre à boitte (Lissage / Bluring )
"""

# ╔═╡ caf20bd3-8a7a-4025-b3df-76f8c9f61458
md"Kernel dim = $(@bind dim Slider(3:2:21,default=3, show_value=true))"

# ╔═╡ ccb00a2d-6c86-4490-af16-b47dd66ff757
md"""
Border = $(@bind border Select(["replicate", "circular", "reflect", "symmetric"]))
"""

# ╔═╡ ed476ab1-ab27-4907-a305-78a9eed21213
let
	img = testimage("chelsea.png")
	K = centered(ones(dim,dim))
	K /=dim^2
	mosaicview(Gray.(img),imfilter(Gray.(img), K, border),nrow=1)
end

# ╔═╡ fbad8364-e23c-4303-90b8-d71181cecdc7
md"""
#### Filtre de Bartlett (filtre bilineaire)
```math
\mathbf{K}_R = \mathbf{K}_B = {1\over4}\begin{pmatrix}
1 & 2 & 1\\
2 & 4 & 2\\
1 & 2 & 1\\
\end{pmatrix}\quad \mathbf{K}_G = {1\over4}\begin{pmatrix}
0 & 2 & 0\\
2 & 4 & 2\\
0 & 2 & 0\\
\end{pmatrix}

```
"""

# ╔═╡ 104cdab1-fd83-4ae4-bd53-2d685156c6c4
md"Zoom on bayer image = $(@bind zoom CheckBox(default=true))"

# ╔═╡ 580d6634-9b7b-4387-9081-ad7e9d6b33be
let
	img = testimage("lighthouse.png")
	KR = KB = 0.25*centered([1 2 1;2 4 2;1 2 1])
	KG = 0.25 * centered([0 1 0; 1 4 1; 0 1 0])
	G = zeros(N0f8,size(img))
	G[1:2:end,2:2:end] = (green.(img))[1:2:end,2:2:end]
	G[2:2:end,1:2:end] = (green.(img))[2:2:end,1:2:end]
	R = zeros(N0f8,size(img))
	B = copy(R)
	B[2:2:end,2:2:end] = (blue.(img))[2:2:end,2:2:end]
	R[1:2:end,1:2:end] = (red.(img))[1:2:end,1:2:end]
	if zoom
		[Gray.(R)[1:5,1:5] Gray.(G)[1:5,1:5] Gray.(B)[1:5,1:5]]
	else
		Ri = imfilter(R, KR, border)
		Gi = imfilter(G, KG, border)
		Bi = imfilter(B, KB, border)
		mosaicview(img,colorview(RGB, R, G, B), colorview(RGB, Ri, Gi, Bi),nrow=1)
	end
end

# ╔═╡ 0485cdd7-89d2-45b0-a0da-e62fb1a78cf6
md"""
#### Filtre de lissage gaussien
```math
h(x,y) = {1\over 2\pi \sigma^2} \exp(-{(x^2+y^2)\over2\sigma^2})
```
"""

# ╔═╡ c437baad-cc24-418a-b176-e295e8a812a6
md"σ = $(@bind σ Slider([0.5,1,1.5,2,2.5,3], default=2, show_value=true))"

# ╔═╡ 2200a9aa-34e1-4844-9f92-1370906b4422
let
	img = HSV.(testimage("monarch_color_256.png"))
	V = imfilter(channelview(img)[3,:,:],Kernel.gaussian(σ),border)
	new_img = HSV.(channelview(img)[1,:,:], channelview(img)[2,:,:], V);
	mosaicview(img,new_img, nrow=1)
end

# ╔═╡ f64d7431-9edb-4872-aa73-181408c9c497
md""" 
#### Accentuation de la netteté de l'image

```math
I_{sharp} = I_{orig} + \gamma(I_{orig} − h_{blur} ∗ I_{orig} )
```
"""

# ╔═╡ d5e6c86f-a341-46a9-9658-e18dda02a52c
md"γ = $(@bind γ Slider(0.1:0.1:1.0, default=0.5, show_value=true))"

# ╔═╡ 7b05dc85-982c-4a94-adb8-e64e77464e12
let
	img = HSV.(testimage("monarch_color_256.png"))
	V = (1+γ) * channelview(img)[3,:,:] .-  γ * imfilter(channelview(img)[3,:,:],Kernel.gaussian(σ),border)
	new_img = HSV.(channelview(img)[1,:,:], channelview(img)[2,:,:], V);
	mosaicview(img,new_img, nrow=1)
end

# ╔═╡ 958c9316-bf30-4069-8cb0-5f63f4bfa9de
md"""
#### Filtre de Sobel
```math
\mathbf{H}_{x,sobel} = \begin{pmatrix}
-1 & 0 & 1\\
-2 & 0 & 2\\
-1 & 0 & 1\\
\end{pmatrix} \quad \mathbf{H}_{y,sobel} = \begin{pmatrix}
-1 & -2 & -1\\
0 & 0 & 0\\
1 & 2 & 1\\
\end{pmatrix}
```
"""

# ╔═╡ bcc1e8c0-824b-469f-8ef8-497151054a02
let
	img = testimage("toucan.png")
	sobel_x, sobel_y = 8 .* Kernel.sobel()
	[imfilter(Gray{Float64}.(img), sobel_x) img imfilter(Gray{Float64}.(img), sobel_y)]
end

# ╔═╡ 116f8969-a810-4c87-a6d2-69b1a510f985
md""" 
#### Filtre LoG
```math
LoG(x,y) = (\frac{x^2+y^2}{\sigma^4} -{2\over\sigma^2})G(x,y,\sigma)
```
"""

# ╔═╡ abcea042-2cde-4413-8443-7df5cfaa58ce
let
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/sudoku.tiff") |> TiffImages.load
	[img 10imfilter(Gray.(img), 10Kernel.LoG(σ))]
end

# ╔═╡ f20e0585-7af6-4faa-9da9-903ceee3ca83
md"""
### Filtrage spectrale 

Le filtrage lineaire se resume en une operation de convolution 2D entre un kernel representant l'effet du filtre et l'image sujette de l'operation de filtrage.
Vue la dualite existante entre l'espace directe de l'image et l'espace reciproque (espace construit par transformation soit de Fourier ou autre), on peut aussi concevoir des filtres lineaires agissant dans l'espace reciproque.
"""

# ╔═╡ 8b08079c-565d-4cd2-b48b-42ab4d50dc46
md"""
#### Representation spectrale des images

Cette representation est obtenue par transformee de Fourier discrete 2D selon la formule ci-apres:

```math
\begin{align*}
I[u,v] &= \sum_{m=0}^{M-1}\sum_{n=0}^{N-1} I[m,n] \exp(-\jmath 2\pi m u)\exp(-\jmath 2\pi n v)\\
I[u,v] &= \sum_{m=0}^{M-1}\left(\sum_{n=0}^{N-1} I[m,n] \exp(-\jmath 2\pi n v)\right)\exp(-\jmath 2\pi m u)\\
I[u,v] &= \sum_{m=0}^{M-1}I[m,v]\exp(-\jmath 2\pi m u)
\end{align*}
```
La separabilite de la transformee de Fourier permet de faciliter l'implementation de l'algorithme de calcul en utilisant le fameux algorithme de la transformee de Fourier Rapide (FFT) pour chaque direction de l'espace directe de l'image.

"""

# ╔═╡ 7b5dcbec-68b9-4829-bb8d-2685c0e8ae08
md"""
#### Filtre a boite spectrale:
```math
H[u,v] = \mathbb{1}_{\{u \leq |D_u|, v \leq |D_v|\}}
```
"""

# ╔═╡ da92ef0c-3123-4a30-bea7-e5a65cf412e7
md"``R_B``=$(@bind rB Slider(0.05:0.05:0.5, default=0.1, show_value=true))"

# ╔═╡ b2ef51f4-81c9-4a61-bd04-d35213ae42ed
let
	img = Gray.(testimage("lighthouse.png"))
	HB = zeros(ComplexF32,size(img))
	C = (size(img,1)÷2, size(img,2)÷2)
	display("Rayon du filtre a boite: $rB")
	for (i,u) in enumerate(range(-0.5,0.5,size(img,1)))
		for (j,v) ∈ enumerate(range(-0.5,0.5,size(img,2)))
			HB[i,j] = (abs(u)<=rB && abs(v) <=rB) ? true : false 
		end
	end
	IMG = fft(Float64.(img))
	img_blur = (fftshift(IMG) .* HB) |> ifftshift |> ifft .|> real .|> Gray
	[img Gray.(abs.(HB)) img_blur]
end

# ╔═╡ 0323d736-82a4-4220-bb80-2317cddb65f2
md"""
#### Filtre a disque spectrale:
```math
H[u,v] = \mathbb{1}_{u^2+v^2 \leq D^2}
```
"""

# ╔═╡ 7a293c05-7bee-44db-8e49-3457c92b38fd
md"``R_D``=$(@bind rD Slider(0.05:0.05:0.5, default=0.1, show_value=true))"

# ╔═╡ 11114219-be37-44a8-9150-66980c261a56
let
	img = Gray.(testimage("lighthouse.png"))
	HD = zeros(Complex,size(img))
	C = (size(img,1)÷2, size(img,2)÷2)
	display("Rayon du filtre a disque: $rD")
	for (i,u) in enumerate(range(-0.5,0.5,size(img,1)))
		for (j,v) ∈ enumerate(range(-0.5,0.5,size(img,2)))
			HD[i,j] = u^2 + v^2 <= rD^2 ? true : false
		end
	end
	IMG = fft(Float64.(img))
	img_blur = (fftshift(IMG) .* HD) |> ifftshift |> ifft .|> real .|> Gray
	[img Gray.(abs.(HD)) img_blur]
end

# ╔═╡ f1b4c9df-2685-457d-a592-28a9c6bae95d
md"""
#### Filtre de Butterworth (Lissage)
```math
H_{bw}[u,v]=\frac{1}{ (1 + (D/D_0)^{2n})^{1\over2}} \text{ avec } D^2=u^2+v^2 \text{ et } n \text{ l'ordre du filtre}
```
"""

# ╔═╡ 39befd8e-d2ca-4bf2-adad-02c22e7a4d0e
md"Ordre du filtre n = $(@bind n Slider(collect(1:9),default=2, show_value=true))"

# ╔═╡ f7693c86-20e3-47cd-8cbd-e28a693ccd0b
md"``R_{bw}`` = $(@bind Rbw Slider(0.05:0.05:0.5, default=0.1, show_value=true))"

# ╔═╡ a535babf-82b8-4228-9bae-c264403ce7a0
let
	img = Gray.(testimage("lighthouse.png"))
	Hbw = zeros(Complex,size(img))
	C = (size(img,1)÷2, size(img,2)÷2)
	display("Diametre du filtre BW: $Rbw")
	for (i,u) in enumerate(range(-0.5,0.5,size(img,1)))
		for (j,v) ∈ enumerate(range(-0.5,0.5,size(img,2)))
			squared_r = Float64(u^2 + v^2)
			Hbw[i,j] = 1/sqrt(1+(squared_r/Rbw^2)^n)
		end
	end
	IMG = fft(Float64.(img))
	img_blur = (fftshift(IMG) .* Hbw) |> ifftshift |> ifft .|> real .|> Gray
	[img Gray.(abs.(Hbw)) img_blur]
end

# ╔═╡ 8b67f7f0-b17a-4172-a67f-b8e575d79a23
md"""
#### Filtre de Butterworth passe-haut ( Detecteur de détails)
```math
H_{bw}[u,v]=\frac{1}{ (1 + (D_0/D)^{2n})^{1\over2}} \text{ avec } D^2=u^2+v^2 \text{ et } n \text{ l'ordre du filtre}
```
"""

# ╔═╡ 4beae30a-1773-4a70-aa71-b54df9c2eafa
md"Ordre du filtre n = $(@bind m Slider(collect(1:9),default=2, show_value=true))"

# ╔═╡ ea2b2e2a-a0ce-4a47-8f3c-a04d7398f344
md"``R_{bw}`` = $(@bind R1bw Slider(0.05:0.05:0.5, default=0.1, show_value=true))"

# ╔═╡ da7305dc-5285-4004-842b-c33ea655d6d6
let
	img = Gray.(testimage("lighthouse.png"))
	Hbw = zeros(Complex,size(img))
	#display("Diametre du filtre BW: $R1bw")
	for (i,u) in enumerate(range(-0.5,0.5,size(img,1)))
		for (j,v) ∈ enumerate(range(-0.5,0.5,size(img,2)))
			squared_r = Float64(u^2 + v^2)
			Hbw[i,j] = 1/sqrt(1+(R1bw^2/squared_r)^m) 
		end
	end
	IMG = fft(Float64.(img))
	img_blur = (fftshift(IMG) .* Hbw) |> ifftshift |> ifft .|> real .|> Gray
	[img Gray.(abs.(Hbw)) 20img_blur]
end

# ╔═╡ a1c17aa7-3b4a-4c92-8ece-3e2135399c0a
md"""
#### Filtre de Butterworth (Passe-bande)
```math
H_{bw}[u,v]=\frac{1}{ (1 + (D_0/D)^{2n})^{1\over2}} \times \frac{1}{(1 + (D/D_1)^{2n})^{1\over2}}  \text{ avec } D_1 > D_0
```
"""

# ╔═╡ db8be1fa-b401-46b5-94cb-22aa3f221435
let
	img = Gray.(testimage("lighthouse.png"))
	Hbw = zeros(Complex,size(img))
	display("Diametre du filtre BW: $Rbw")
	for (i,u) in enumerate(range(-0.5,0.5,size(img,1)))
		for (j,v) ∈ enumerate(range(-0.5,0.5,size(img,2)))
			squared_r = Float64(u^2 + v^2)
			Hbw[i,j] = 1/sqrt(1+(Rbw^2/squared_r)^n) * 1/sqrt(1+(squared_r/(1.15Rbw)^2)^n)
		end
	end
	IMG = fft(Float64.(img))
	img_blur = (fftshift(IMG) .* Hbw) |> ifftshift |> ifft .|> real .|> Gray
	[img Gray.(abs.(Hbw)) 10img_blur]
end

# ╔═╡ 64b4355a-bc5d-4cfb-ba18-9cc4e4739147
md"""
#### Filtre Gaussian Spectrale
```math
G[u,v] = \frac{1}{2\pi\sigma^2} \exp\left(-\frac{r^2}{2\sigma^2}\right) \text{ avec } r^2=u^2+v^2
```
"""

# ╔═╡ 4b3d8789-6305-4c07-a71c-1a00cff0ccbb
md"Variance du filtre gaussien = $(@bind σ1 Slider(0.1:0.1:0.5, default=0.1,show_value=true))"

# ╔═╡ 45bb2640-f25e-4fd9-9112-ef5c5ef17d64
let
	img = Gray.(testimage("lighthouse.png"))
	G = zeros(Complex,size(img))
	C = (size(img,1)÷2, size(img,2)÷2)
	display("variance du filtre Gaussien: $(round(σ1^2,digits=3))")
	for (i,u) in enumerate(range(-0.5,0.5,size(img,1)))
		for (j,v) ∈ enumerate(range(-0.5,0.5,size(img,2)))
			square_r = Float64(u^2 + v^2)
			G[i,j] = exp(-square_r/(2σ1^2))
		end
	end
	#G/=sum(real.(G))
	IMG = fft(Float64.(img))
	img_blur = (fftshift(IMG) .* G) |> ifftshift |> ifft .|> real .|> Gray
	[img Gray.(abs.(G)) img_blur]
	#sum(real.(G))
end

# ╔═╡ e3e40fcc-43d6-42c1-bb30-3cc42b3835c0
md"### Filtrage non-lineaire"

# ╔═╡ bb388c66-9de1-4498-8818-79914a47b12c
md"#### Exemple 1: denoising an image with linear filter"

# ╔═╡ bdef3a5c-982a-467c-ade4-6dcb2ead59b5
md"""
S&P density = $(@bind density Slider([5,10,15,20],default=5, show_value=true))  ``\quad`` Filter Type = $(@bind ftype Select(["Box" => "📦", "Gaussian" => "🔔"]))
"""

# ╔═╡ 7a75ffb8-a863-451e-956c-34342391195f
let 
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/livingroom.tif") |> TiffImages.load
	img_sp = load("cv_pics/livingroom_sp_d$(density).png")
	kern_dim = 5
	if ftype == "Box"
		filter = centered(ones(kern_dim, kern_dim)/kern_dim^2)
	else
		filter = Kernel.gaussian(1)
	end
	img_blur = imfilter(img_sp, filter)
	display("PSNR : $(round(assess_psnr(img_blur, img),digits=2))")
	[img img_sp img_blur]
end

# ╔═╡ 588c45c5-d649-4ac7-b37a-74f0e04c7a71
md""" 
#### Filtre non lineaire

##### Filtre d'ordre : Mediane

"""

# ╔═╡ bc64694c-0914-43f0-9df0-78b946b86945
md"""
Median window size = $(@bind K Slider([1,2,3,4,5], default=2, show_value=true))

+ Shot noise
"""

# ╔═╡ a936c86b-d683-42ca-9ba2-45eb92b2d75c
let 
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/livingroom.tif") |> TiffImages.load
	img_sp = load("cv_pics/livingroom_sp_d$(density).png")
	img_sp_ext = [repeat([img_sp[1,1]],K,K) repeat(img_sp[1,:],1,K)' repeat([img_sp[1,end]],K,K)
	repeat(img_sp[:,1],1,K) img_sp repeat(img_sp[:,end],1,K)
	repeat([img_sp[end,1]],K,K) repeat(img_sp[end,:],1,K)' repeat([img_sp[end,end]],K,K)]
	img_filt = similar(img)
	for i=K+1:size(img_sp_ext,1)- K
		for j=K+1:size(img_sp_ext,2)-K
			V=img_sp_ext[i-K:i+K,j-K:j+K]
			img_filt[i-K,j-K] = median(sort(V[:]))
		end
	end
	display("PSNR : $(round(assess_psnr(img_filt, img),digits=2))")
	[img img_sp img_filt]
end

# ╔═╡ 3ce0535c-e231-4a06-8249-e2fe89249f60
md"+ Gaussian noise"

# ╔═╡ 9453f9aa-8e0e-4855-8b10-4d4176d297aa
let 
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/livingroom.tif") |> TiffImages.load
	img_sp = load("cv_pics/livingroom_gn.png")
	img_sp_ext = [repeat([img_sp[1,1]],K,K) repeat(img_sp[1,:],1,K)' repeat([img_sp[1,end]],K,K)
	repeat(img_sp[:,1],1,K) img_sp repeat(img_sp[:,end],1,K)
	repeat([img_sp[end,1]],K,K) repeat(img_sp[end,:],1,K)' repeat([img_sp[end,end]],K,K)]
	img_filt = similar(img)
	for i=K+1:size(img_sp_ext,1)- K
		for j=K+1:size(img_sp_ext,2)-K
			V=img_sp_ext[i-K:i+K,j-K:j+K]
			img_filt[i-K,j-K] = median(sort(V[:]))
		end
	end
	display("PSNR : $(round(assess_psnr(img_filt, img),digits=2))")
	[img img_sp img_filt]
end

# ╔═╡ 89b15473-61e3-41ec-b59b-93d194c5e423
md"+ α-trimming mean algorithm"

# ╔═╡ 37feb1cc-bdf7-4097-805e-e6eca7b8801a
md"α = $(@bind α Slider(0.1:0.1:0.5, default=0.2, show_value=true))"

# ╔═╡ 60791d88-a51e-4206-ab9f-848c4e93f322
let 
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/livingroom.tif") |> TiffImages.load
	img_sp = load("cv_pics/livingroom_sp_d$(density).png")
	img_sp_ext = [repeat([img_sp[1,1]],K,K) repeat(img_sp[1,:],1,K)' repeat([img_sp[1,end]],K,K)
	repeat(img_sp[:,1],1,K) img_sp repeat(img_sp[:,end],1,K)
	repeat([img_sp[end,1]],K,K) repeat(img_sp[end,:],1,K)' repeat([img_sp[end,end]],K,K)]
	img_filt = similar(img)
	T = ceil(Int,α*(2K+1))
	for i=K+1:size(img_sp_ext,1)- K
		for j=K+1:size(img_sp_ext,2)-K
			V=img_sp_ext[i-K:i+K,j-K:j+K]
			img_filt[i-K,j-K] = mean(sort(V[:])[T:end-T])
		end
	end
	display("PSNR : $(round(assess_psnr(img_filt, img),digits=2))")
	[img img_sp img_filt]
end

# ╔═╡ 308b820b-f354-4822-af35-c5f12ca3a62a
md"+ Adaptative Median Filter"

# ╔═╡ 1ee4fb66-b23a-4aa5-b5b0-f2f734d3b242
let 
	img = download("https://raw.githubusercontent.com/JuliaImages/TestImages.jl/images/images/livingroom.tif") |> TiffImages.load
	img_sp = load("cv_pics/livingroom_sp_d$(density).png")
	img_filt = copy(img_sp)
	for K=1:5
		img_sp_ext = [repeat([img_sp[1,1]],K,K) repeat(img_sp[1,:],1,K)' repeat([img_sp[1,end]],K,K); repeat(img_sp[:,1],1,K) img_sp repeat(img_sp[:,end],1,K); repeat([img_sp[end,1]],K,K) repeat(img_sp[end,:],1,K)' repeat([img_sp[end,end]],K,K)]
		for i=K+1:size(img_sp_ext,1)- K
			for j=K+1:size(img_sp_ext,2)-K
				if img_filt[i-K,j-K] != img_sp[i-K,j-K]
					continue
				else
					V=img_sp_ext[i-K:i+K,j-K:j+K]
					Med= median(sort!(V[:]))
					σ1 = Float64.(Med) - Float64.(V[1])
					σ2 = Float64(Med) - Float64(V[end])
					if σ1 < 0 || σ2 < 0
						img_filt[i-K,j-K]=Med
					end
				end
			end
		end
	end
	display("PSNR : $(round(assess_psnr(img_filt, img),digits=2))")
	[img img_sp img_filt]
end

# ╔═╡ 13d0dbd0-9666-488a-bf5e-3282692b4032
# ╠═╡ disabled = true
# ╠═╡ skip_as_script = true
#=╠═╡
begin
	im = load("cv_pics/cap19.png") .|> Gray
	imb = @. Gray(im<0.5)
end;
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─962e4274-38fd-11ed-3890-bfb3c8d0f09d
# ╟─2fb557e2-258f-4d47-81a3-495ef846b080
# ╟─b4a6ae81-bd2c-4245-b22b-47c5c317d417
# ╟─1cbf1b09-c594-4415-8301-7641944cc599
# ╟─63be3bce-7a1d-41b7-918d-00cccbc575df
# ╟─3308305f-a842-4dd4-b03a-500bd7fa4f8b
# ╟─9edd311a-b2a2-46f8-89d2-e21a71edcaef
# ╟─f1b20ca8-2e71-49a4-98a0-da4d7ac9ed1a
# ╟─9eb8287e-b81b-48c7-8b97-72cdcdc86f20
# ╟─98203216-6988-4e62-b1d5-eaf355cb6c1a
# ╟─bf3c567d-3cc5-472c-a4c0-bb98bebfa1d6
# ╟─caf20bd3-8a7a-4025-b3df-76f8c9f61458
# ╟─ccb00a2d-6c86-4490-af16-b47dd66ff757
# ╟─ed476ab1-ab27-4907-a305-78a9eed21213
# ╟─fbad8364-e23c-4303-90b8-d71181cecdc7
# ╟─104cdab1-fd83-4ae4-bd53-2d685156c6c4
# ╟─580d6634-9b7b-4387-9081-ad7e9d6b33be
# ╟─0485cdd7-89d2-45b0-a0da-e62fb1a78cf6
# ╟─c437baad-cc24-418a-b176-e295e8a812a6
# ╟─2200a9aa-34e1-4844-9f92-1370906b4422
# ╟─f64d7431-9edb-4872-aa73-181408c9c497
# ╟─d5e6c86f-a341-46a9-9658-e18dda02a52c
# ╟─7b05dc85-982c-4a94-adb8-e64e77464e12
# ╟─958c9316-bf30-4069-8cb0-5f63f4bfa9de
# ╟─bcc1e8c0-824b-469f-8ef8-497151054a02
# ╟─116f8969-a810-4c87-a6d2-69b1a510f985
# ╠═abcea042-2cde-4413-8443-7df5cfaa58ce
# ╟─f20e0585-7af6-4faa-9da9-903ceee3ca83
# ╟─8b08079c-565d-4cd2-b48b-42ab4d50dc46
# ╟─0e2fef44-65e4-45bd-a126-1c15c598121c
# ╟─7b5dcbec-68b9-4829-bb8d-2685c0e8ae08
# ╟─da92ef0c-3123-4a30-bea7-e5a65cf412e7
# ╟─b2ef51f4-81c9-4a61-bd04-d35213ae42ed
# ╟─0323d736-82a4-4220-bb80-2317cddb65f2
# ╟─7a293c05-7bee-44db-8e49-3457c92b38fd
# ╟─11114219-be37-44a8-9150-66980c261a56
# ╟─f1b4c9df-2685-457d-a592-28a9c6bae95d
# ╟─39befd8e-d2ca-4bf2-adad-02c22e7a4d0e
# ╟─f7693c86-20e3-47cd-8cbd-e28a693ccd0b
# ╟─a535babf-82b8-4228-9bae-c264403ce7a0
# ╟─8b67f7f0-b17a-4172-a67f-b8e575d79a23
# ╟─4beae30a-1773-4a70-aa71-b54df9c2eafa
# ╟─ea2b2e2a-a0ce-4a47-8f3c-a04d7398f344
# ╟─da7305dc-5285-4004-842b-c33ea655d6d6
# ╟─a1c17aa7-3b4a-4c92-8ece-3e2135399c0a
# ╟─db8be1fa-b401-46b5-94cb-22aa3f221435
# ╟─64b4355a-bc5d-4cfb-ba18-9cc4e4739147
# ╟─4b3d8789-6305-4c07-a71c-1a00cff0ccbb
# ╟─45bb2640-f25e-4fd9-9112-ef5c5ef17d64
# ╟─e3e40fcc-43d6-42c1-bb30-3cc42b3835c0
# ╟─bb388c66-9de1-4498-8818-79914a47b12c
# ╟─bdef3a5c-982a-467c-ade4-6dcb2ead59b5
# ╟─7a75ffb8-a863-451e-956c-34342391195f
# ╟─588c45c5-d649-4ac7-b37a-74f0e04c7a71
# ╟─bc64694c-0914-43f0-9df0-78b946b86945
# ╟─a936c86b-d683-42ca-9ba2-45eb92b2d75c
# ╟─3ce0535c-e231-4a06-8249-e2fe89249f60
# ╟─9453f9aa-8e0e-4855-8b10-4d4176d297aa
# ╟─89b15473-61e3-41ec-b59b-93d194c5e423
# ╟─37feb1cc-bdf7-4097-805e-e6eca7b8801a
# ╟─60791d88-a51e-4206-ab9f-848c4e93f322
# ╟─308b820b-f354-4822-af35-c5f12ca3a62a
# ╟─1ee4fb66-b23a-4aa5-b5b0-f2f734d3b242
# ╟─13d0dbd0-9666-488a-bf5e-3282692b4032
# ╟─de935179-fd70-4a68-8f7f-e23b2b9187aa
