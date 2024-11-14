### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ e2f052ba-7de7-4ff9-a7e2-95d8f8a8d25d
begin
	import Pkg
	Pkg.activate(pwd())
end

# ╔═╡ d191e543-f478-4486-8452-39e1b9328745
begin
	using TikzPictures
	using Images
end

# ╔═╡ 031e796e-2d49-11ed-3566-1fc673c08edd
md"""
# Vision par ordinateur
## Filiere GSEII - S5
### Preparee par : Pr. H.Belkebir
#### Annee universitaire: 2022/2023
"""

# ╔═╡ bd1b4c8e-2cc1-4eef-a487-00a80b5b502a
md"""
## Table de matiere

### I. Introduction generale
### II. Formation et representation des images.
### III. Pre-traitements des images numeriques.
### V. Transformations geometrique
### VI. Segmentation des images.
"""

# ╔═╡ 8ad9b580-5e06-4653-a58b-dd965a79892c
md"""
## Remarques d'ordre generale:

__Le long de ce cours nous allons presenter les notions globales et les techniques de base concernant chaque partie. Nous allons aussi accompagner ces notions par des exemples qui illustrent comment mettre en pratique la partie theorique du cours.
Pour ce faire nous allons faire appel a differentes bibliotheque de traitement d'image telle que OpenCV, PILLOW, Images et skImage.__
"""

# ╔═╡ 6959b3ba-751a-45e9-a686-2a071b553284
md"### I. Introduction generale:"

# ╔═╡ d16bb1a5-8e66-48fd-a12a-72072de645f2
md"""
#### 1. Qu'est ce que la vision par ordinateur ?
"""

# ╔═╡ 76cf9520-6993-45d2-bb2c-57d8ce253853
TikzPicture(L"""
\tikzset{every node/.style={rectangle, rounded corners, inner sep=12pt, font=\Large\bfseries, text=orange!75!black, align=center}}
\node[draw](0,0)(V){Vision};
\node[draw, below=2cm of V, text width=4cm](EV){Perception du monde a travers la lumiere};
\node[style={},right = 2.5cm of V, above]{par};
\node[draw, right = 5cm of V](C){Ordinateur};
\node[draw, below=2cm of C, text width=8cm](EC){Calculateur electronique utilisant des algorithmes numeriques pour executer des taches complexes};
\draw (V) to (C);
\draw[dashed, gray] (V) to (EV);
\draw[dashed, gray] (C) to (EC); 
""", options="scale=1.0", preamble="\\usetikzlibrary{positioning}")

# ╔═╡ 842d5c7e-49d1-4002-9b4d-14f7c11ce49f
md"""
#### 2. Que cherche a accomplir la vision artificielle ?
> **In computer vision, we are trying to do the inverse, i.e., to describe the world that we see in one or more images and to reconstruct its properties, such as shape, illumination, and color distributions.** (*Richard Szeliski - Computer Vision: Algorithms and Applications*)
"""

# ╔═╡ e3e19cd0-3132-44af-b483-8423964e2f14
md"""
##### 1. Exemples de vision artificielle
- OCR (optical caracter recognition).

- Machine Inspection.

- Retail: object recognition for automated checkout lanes and fully automated stores;

- Warehouse logistics: autonomous package delivery and pallet-carrying “drives” and parts picking by robotic manipulators;

- Medical imaging: registering pre-operative and intra-operative imagery or performing long-term studies of people’s brain morphology as they age;

- Self-driving vehicles: capable of driving point-to-point between cities as well as autonomous flight;

- 3D model building (photogrammetry): fully automated construction of 3D models from aerial and drone photographs;

- Match move: merging computer-generated imagery (CGI) with live action footage by tracking feature points in the source video to estimate the 3D camera motion and shape of the environment. Such techniques are widely used in Hollywood, e.g., in movies such as Jurassic Park ; they also require the use of precise matting to insert new elements between foreground and background elements.

- Motion capture (mocap): using retro-reflective markers viewed from multiple cameras or other vision-based techniques to capture actors for computer animation;

- Surveillance: monitoring for intruders, analyzing highway traffic and monitoring pools for drowning victims;

- Fingerprint recognition and biometrics: for automatic access authentication as well as forensic applications.
"""

# ╔═╡ 6ecf0bc4-7184-4c5d-bdbf-92a7d6aaf826
md"""
##### 2. Timeline of computer vision
$(load("caps/computer_vision_timeline.png"))
"""

# ╔═╡ ca75d1e5-85e6-40f1-90d6-26254449c469
html"""
<h4> 3. Relation avec le traitement d'image.</h4>
<h5>1. Qu'est ce que le traitement d'image:</h5>
<br>
<span style="color:blue; font-wight:bold">C'est l'ensemble des techniques et methodes utilisees pour faciliter et concretiser le processus d'exploitation de l'information transportee par l'image. Elles sont organizee en trois familles:
<ul>
<li> Traitement de bas niveau (debruitage, amelioration du contraste et des detailles de l'image ...)
<li> Traitement de moyen niveau (Segmentation et classification des objects)
<li> Traitement de haut niveau ( Description semantique de l'image)
</ul>
</span>
"""

# ╔═╡ 986ca0d8-5736-48d4-9c5e-be3b017b023e
TikzPicture(L"""
\tikzset{noeuds/.style={draw,rectangle,thick,blue!75!black,text=black!75!white,rounded corners,font=\sf\scriptsize,inner sep=6pt,align=left,text width=0.125\textwidth}}
\tikzset{implies/.style={blue!75!gray,implies-implies, double distance=3pt, thick}}
\node[noeuds,font=\Large\bfseries,text width=0.3\linewidth,inner ysep=60pt, inner xsep=12pt,align=center] (BK)at (0,0) {Base de connaissance};
\node[noeuds](aq) at (-128pt,-70pt) {Image acquisitions};
\node[noeuds](IE) at (-128pt,-23pt) {Image Enhancing};
\node[noeuds](IR) at (-128pt,24pt) {Image restauration};
\node[noeuds](CIP) at (-128pt,70pt) {Color image processing};
\node[noeuds](MRP) at (-44pt,116pt) {Multi-resolution processing};
\node[noeuds,inner ysep=10] (ICOM) at (40pt,116pt) {Image compression};
\node[noeuds](IOR) at (128pt,-70pt) {Object Recognition};
\node[noeuds](IRD) at (128pt,-23pt) {Representation and description};
\node[noeuds](IS) at (128pt,24pt) {Segmentation};
\node[noeuds](IMP) at (128pt,70pt) {Morphological processing};
\draw[implies] (aq.east) to  (-65pt,-70pt);
\draw[implies] (IE.east) to  (-65pt,-23pt);
\draw[implies ](IR.east) to  (-65pt,24pt);
\draw[implies] (CIP.east) to  (-65pt,70pt);
\draw[implies] (MRP.south) to  (-44pt,78pt);
\draw[implies ](ICOM.south) to  (40pt,78pt);
\draw[implies] (IOR.west) to  (65pt,-70pt);
\draw[implies] (IRD.west) to  (65pt,-23pt);
\draw[implies ](IS.west) to  (65pt,24pt);
\draw[implies ](IMP.west) to  (65pt,70pt);
\node[font=\Large\bfseries,orange!75!black] at (0,-95pt){Image Processing Synopsis};
""",options="scale=1.05", preamble="\\usetikzlibrary{arrows,decorations.pathmorphing,backgrounds,positioning,fit,matrix}")

# ╔═╡ 816bba41-8fe7-47d2-96be-085c9be0e0e9
md"""
##### 2. Lien avec la vision artificielle.
"""

# ╔═╡ 2757f321-5c90-4cf5-ba10-aa5a498c6967
TikzPicture(L"""
\tikzset{every node/.style={rectangle, rounded corners, inner sep=12pt, font=\Large\bfseries, text=orange!75!black, align=center}}
\node[draw](0,0)(ImP){Traitement d'image};
%\node[style={},right = 2.5cm of ImP, above]{};
\node[draw, above = 2cm of ImP](C){Vision par Ordinateur};
\draw[>=stealth,double, ->, gray, thick](ImP) to (C);
""", options="scale=1.0", preamble="\\usetikzlibrary{positioning, arrows, arrows.meta}")

# ╔═╡ Cell order:
# ╟─e2f052ba-7de7-4ff9-a7e2-95d8f8a8d25d
# ╟─d191e543-f478-4486-8452-39e1b9328745
# ╟─031e796e-2d49-11ed-3566-1fc673c08edd
# ╟─bd1b4c8e-2cc1-4eef-a487-00a80b5b502a
# ╟─8ad9b580-5e06-4653-a58b-dd965a79892c
# ╟─6959b3ba-751a-45e9-a686-2a071b553284
# ╟─d16bb1a5-8e66-48fd-a12a-72072de645f2
# ╟─76cf9520-6993-45d2-bb2c-57d8ce253853
# ╟─842d5c7e-49d1-4002-9b4d-14f7c11ce49f
# ╟─e3e19cd0-3132-44af-b483-8423964e2f14
# ╟─6ecf0bc4-7184-4c5d-bdbf-92a7d6aaf826
# ╟─ca75d1e5-85e6-40f1-90d6-26254449c469
# ╟─986ca0d8-5736-48d4-9c5e-be3b017b023e
# ╟─816bba41-8fe7-47d2-96be-085c9be0e0e9
# ╟─2757f321-5c90-4cf5-ba10-aa5a498c6967
