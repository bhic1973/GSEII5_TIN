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

# ╔═╡ bf73d046-36a6-46f3-9b62-b8815920671b
begin
	import Pkg
	Pkg.activate(".")
end

# ╔═╡ e8cd9497-fa73-46bc-989e-8c1816b2a7c7
begin
	using PlutoUI , LaTeXStrings
	using Plots, Images
	using DSP
	using TikzPictures
	using ReferenceFrameRotations
end

# ╔═╡ 6da4b103-e1c5-465b-aa92-47aa2cab43a2
using PyCall

# ╔═╡ d048d350-89ab-11eb-2ff5-038da6d710e3
# ╠═╡ disabled = true
#=╠═╡
using FFTW
  ╠═╡ =#

# ╔═╡ adc8674a-89a2-11eb-2caf-1915cd540d9e
md"# Formation et representation numerique d'une image"

# ╔═╡ b6409d51-9730-4fe3-b829-a13894103ba1
#im_file = download("https://encrypted-tbn0.gstatic.com/images?q=tbn:ANd9GcTXyb_LDJ4wTVusoWyI53Nd-HeIfrDacsqHGw&usqp=CAU");
imFile = "tin_pics/cap01.png";

# ╔═╡ fd1863f7-dd7c-4d9b-af44-8c9488be79cd
md" ## A. Formation des images"

# ╔═╡ 915aee51-5be0-4ba2-81ea-15639f6b30c5
TikzPicture(L"""
\tikzset{
    every node/.style={
		draw,
		text=black!75!gray, 
		font=\large\bfseries, 
		inner sep=9pt, 
		rectangle, 
		%rounded corners,
		fill = orange!15,
		solid
	},
	every path/.style={
		dashed,
	}
}
\draw[thick, blue!75!black](0cm,0cm)node[]{Monde} to++(20cm,0cm)node[above=0.25cm, inner ysep=6pt, anchor=south]{Capteur} to (40cm,0cm)node[]{Image};
"""
, options="scale=0.25", preamble="\\usetikzlibrary{positioning}")

# ╔═╡ 14100cb9-52fa-425c-97a1-0477c4c70654
md"""
$(load(joinpath(@__DIR__(),"cv_pics/cap01.png")))

```math
	im(x,y) = i(x,y) \times r(x,y)
```
avec:

``i(x,y)``: Fonction d'illumination ``(0 \rightarrow \infty)``

``r(x,y)``: Fonction de reflectance ``(0 \rightarrow 1)``
"""

# ╔═╡ 5a116364-af9e-4571-ae6e-2df7450c92bb
md"""
Cannaux de couleur de l'image: $(@bind c Select(["1"=>"Red", "2"=>"Green", "3"=>"Blue"]))
"""

# ╔═╡ 504b6ae5-fa6e-4be4-85da-e36e614514d8
md"""
### I. Composantes principales du processus de formation des images 
$(load("cv_pics/cap02.png"))
"""

# ╔═╡ 736fbc4f-4559-4131-8674-ecdd2ee8ce8f
md""" #### 1. Transformation geometrique 
 
Une image est le fruit de l'interaction de quatre elements :
- Les condition de l'illumination scenique.
- La geometrie de la scene.
- Les proprietes surfacique des composantes de la scene
- L'optique et l'electronique du capteur

"""

# ╔═╡ 33f6b6be-cc82-476a-85da-2a6f98505a6a
md""" 
##### 1.1 Primitifs geometriques

Pour pouvoir decrire les composantes d'une scene mathematiquement, on a besoin d'utiliser les primitifs geometrique 2D et 3D a savoir:
- Les points;
- Les droites;
- Les plans;
- Les lignes courbes;
- Les surfaces courbes;
- ``\dots``

"""

# ╔═╡ 1aee84e1-eadc-48e2-a702-0b31ee16e0b2
md"""
+ __Points 2D__

Un point est decrit dans l'espace 2D par une pair de valeurs reelles:
```math
P \in \mathbb{R}^2\quad \text{tel que } \mathbf{P}=\begin{bmatrix} x \\ y\\ \end{bmatrix}
```
Parfois on utilise aussi la representation en coordonnees homogenes:
```math
\tilde{P} \in \mathcal{P}^2 \quad \text{tel que } \mathbf{\tilde{P}} = [\tilde{x}, \tilde{y}, \tilde{w}]^T
```
La relation entre les deux representations est donnee par:

```math
\mathbf{\tilde{P}} =\tilde{w}\times[x, y, 1]^T
``` 
Avec ``\bar{x}=[x, y, 1]^T`` est le vecteur augmentee.

"""

# ╔═╡ 60358d7e-8cbc-4d53-9167-2ab36e660ceb
md"""
+ __Droite 2D__

Une droite est entierement definis par:

- Deux points ``\mathbf{\bar{x}}`` et ``\mathbf{\bar{y}}``;
- Un point ``\mathbf{\bar{x}}`` et la normale ``\mathbf{\tilde{n}}`` a ``D`` en ce point.

soit ``\mathbf{\tilde{u}}=[a, b, c]^T`` un vecteur quelconque et soit ``\mathbf{\bar{x}}`` un point du plan. On peut definir une droite ``D`` par l'annulation du produit scalaire suivant:
```math
\mathbf{\tilde{u}}\cdot\mathbf{\bar{x}} = a\,x\; +\; b\,y\;+\;c=0
```
La normale de ``D`` est donnee par:
```math
\mathbf{\tilde{n}} = [\tilde{n}_x, \tilde{n}_y, d]^T = \mathbf{\tilde{u}}/\sqrt{a^2+b^2}
```
$(load("cv_pics/cap03.png"))

Il est aussi possible de representer la normale ``\mathbf{\tilde{n}}``  par l'angle de rotation ``\theta`` par rapport a ``(Ox)`` 
```math
\Longrightarrow \mathbf{\tilde{n}} = [\cos(\theta), \sin(\theta), d]^T.
```
>_NB: En representation homogene, l'intersection __x__ de deux droites est donnee par le produit vectorielle suivant:_
```math
\mathbf{\tilde{x}} = \mathbf{\tilde{u}}_1 \wedge \mathbf{\tilde{u}}_2
```
_Question: La ligne definie par_ ``\mathbf{\tilde{x}}_1 \wedge \mathbf{\tilde{x}}_2``_, qu'est ce qu'elle represente ?_

"""

# ╔═╡ 2aff139c-324e-402c-826c-31f17af7a671
md""" 
+ __Pour s'exercer :__
  
  Least squares intersection point and line fitting—advanced. 
  
The last Equation in the previous section shows how the intersection of two 2D lines can be expressed as their cross product, assuming the lines are expressed as homogeneous coordinates.

1. If you are given more than two lines and want to find a point  ``\tilde{x}`` that minimizes the sum of squared distances to each line:

```math
D = \sum_i( \tilde{x} · \tilde{I}_i)^2
```
how can you compute this quantity?
  
2. To fit a line to a bunch of points, you can compute the centroid (mean) of the points and the covariance matrix of the points around this mean. Show that the line passing through the centroid along the major axis of the covariance ellipsoid (largest eigenvector) minimizes the sum of squared distances to the points.

3. These two approaches are fundamentally different, even though projective duality tells us that points and lines are interchangeable. Why are these two algorithms so apparently different? Are they actually minimizing different objectives?
"""


# ╔═╡ 9d84ac78-f151-49f8-9eea-6b9c09acc10f
md"""
+ __Point 3D__

Un point ``\mathbf{\tilde{x}}`` en trois dimensions est defini dans la representation homogene par le quadruplet: 
```math
\tilde{w}[x, y, z, 1]^T \in \mathcal{P}^3
```
+ __Plan 3D__

Soit ``\mathbf{\bar{x}} = [x,y,z,1]^T`` un point 3D et soit ``\mathbf{\tilde{u}}=[a,b,c,d]^T`` un vecteur 3D. On construit le plan ``P``  en utilisant l'equation suivante:
```math
\mathbf{\tilde{x}} \cdot \mathbf{\tilde{u}} =  a\,x\;+\;b\,y\;+\;c\,z\;+\;d=0
```
Nous pouvons aussi normaliser l'equation de ``P`` et extraire la normale ``\mathbf{\tilde{n}} = [\tilde{n}_x, \tilde{n}_y, \tilde{n}_z, e]`` en utilisant le meme approche elaboree pour le cas 2D. De meme pour la formulation spherique en utilisant les angles ``\phi`` et ``\theta``: 
```math
\mathbf{\tilde{n}} = [\cos(\phi)\cos(\theta), \sin(\phi)\cos(\theta), \sin(\theta), e]
```
+ __Lignes 3D__

On definit une ligne droite ``D`` dans un espace 3D en utilisant deux points ``\mathbf{\tilde{x}}`` et ``\mathbf{\tilde{y}}`` appartenant a ``D`` selon l'equation ci-apres:
```math
\mathbf{\tilde{M}} = \lambda\, \mathbf{\tilde{x}}\, +\, (1\,-\,\lambda)\,\mathbf{\tilde{y}}
```
"""

# ╔═╡ 719bdf63-972b-4c3a-88f7-4c6fa0f7e054
md"""
##### 1.2 Transformation 2D

+ __Translation__
  
  Une trnslation est definie par le vecteur de translation ``\mathbf{t}`` et le point de depart ``\mathbf{x}_0`` tel que ``\mathbf{x} = \mathbf{x}_0 +\mathbf{t}``
  
  La representation homogene permet de formuler cette operation en produit matrice-vecteur:
```math
\mathbf{\tilde{x}} = \mathcal{T} \mathbf{\tilde{x}}_0 = \begin{bmatrix} 
1 & 0 & t_x  \\
0 & 1 & t_y \\
0 & 0 & 1 \\
\end{bmatrix} \mathbf{\tilde{x}}_0
```
+ __Transformation euclidienne: rotation et translation__
La rotation d'un angle ``\theta`` est definie par l'expression suivante:

```math
\mathbf{x} = \mathbf{R} \mathbf{x}_0 + \mathbf{t}\quad \text{avec }\mathbf{R}=\begin{bmatrix} \cos(\theta) & -\sin(\theta)\\\sin(\theta) & \cos(\theta)\\\end{bmatrix}
```
Ce qui permet d'obtenir en representation homogene
```math
\mathbf{\tilde{x}} = \mathcal{T} \mathbf{\tilde{x}}_0 = \begin{bmatrix} 
\cos(\theta) & -\sin(\theta) & t_x  \\
\sin(\theta) & \cos(\theta) & t_y \\
0 & 0 & 1 \\
\end{bmatrix} \mathbf{\tilde{x}}_0
```
+ __Transformation de Similaritee : rotation + Changement d'echelle__
Le changement d'echelle d'un facteur ``s`` est realisee en appliquant la relation suivante:
```math
\mathbf{x} = s\mathbf{R} \mathbf{x}_0 + \mathbf{t}
```

+ __Transformation affine:__

  Soit ``A`` une matrice ``2\times3`` a coefficients reels, la transformation affine est definie en representation homogene comme suit:
```math
\mathbf{\tilde{x}} = \mathcal{A} \mathbf{\tilde{x}}_0 = \begin{bmatrix} 
a_1 & a_2 & a_3  \\
a_4 & a_5 & a_6 \\
0 & 0 & 1 \\
\end{bmatrix} \mathbf{\tilde{x}}_0
```

+ __Transformation projective (Homographie):__

  Une transformation projective est realise par une matrice ``3\times3`` à coefficient réels:
```math
\mathbf{\tilde{x}} = \mathcal{H} \mathbf{\tilde{x}}_0 = \begin{bmatrix} 
a_1 & a_2 & a_3  \\
a_4 & a_5 & a_6 \\
a_7 & a_8 & 1 \\
\end{bmatrix} \mathbf{\tilde{x}}_0
```
$(load("cv_pics/cap04.png"))
"""

# ╔═╡ 4146099c-7530-461f-967d-ca2b60ebda24
md"""
+ __Exercez vous :__

Write a program that lets you interactively create a set of rectangles and then modify their "pose" (2D transform). You should implement the following steps:

1. Open an empty window (“canvas”).

2. Shift drag (rubber-band) to create a new rectangle.

3. Select the deformation mode (motion model): translation, rigid, similarity, affine, or perspective.
4. Drag any corner of the outline to change its transformation.

This exercise should be built on a set of pixel coordinate and transformation classes, either implemented by yourself or from a software library. Persistence of the created representation (save and load) should also be supported (for each rectangle, save its transformation).

Suggestion: use python turtle or opencv module
"""

# ╔═╡ 1dc7eb8d-a04d-4f24-bd8d-97c6cf67a213
md"""
##### 1.3 Transformation 3D

Les definitions etablies en 2D restent valables pour la 3D avec un traitement adaptee pour la rotation 3D.
 
+ __Formule de Rodriguez:__

  Considerons une rotation autour d'un axe ``\mathbf{n}`` d'un angle ``\theta``. La matrice de rotation est donnee par:

```math
\mathbf{R}(\mathbf{n}, \theta) = \mathbf{I} +\sin(\theta)[\mathbf{n}]_\wedge + (1-\cos(\theta))[\mathbf{n}]^2_\wedge \text{ avec } [\mathbf{n}]_\wedge=\begin{bmatrix}
0 & -n_z & n_y \\
n_z & 0 & -n_x \\
-n_y & n_x & 0 \\
\end{bmatrix}
```
+ __Quaternions unitaires:__

  Un quaternion est un nombre hypercomplexe au sens mathematique de dimension quatre. Il est utilisee pour representer les rotation en 3D et il s'ecrit comme suite :

```math
q = q_0 + q_1\,\imath + q_2\,\jmath + q_3\,\mathit{k}\quad \text{ avec } \imath^2=\jmath^2 = \mathit{k}^2 = \imath\jmath\mathit{k} = -1
```
On parle de quaternion unitaire lorsque ``\mid\mid{q}\mid\mid=1``. Il peut aussi etre derivee de la representation angle/axe de la maniere ci-apres:
```math
q(w,\mathbf{v}) = (\cos(\theta/2), \sin(\theta/2)\times\mathbf{n})
```
Nous pourions alors transformer la formule de Rodriguez pour obtenir la matrice de rotation 3D a partir du quaternion ``q(w,[x,y,z])``:
```math
\mathbf{R}(\mathbf{q})=\begin{pmatrix}
1-2(y^2+z^2) & 2(xy - wz) & 2(xz+wy)\\
2(xy+wz) & 1-2(x^2+z^2) & 2(yz -wx)\\
2(xz - wy) & 2(yz+wx) & 1-2(x^2+y^2)\\
\end{pmatrix}
```
+ __Exemple:__

  Supposant qu'on realise une rotation autour de l'axe ``OY`` d'un angle ``\theta=60^{\circ}``. Le quaternion representant cette rotation se calcule a partir de la formule deja donnee:

q=$(Quaternion(cosd(30),sind(30)*[0,1,0]))

La matrice de rotation equivalente est :

$(display(quat_to_dcm(Quaternion(cosd(30),sind(30)*[0,1,0])));)
"""

# ╔═╡ b6257893-55f2-4434-ab23-04725d84ba39
md"""
##### 1.4 Projection 3D vers 2D
Connaissant maintenant comment représenter et transformer spatialement les primitifs géométriques de base, il devient alors possible de s'attaquer à la question de leur projection dans le plan de l'image.

+ __Projection Orthographique__

  Elle consiste à supprimer la composante de profondeur de la primitive 3D pour retrouver son correspondant 2D dans le plan de l'image.
  
```math
\mathbf{x}_{2D} = [\mathbf{I}_{2\times2}|\mathbf{0}] \mathbf{x}_{3D}
```
  Ou en représentation homogène :


```math
\mathbf{\tilde{x}}_{2D} = \begin{pmatrix}
1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & 0 & 0 & 0\\
0 & 0 & 0 & 1\\ 
\end{pmatrix} \mathbf{\tilde{x}}_{3D}
```

> _NB: Dans cette représentation, on supprime la contribution de la profondeur, mais on conserve l'effet de l'échelle représenté par le paramètre w_

  Dans la pratique, nous devons ajouter un paramètre d'échelle à ce type de transformation pour adapter les dimensions de la scène aux dimensions du plan du capteur de l'image.

```math
\mathbf{x}_{2D} = [s\mathbf{I}_{2\times2}|\mathbf{0}] \mathbf{x}_{3D}
```
$(load("cv_pics/cap05.png"))
"""

# ╔═╡ 964f2417-f434-4390-ad03-f75fe7b1473d
md"""
+ __Projection perspective:__

En vision par ordinateur, on utilise couremment la projection perspective pour representer les elements 3D d'une scene dans un plan:

```math
\mathbf{\bar{x}} = \mathcal{P}\left\{\begin{bmatrix}x\\y\\z\end{bmatrix}\right\} = \begin{bmatrix}x/z\\y/z\\1\end{bmatrix}
```

ou encore en representation homogene:

```math
\mathbf{\tilde{x}} = \mathcal{P}\{\mathbf{\bar{p}}\}=\begin{pmatrix}
1 & 0 & 0 & 0\\
0 & 1 & 0 & 0\\
0 & 0 & 1 & 0\\
\end{pmatrix}\mathbf{\bar{p}}
```
$(load("cv_pics/cap06.png"))
"""

# ╔═╡ 9bc6ead4-9b52-4e9d-8421-11ddf8546477
md"""
##### 1.5 PinHole Camera Model
$(load("cv_pics/cap07.png"))

C'est un modèle simplifié de la caméra utilisée pour la capture de l'image. Il sert à créer un mappage entre les coordonnées de l'image et les coordonnées 3D de la scène. Il est complètement défini si on connait les paramètres intrinsèques (centre de la caméra, focale, ...) et extrinsèques (position et orientation / référentiel de la scène) de la caméra.

+ __Paramètres intrinsèques de la caméra :__

  Ils sont regroupés dans la matrice de calibrage ``\mathbf{K}`` donnée par :

```math
\mathbf{K} =\begin{pmatrix}f_x & s & c_x\\0 & f_y & c_y\\0 & 0 & 1\\\end{pmatrix}
```

+ __Paramètres extrinsèques :__

Ils sont représentés par la matrice de rotation et le vecteur de translation du référentiel de la caméra / au référentiel de la scène 3D.

```math
[\mathbf{R}|\mathbf{t}]
```
+ __Matrice de la caméra pinhole ``\mathbf{P}``:__

```math
\mathbf{P}=\mathbf{K}[\mathbf{R}|\mathbf{t}]
```

+ __Effet de la distorsion radiale :__


  Les images formées par des dispositifs optiques à symétrie radiale tels que les lentilles introduisent une distorsion au niveau de l'image qui se manifeste par la courbature plus au moins prononcées des éléments droits de la scène.  
  
  Soit ``(x_c, y_c)`` les coordonnés d'un point obtenu après la transformation extrinsèque. Ils ont pour expression :
```math
\begin{align*}
x_c &= \frac{\mathbf{r}_x \cdot \mathbf{p} + t_x}{\mathbf{r}_z \cdot \mathbf{p} + t_z}\\
y_c &= \frac{\mathbf{r}_y \cdot \mathbf{p} + t_y}{\mathbf{r}_z \cdot \mathbf{p} + t_z}\\
\end{align*}
```

En adoptant un modèle quadratique de la distorsion, il devient aisé de compenser ses effets :
```math
\begin{align*}
\hat{x}_c &= x_c(1 + k_1 r_c^2 + k_2 r_c^4)\\
\hat{y}_c &= y_c (1 + k_1 r_c^2 + k_2 r_c^4)\\ 
\end{align*}
```
Avec ``r_c^2 = x_c^2 + y_c^2``
"""

# ╔═╡ 55f1f9b2-4ff8-46d1-9205-ba514050abd6
md"""
#### 2. Photométrie de l'image
$(load("cv_pics/cap08.png"))

##### 2.1 Illumination

Pour produire une image, la scène doit être éclairée par une source lumineuse primaire/secondaire, ponctuelle/étendue monochromatique ou polychromatique. Toute source lumineuse est caractérisée par la distribution spectrale de son intensité ``L(\lambda)``. Cette dernière décroît selon une loi quadratique avec la distance séparant la source de l'objet illuminé. 

##### 2.2 Diffusion et ombrage

Lorsque les objets d'une scène donnée sont éclairés par une certaine source lumineuse, ils diffusent la lumière reçue selon différent scenario selon la nature de leurs surfaces :
+ Réflexion ou diffusion spéculaire ;
+ Diffusion ;
+ Absorption.

Pour modéliser l'effet de la diffusion, on utilise la __BRDF__ (Fonction de distribution bidirectionnelle de la réflectance) et on se limite aux matériaux isotropes :
```math
fr(\mathbf{v}_i, \mathbf{v}_r , \mathbf{n}; \lambda) 
```
Pour calculer la quantitee de lumieure qui sera diffusee du point ``\mathbf{p}`` de la surface d'un objet dans la direction ``\mathbf{v}_r`` on utilise la formule suivante:

```math
L(\mathbf{v}_r,\lambda) = \int L(\mathbf{v}_i,\lambda) f_r(\mathbf{v_i}, \mathbf{v}_r, \mathbf{n}, \lambda) \cos^{+}(\theta_i)d\mathbf{v}_i
```
##### 2.3 L'optique de la caméra :
Principalement formée par des miroirs, des prismes et des lentilles. Pour comprendre leur fonctionnement, se référer au cours de l'optique géométrique.
"""

# ╔═╡ 55ceabfa-df49-4d68-bd14-c69a119052d0
md"""
#### 3. Images numeriques

Apres avoir fais le survole des differents elements composant le processus de formation de l'image nous aboutissons a notre finale destination. C'est la phase ou l'energie lumineuse capturee par la camera est transformee en information numerique.

$(load("cv_pics/cap09.png"))
 
##### 3.1 Capteurs electroniques (CCD/CMOS):

L'énergie lumineuse tombant sur les cellules photosensibles du capteur durant le temps d'exposition génèrent une quantité de charges qui sera transformées en tension électrique et convertit après en valeur numérique par un CAN de résolution binaire de 8, 10, 16 ou même 24 bits dans certaines cas. 

$(load("cv_pics/cap10.png"))


+ __Temps d'exposition :__ 

  Il influence aussi la qualité de la prise réalisée. Ainsi, un temps d'exposition très court élimine les effets de floutage généré par le mouvement rapide des éléments de la scène.


+ __Pas d'échantillonnage :__

  Il est en relation directe avec la résolution spatiale de l'image. De plus en plus, il est petit, mieux c'est pour la qualité de l'image et le niveau de détail qu'elle contient. Dans les capteurs courant, il est de l'ordre du micromètre. 


+ __Facteur de remplissage :__

  C'est le ratio entre la surface active de la cellule photosensible et sa taille théorique. Plus ce facteur s'approche de un plus l'énergie lumineuse capturée par l'optique de la caméra sera complètement convertis en charges électriques.


+ __Taille du capteur électronique :__


+ __Gain Analogique :__

  Dans les caméras vidéo, le gain de ces amplificateurs était traditionnellement contrôlé par l'AGC, qui ajustait ces valeurs pour obtenir une exposition globale de bonne qualité. Dans les appareils photo numériques récents, l'utilisateur dispose désormais d'un contrôle supplémentaire sur ce gain grâce au réglage ISO, qui est généralement exprimé en unités standard telles que 100, 200 ou 400 ISO. 
  
  Étant donné que le contrôle d'exposition automatisé de la plupart des appareils photo ajuste également l'ouverture et la vitesse d'obturation, le réglage manuel de l'ISO supprime un degré de liberté du contrôle de l'appareil photo, tout comme la spécification manuelle de l'ouverture et de la vitesse d'obturation. 

  En théorie, un gain plus élevé permet à l'appareil photo de mieux fonctionner dans des conditions de faible luminosité (moins de flou de mouvement dû à de longs temps d'exposition lorsque l'ouverture est déjà au maximum). Cependant, des réglages ISO plus élevés amplifient généralement le bruit du capteur

+ __Bruit du capteur :__

"""

# ╔═╡ 0843221c-97d5-4ed5-8163-5337d83530c0
md"""
##### 3.2 La couleur
+ __Perception naturelle :__

$(plot([load("cv_pics/cap11.png") load("cv_pics/cap12.png")], framestyle=:none))

###### 3.2.1 Classe des espaces de couleurs :
  
+ __Synthèse additive :__ (Superposition/Juxtaposition/altération temporelle) des couleurs primaires (R, G, B)
```math
\{C\} = r\,\{R\} + g\,\{G\} + b\,\{B\} \text{ avec } x =\int_0^\infty x(\lambda) L(\lambda)\,d\lambda\;, x=r,\,g,\, b 
```
$(load("cv_pics/cap13.png"))
   
+ __Synthèse soustractive :__ (Addition des couleurs secondaires (C, M, Y, K))
```math
\{C\} = c\,(\{N\} - \{R\}) + m\,(\{N\} - \{G\} + y\,(\{N\}-\{B\}) 
```
$(load("cv_pics/cap14.png"))

##### 3.3 Espace CIE RGB et XYZ

+ __Distribution spectrale des composantes tri-colorimétrique :__

$(load("cv_pics/cap15.png"))

+ __Matrice de passage ``RGB \rightarrow XYZ`` :__
```math
\begin{bmatrix}X\\Y\\Z\\ \end{bmatrix} ={1\over0.17697} 
\begin{bmatrix*}[l] 
0.49 & 0.31 & 0.20\\
0.17697 & 0.81240 & 0.01063\\
0.00 & 0.01 & 0.99 \\
\end{bmatrix*}
\begin{bmatrix} 
R \\ G \\ B\\
\end{bmatrix}
```
+ __Diagramme de chromacite :__
Les composantes trichomatiques de l'espace ``XYZ`` peuvent etre normalise par le biais de la relation suivante:

```math
x  = \frac{X}{X+Y+Z}\,, y  = \frac{Y}{X+Y+Z}\,, z  = \frac{Z}{X+Y+Z} \text{ tel que } x+y+z = 1
```

On appelle coordonnee de chromacitee le triplet ``(x,y,z)`` et vue que ces derniers ne sont pas independants on represente alors la couleur par deux coordonnnes seulement a savoir la pair (x,y) ce qui nou permet de tracer le diagramme de chromacite de l'espace ``XYZ`` 

$(load("cv_pics/cap16.png"))

##### 3.4 Espace de couleur L.a.b

Meme si le CIE XYZ permet de bien representer les couleurs naturelles tout en offrant la possibilite de separer les composantes radiometriques des composantes chromatiques mais il presente aussi certains defaut qui ont poussee le CIE pour developper un autre espace de couleur qui tient en compte des non-linearites du processus de perception de la couleur.

```math
L = 116f( {Y\over Y_n}) \text{ avec } Y_n \text{est la luminance du neutre.}
```
La fonction ``f`` est donnee par 
```math
 f(t) = \left\{\begin{matrix*}[l]
t^{1/3} & t > \delta^3\\
t/(3\delta^2) + 2\delta/3 & autrement\\
\end{matrix*}\right. \text{ avec } \delta = 6/29
```
Les composantes ``a`` et ``b`` sont donnees par les relations ci-apres:

```math
a = 500\left[f(X/X_n) -f(Y/Y_n)\right], \quad b = 200\left[f(Y/Y_n) -f(Z/Z_n)\right]
```
"""

# ╔═╡ a3c6a2a7-56cc-41b7-b64b-26648b3b39b1
md"""
+ __Exemples :__

Couleurs principales du sRGB : $([RGB{N0f8}(1,0,0) RGB{N0f8}(0,1,0) RGB{N0f8}(0,0,1)])

Couleurs secondaires du sRGB : $([RGB(1,1,0) RGB(0,1,1) RGB(1,0,1)])

Couleurs Principales du système Lab : $([Lab(50,0,0) Lab(0,0.5,0) Lab(0,0,0.5)])

Couleurs Principales du système HSV : $([HSV(0,100,100) HSV(120,100,100) HSV(240,0.25,1.0)])
"""

# ╔═╡ d3f46d9e-54f4-4c63-96e1-ca6239a96922
md"""
##### 3.5 Comment synthétiser la couleur ?
En théorie, les cellules photosensibles du capteur électronique doivent disposer de trois zones spectralement sensibles, l'une juxtaposée à côté de l'autre, pour capturer l'intensité lumineuse de chaque couleur principale :

```math
\begin{align*}
R &= \int_0^{\infty} L(\lambda) S_r(\lambda)\,d\lambda\\
G &= \int_0^{\infty} L(\lambda) S_g(\lambda)\,d\lambda\\
B &= \int_0^{\infty} L(\lambda) S_b(\lambda)\,d\lambda\\
\end{align*}
```
Dans le cas pratique, les capteurs disposent de cellules à large spectre, d'où l'impossibilité de générer directement les valeurs des composantes trichromatiques (R, G, B). Pour palier a ce problème la surface du capteur est couverte d'un filtre de couleur spéciale : Le filtre de Bayer.
$(load("cv_pics/cap17.png"))
"""

# ╔═╡ 8a4b167e-8c0b-4126-9172-0abdc9f2d52c
md" ## B. Representation des image"

# ╔═╡ 06c4ba8e-267e-4663-ab6e-a507bf13da54
md"""
Quel que soit le système de couleur retenu, une image est composée de pixels distribués uniformément sur une grille spatiale rectangulaire. Chaque pixel y est repéré par son indice de ligne et de colonne et caractérisé au moins par une valeur (dans le cas des images de luminescence), trois valeurs dans le cas des images naturelles, voir plus dans d'autres cas.

Dans ce qui suit, nous allons introduire les représentations d'images couramment utilisées en vision artificielle en mettant l'accent sur des frameworks spécifiques tel que :
- Image Processing toolbox de Matlab/octave;

- Pillow de Python;

- JuliaImages echosystem;

- OpenCV;
"""

# ╔═╡ a2b52599-4395-43d3-8532-76d7ad0cdbc3
md"""
### 1. Matlab/octave Image processing toolbox
Matlab/octave utilise un tableau multi-dimensionnel pour charger les pixels de l'image dans son environement de travail:

```matlab
image_file_path = "~/Documents/cv_pics/cap17.png";
im = imread(image-file_path);
size(im)
exit()
```
"""

# ╔═╡ 94184d6c-aaa8-47c5-810c-7f45035142c9
with_terminal() do 
	run(`flatpak run org.octave.Octave -qf scripts/im_data_struct.m`);
end

# ╔═╡ eed62ec4-1b91-4529-970b-ddf652917e62
md"### 2. Pillow Python image toolbox"

# ╔═╡ 73c3901a-ec96-4ed3-9b9f-49ee9bafbc10
begin
	IMG = pyimport("PIL.Image")
	np = pyimport("numpy");
	im = IMG.open("cv_pics/cap17.png");
	im_npy = np.array(im);
	py"""
	import PIL.Image as Image
	def f():
		with Image.open('cv_pics/cap17.png') as im:
			return list(im.getdata())[:10]
	"""
	with_terminal() do
		println("the size of the loaded image is : $(im.size)")
		println("Get 10 pixels from the object im instantiated from the class Image from Pillow:\n")
		for pix in py"f"()
			@show pix
		end
		println("\nEvery pixel is a tuple that contains the values of image canals at that pixel")
		println("\nSize of image after conversion to numpy.ndarray: $(size(im_npy))")
	end
end

# ╔═╡ 3b88fa8c-eb26-4847-8d33-c9e0bb26294e
md"### 3. ImagesJulia ToolBox"

# ╔═╡ cdb85e3d-8b85-4e54-aa41-2b126d7f59e2
let
	with_terminal() do
		im = load("cv_pics/cap17.png");
		println("The size of the loaded image is : $(size(im))")
		@show typeof(im)
		@show eltype(im)
		count = 0
		for pix in im
			count > 10 && break
			@show (pix.r, pix.g, pix.b)
			count+=1
		end
	end
end

# ╔═╡ 51335be9-2589-4ee3-8611-f9b476852d98
md"""
### 4. OpenCV ToolBox
#### 4.1 Implementation Python
"""

# ╔═╡ e59f8aa3-e61d-417d-b7cb-0f1910957af3
with_terminal() do
	println("the size of the image is:")
	code = """
	import cv2 as cv;
	im = cv.imread('cv_pics/cap17.png');
	print(f\"{im.shape}\\nThe type of the image is: {type(im)}\");
	print(f\"The type of the image pixels is: {im.dtype}\");
	nrows, ncols, ncanals = im.shape;
	for i in range(10):
		print(f"The pixel {i} is {im[i,1,:]}");
	""";
	run(`python -c $code`)
end

# ╔═╡ 468eec4f-8783-4575-824e-8b203c8fb96b
md"#### 4.2 Implementation c++ "

# ╔═╡ d5b12181-b245-42d8-883a-b9229ac7a705
with_terminal() do
	run(`/home/belkebir/Documents/Formations/OpenCV/prog1`)
end

# ╔═╡ 196159da-9dcf-483b-be89-e3f088d82869
# ╠═╡ disabled = true
#=╠═╡
Int.(rawview(channelview(im)[parse(Int,c),:,:]))
  ╠═╡ =#

# ╔═╡ 47cedbd3-2af0-43ba-bde9-ef89d32baef2
# ╠═╡ disabled = true
#=╠═╡
md"""
x = $(@bind x_ Slider(2:size(im,2)-1, default=120, show_value=true))
y = $(@bind y_ Slider(2:size(im,1)-1, default=170, show_value=true))

Type de Voisinage = $(@bind vtype TextField(default="v4"))
"""
  ╠═╡ =#

# ╔═╡ 03ecd740-5fed-4492-8c85-cead1955fc97
#=╠═╡
begin
	gr()
	if vtype == "v4"
		title = L"V_4"
		V = zeros(Bool,3,3)
		V[2,1:3] = [1,1,1]
		V[1:3,2] = [1,1,1]
		V = V .* im[y_-1:y_+1,x_-1:x_+1]
		V[2,2]=RGB(0,1,0)
	elseif vtype=="vd"
		title = L"V_D"
		V = ones(Bool,3,3)
		V[2,1:3] = zeros(Bool,3)
		V[1:3,2] = zeros(Bool,3)
		V = V .* im[y_-1:y_+1,x_-1:x_+1]
		V[2,2]=RGB(0,1,0)
	elseif vtype == "v8"
		title = L"V_8"
		V = ones(Bool,3,3)
		V = V .* im[y_-1:y_+1,x_-1:x_+1]
		V[2,2]=RGB(0,1,0)
	else
		V = zeros(RGB{N0f8}, 3, 3)
		title="Voisinage non reconue"
	end
	plot(V, t=:surface,title=title, framestyle=:none, size=(192,192))
end
  ╠═╡ =#

# ╔═╡ 9532bbfb-fa98-41f3-8f26-6a7473875529
# ╠═╡ disabled = true
#=╠═╡
md"""
#### Choose the distance between this :
1. Euclidienne
2. Manhattan
3. Echequier

x = $(@bind X Slider(2:size(im,2)-1, default=120, show_value=true))
y = $(@bind Y Slider(2:size(im,1)-1, default=170, show_value=true))

Type de distance = $(@bind dtype TextField(default="Euclidienne"))
"""
  ╠═╡ =#

# ╔═╡ c2576768-0411-4e27-ae4f-53bac5d46bf7
#=╠═╡
let
	M, N = size(im)
	dmap = zeros(M,N)
	if dtype == "Euclidienne"
		for i in 1:M
			for j in 1:N
				dmap[i,j]=sqrt((X - j)^2+(Y - i)^2)
			end
		end
	elseif dtype == "Manhattan"
		for i in 1:M
			for j in 1:N
				dmap[i,j]=abs(X - j)+abs(Y - i)
			end
		end
	elseif dtype == "Echequier"
		for i in 1:M
			for j in 1:N
				dmap[i,j]=max(abs(X - j),abs(Y - i))
			end
		end
	else
		nothing
	end
	dmap/=maximum(dmap)
	plot(im,framestyle=:none)
	plot!(dmap,t=:heatmap,c=:jet,framestyle=:none,alpha=0.55)
end
  ╠═╡ =#

# ╔═╡ 9ec90f0e-89a9-11eb-1597-471e29146d17
# ╠═╡ disabled = true
#=╠═╡
let
	if parse(Bool,filt_stat)
		kern = centered(ones(2*(step-1)+1,2*(step-1)+1))/(2*(step-1)+1)^2
		#kern = centered(ones(3,3))/9
		imf = imfilter(im, kern)
	else
		imf = im
	end
	im_decim = imf[1:step:end, 1:step:end]
	M, N = size(im_decim)
	y = collect(1:M)
	x = collect(1:N)
	plot(plot(im, st=:surface, title="Original"), plot(im_decim, st=:surface, title="Decimated"), framestyle=:none)
end
  ╠═╡ =#

# ╔═╡ 0afb1e68-89ac-11eb-31d1-11bed6098935
# ╠═╡ disabled = true
#=╠═╡
begin 
	Im = abs.(fftshift(fft(Float64.(channelview(im)[1,:,:]))))
	Im = Im/maximum(Im)
	u = range(-0.5, stop=0.5, length=N)
	v = range(-0.5, stop=0.5, length=M)
	plot(u, v, amp2db.(Im), st=:heatmap, c=:coolwarm, title="Spectre en [dB]")
end
  ╠═╡ =#

# ╔═╡ 5f53543a-89c5-11eb-394a-ff71e46bee7b
# ╠═╡ disabled = true
#=╠═╡
md"""
#### let displays some pixels
Line = $(@bind line Slider(1:M, default=50, show_value=true))
"""
  ╠═╡ =#

# ╔═╡ aa37cd62-89c4-11eb-36b5-07e8a328f929
#=╠═╡
begin
	md"""
	__Pixels code__ = $(eltype(im))
	
	$(im[line,75:81])
	
	__Color value of pixel at the line $(line) and coloumn 75 is__:
	* Red = $(im[line, 75].r)
	* Green = $(im[line, 75].g)
	* Blue = $(im[line, 75].b)
	"""

end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─bf73d046-36a6-46f3-9b62-b8815920671b
# ╟─e8cd9497-fa73-46bc-989e-8c1816b2a7c7
# ╟─adc8674a-89a2-11eb-2caf-1915cd540d9e
# ╟─b6409d51-9730-4fe3-b829-a13894103ba1
# ╟─fd1863f7-dd7c-4d9b-af44-8c9488be79cd
# ╟─915aee51-5be0-4ba2-81ea-15639f6b30c5
# ╟─14100cb9-52fa-425c-97a1-0477c4c70654
# ╟─5a116364-af9e-4571-ae6e-2df7450c92bb
# ╟─504b6ae5-fa6e-4be4-85da-e36e614514d8
# ╟─736fbc4f-4559-4131-8674-ecdd2ee8ce8f
# ╟─33f6b6be-cc82-476a-85da-2a6f98505a6a
# ╟─1aee84e1-eadc-48e2-a702-0b31ee16e0b2
# ╟─60358d7e-8cbc-4d53-9167-2ab36e660ceb
# ╟─2aff139c-324e-402c-826c-31f17af7a671
# ╟─9d84ac78-f151-49f8-9eea-6b9c09acc10f
# ╟─719bdf63-972b-4c3a-88f7-4c6fa0f7e054
# ╟─4146099c-7530-461f-967d-ca2b60ebda24
# ╟─1dc7eb8d-a04d-4f24-bd8d-97c6cf67a213
# ╟─b6257893-55f2-4434-ab23-04725d84ba39
# ╟─964f2417-f434-4390-ad03-f75fe7b1473d
# ╟─9bc6ead4-9b52-4e9d-8421-11ddf8546477
# ╟─55f1f9b2-4ff8-46d1-9205-ba514050abd6
# ╟─55ceabfa-df49-4d68-bd14-c69a119052d0
# ╟─0843221c-97d5-4ed5-8163-5337d83530c0
# ╟─a3c6a2a7-56cc-41b7-b64b-26648b3b39b1
# ╟─d3f46d9e-54f4-4c63-96e1-ca6239a96922
# ╟─8a4b167e-8c0b-4126-9172-0abdc9f2d52c
# ╟─06c4ba8e-267e-4663-ab6e-a507bf13da54
# ╟─a2b52599-4395-43d3-8532-76d7ad0cdbc3
# ╟─94184d6c-aaa8-47c5-810c-7f45035142c9
# ╟─6da4b103-e1c5-465b-aa92-47aa2cab43a2
# ╟─eed62ec4-1b91-4529-970b-ddf652917e62
# ╟─73c3901a-ec96-4ed3-9b9f-49ee9bafbc10
# ╟─3b88fa8c-eb26-4847-8d33-c9e0bb26294e
# ╟─cdb85e3d-8b85-4e54-aa41-2b126d7f59e2
# ╟─51335be9-2589-4ee3-8611-f9b476852d98
# ╟─e59f8aa3-e61d-417d-b7cb-0f1910957af3
# ╟─468eec4f-8783-4575-824e-8b203c8fb96b
# ╟─d5b12181-b245-42d8-883a-b9229ac7a705
# ╟─196159da-9dcf-483b-be89-e3f088d82869
# ╟─47cedbd3-2af0-43ba-bde9-ef89d32baef2
# ╟─03ecd740-5fed-4492-8c85-cead1955fc97
# ╟─9532bbfb-fa98-41f3-8f26-6a7473875529
# ╟─c2576768-0411-4e27-ae4f-53bac5d46bf7
# ╟─9ec90f0e-89a9-11eb-1597-471e29146d17
# ╟─d048d350-89ab-11eb-2ff5-038da6d710e3
# ╟─0afb1e68-89ac-11eb-31d1-11bed6098935
# ╟─5f53543a-89c5-11eb-394a-ff71e46bee7b
# ╟─aa37cd62-89c4-11eb-36b5-07e8a328f929
