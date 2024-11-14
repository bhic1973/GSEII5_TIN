# Génie des Systèmes Embarqués et Informatique Industrielle S5
## Traitement d'image
### Procédure pour visualiser les notesbooks contenant le support de cours interactif
1. **Installer Julia LTS 1.10.6**
   Pour Installer Julia dans sa dernière version LTS utiliser les liens suivants:
   - Linux :
     ```bash
     wget -r https://julialang-s3.julialang.org/bin/linux/aarch64/1.10/julia-1.10.6-linux-aarch64.tar.gz
     ```
   - Windows
     [Julia Installer link](https://julialang-s3.julialang.org/bin/winnt/x64/1.10/julia-1.10.6-win64.exe)

2. Une fois L'installation est terminé, créer un répertoire de travail (Ex: TIN) et tapez les commandes suivantes dans le terminal :
   ```cmd
   cd TIN
   julia --project
   ```
   Dans le REPL de julia executez les instructions suivantes :
   ```julia
   ]
   activate .
   instantiate
   ```
   *NB: Vérifier que vous avez télécharger les fichiers Manifest.toml et Project.toml depuis ce repo avant de lancer les commandes ci-dessus.*
3. Une fois c'est terminé, Télécharger les notebooks en locale et placez les dans le repertoire TIN.
4. Relancer julia
5. Executez les commandes suivantes:
   ```julia
   using Pluto
   Pluto.run()
   ```
6. Votre navigateur va se lancer et la fenêtre de Pluto va apparaître que vous pouvez maintenant manipuler sans aucune difficulté. 
