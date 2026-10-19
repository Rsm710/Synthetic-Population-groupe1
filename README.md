# Population synthétique : greedy colonne par colonne

On reconstruit N individus décrits par K attributs catégoriels. La population doit respecter **exactement** des
tables de marges, avec le moins possible de clones (individus identiques). Tout est lu depuis un CSV
`cell_id,table_id,arity,scope,cell,target` : K, N, les domaines, les tables et leurs arités.

La méthode est un **greedy** : une construction, puis un post-traitement local qui garde les marges exactes. CP-SAT
sert uniquement de référence, pour connaître l'optimum.

## En bref

Instance `constraint_cells.csv` : N = 500 individus, K = 12 attributs à 3 valeurs, 45 tables d'arité 1 à 4.

| | Personas distinctes (moy. / meilleure) | Temps / run |
|---|---|---|
| **défaut** | **302,5 / 307** | 2,2 s |
| défaut, dé-clonage 6 s | 307,9 / 311 | 6,6 s |
| optimum prouvé (CP-SAT) | 331 | 900 s |

- Tous les runs sont exacts : erreur nulle sur toutes les cellules, revérifiée par `Evaluation`.
- On atteint **91 % de l'optimum en moyenne** et 94 % au meilleur run. Les marges imposent à elles seules 169 copies
  en trop.
- Les ressemblances entre individus restent proches d'une population réaliste (voir
  [Référence de réalisme](#référence-de-réalisme)).

## Démarrage rapide (PowerShell, JDK 21)

```powershell
javac -encoding UTF-8 -d out (Get-ChildItem src\*.java).FullName
java -cp out GreedyMain constraint_cells.csv --runs 10                 # écrit population.csv (meilleur run)
java -cp out GreedyMain constraint_cells.csv --runs 20 --quiet         # statistiques sur 20 runs
java -cp out CheckPopulation constraint_cells.csv population.csv       # revérifie un fichier écrit
```

Pour des accents corrects dans la console : `[Console]::OutputEncoding = [Text.Encoding]::UTF8` et
`java -Dstdout.encoding=UTF-8 ...`.

Autres commandes :

```powershell
java -cp out GreedyMain constraint_cells.csv --runs 20 --quiet --declone 0   # construction seule
java -cp out GreedyMain constraint_cells.csv --runs 20 --quiet --deg4 0      # sans doubles échanges
java -cp out Deg4Diag constraint_cells.csv 1 population.csv                  # mouvements de degré 4 améliorants
java -cp out Bounds constraint_cells.csv population.csv                      # bornes par propagation
java -cp out Main constraint_cells.csv --runs 2 --time 8000                  # ancienne référence

# optimum exact (Python 3.12 : pip install ortools pandas) ; '-' au lieu de population.csv = sans solution de départ
python tools\reference_cpsat.py constraint_cells.csv 900 population.csv cpsat_population.csv
# référence de réalisme (numpy) : échantillon de N individus tiré selon P*
python tools\reference_ipf.py constraint_cells.csv
```

## Fichiers

| Fichier | Rôle |
|---|---|
| `src/GreedyMain.java` | point d'entrée : runs, statistiques, écriture de `population.csv` |
| `src/ColumnGreedy.java` | construction colonne par colonne (greedy 1, par défaut) |
| `src/ColumnGreedy2.java` | variante : places uniques d'abord, doublons concentrés (`--greedy 2`) |
| `src/Declone.java` | dé-clonage : échanges qui gardent les marges exactes |
| `src/Evaluation.java` | vérification indépendante (seul juge des résultats) et métriques |
| `src/CheckPopulation.java` | revérifie un fichier population |
| `src/Deg4Diag.java`, `Bounds.java`, `Lns.java` | outils d'analyse (non utilisés par `GreedyMain`) |
| `src/Main.java`, `GreedyBuilder.java`, `LocalSearch.java`, `Population.java` | ancienne référence : greedy individuel + recuit |
| `tools/reference_cpsat.py` | optimum exact par CP-SAT |
| `tools/reference_ipf.py` | référence de réalisme : loi d'entropie maximale P* (IPF) |

## Méthode

Défaut : `--entropy full --dentropy loo --deg4 0.2`. Un run = construction, puis 1,5 s de dé-clonage.

Vocabulaire :
- une **persona** est une combinaison de valeurs des K attributs ;
- deux individus sont **jumeaux** s'ils sont identiques sur les colonnes déjà fixées ;
- la **distance** entre deux individus est le nombre d'attributs où ils diffèrent (distance 0 : clones ; distance 1 :
  quasi-clones).

### 1. Construction (`ColumnGreedy`)

On fixe un attribut à la fois, pour toute la population. Une colonne exacte n'est jamais cassée ensuite.

1. **Ordre des attributs.** On prend d'abord les étapes « parfaites » : parmi les projections connues des tables de
   l'attribut, l'une contient toutes les autres, donc la colonne est toujours faisable. Sinon, on prend l'attribut
   le plus relié aux colonnes déjà fixées.
2. **Remplir la colonne a.** Chaque table qui contient a impose des comptes exacts sur sa projection (attributs connus
   et a). On traite d'abord les individus les plus contraints. Une valeur est possible si :
   - elle a du résiduel dans toutes ses cellules ;
   - elle respecte le **look-ahead entier** : pour chaque attribut futur c, un groupe d'individus (mêmes valeurs sur
     les attributs liés à c) contient au plus cap(g) = Σ_w min_t cible_t(g, w) individus. C'est la borne u_x du cours
     sur les personas partielles.
3. **Choisir parmi les valeurs possibles.** On sépare les jumeaux : la valeur v est préférée si elle a reçu peu de
   jumeaux de l'individu. En mode `--entropy full` (défaut), le score est 5·10⁴ − δ(n_{g,v}) + α ln(mr), avec :
   - n_{g,v} le nombre de jumeaux déjà envoyés sur v, et δ(n) = (n+1) ln(n+1) − n ln n ;
   - mr le plus petit résiduel des cellules de v ;
   - α = 0,3.

   Chaque colonne maximise ainsi de façon gloutonne l'entropie conditionnelle H(X_a | X_connus). Le terme α ln(mr) est
   indispensable : sans lui, le choix est myope et force des paquets aux colonnes suivantes.
4. **En cas d'impasse.** On fait une micro-réparation de la colonne, puis de nouveaux tirages. Ensuite, on revient en
   arrière sur 1 puis 2 colonnes déjà fixées (les attributs des tables violées). Enfin, on relance toute la
   construction, au plus 10 fois.

### 2. Dé-clonage (`Declone`)

On part d'une population exacte et on applique des mouvements qui **gardent toutes les marges exactes**.

**Échange simple.** Deux individus i et j échangent leurs valeurs sur un ensemble d'attributs S. Ce n'est permis que
s'ils coïncident sur la frontière B(S) : les attributs hors de S qui partagent une table avec S. Chaque table est alors
incluse dans S ∪ B(S) ou disjointe de S, donc ses comptes ne changent pas.

Ce voisinage contient **tous** les mouvements à deux individus. En effet, si i, j deviennent i′, j′ avec les mêmes
projections sur toutes les tables, alors sur chaque table i′ prend ses valeurs entièrement chez i ou entièrement
chez j. Mais le graphe des tables est dense : seuls 85 ensembles S sur 2047 laissent un attribut hors de S ∪ B(S). Ce
voisinage plafonne donc vers 292 distinctes.

**Double échange compensé (`--deg4`).** C'est un mouvement à 4 individus :
- (i, j) échangent S alors qu'ils diffèrent sur un attribut f de la frontière (plus généralement un ensemble F,
  1 ≤ |F| ≤ `--dfmax`). Seul, ce mouvement casse les tables Tf qui touchent F et S.
- Une seconde paire (k, l) échange le même S et annule exactement l'écart. Avec U = S ∩ ∪Tf et C = ∪Tf \ (S ∪ F), il
  faut :
  - k = (valeurs de i sur F, de j sur U, de i sur C) ;
  - l = (valeurs de j sur F, de i sur U, de i sur C) ;
  - k et l égaux sur B(S) \ F.

Pour chaque table de Tf, les cellules ajoutées par (i, j) sont celles que retire (k, l), et inversement. Toute autre
table est dans S ∪ (B(S) \ F) ou disjointe de S pour chaque paire. Ces mouvements existent parce que le graphe des
tables a des cycles qu'aucune table ne couvre (A3–A9–A11, A0–A2–A4…). Le modèle n'est pas décomposable : sa base de
Markov contient des mouvements de degré ≥ 4.

**Recherche.** On tire un individu i en doublon et un S :
- échange simple : on prend le meilleur partenaire j ;
- double échange (probabilité 0,2) : on essaie jusqu'à 4 partenaires j, puis on cherche toutes les paires (k, l) en un
  passage sur la population.

Les gains positifs et nuls sont acceptés, et on garde la meilleure population rencontrée.

**Score lexicographique.** On maximise d'abord le nombre de distinctes, puis, à égalité, Σ_b H(X_{−b}) (`--dentropy
loo`). Ici X_{−b} est la population privée de l'attribut b, et H(X_A) = ln N − (1/N) Σ_c n_c ln n_c, en nats.
- H(X_{−b}) = ln N si et seulement si tous les individus restent distincts quand on ignore b.
- Le déficit D1 = K·ln N − Σ_b H(X_{−b}) est donc nul si et seulement si aucune paire n'est à distance ≤ 1.
- Les comptes sont tenus à jour de façon incrémentale.

### Pourquoi l'entropie

Soit P̂ une loi empirique qui respecte exactement les tables, et P* la loi d'entropie maximale sous ces tables (IPF).
Alors KL(P̂ ‖ P*) = H(P*) − H(P̂) : maximiser H(P̂) rapproche de la loi la moins structurée compatible avec les marges.
Les entropies des projections incluses dans une table sont imposées ; seules les projections non couvertes peuvent
changer.

### Mesures (`Evaluation`)

`Evaluation.evaluate` recalcule tout depuis zéro. Il donne :
- l'erreur L1 et le nombre de personas distinctes ;
- H(X) et les déficits D1 (rayon 1) et D2 (rayon 2) ;
- Σ_b H(X_b | X_{−b}) ;
- les paires d'individus à distance 0, 1 et 2 ;
- Σ I(X_a ; X_b) sur les 19 paires d'attributs couvertes par aucune table.

## Résultats

20 runs (seeds 0–19), dé-clonage 1,5 s sauf mention contraire, tous exacts.

| Configuration | Distinctes moy. | Meilleure | Temps / run |
|---|---|---|---|
| construction seule (greedy 1) | 251,9 | 266 | 0,7 s |
| construction seule (greedy 2) | 258,8 | 270 | 1,1 s |
| + dé-clonage par échanges (`--entropy off --dentropy off --deg4 0`) | 292,5 | 297 | 2,6 s |
| + entropie (`--entropy full --dentropy loo --deg4 0`) | 292,8 | 299 | 2,2 s |
| **+ doubles échanges compensés (défaut)** | **302,5** | 307 | 2,2 s |
| défaut, dé-clonage 6 s (10 runs) | 307,9 | 311 | 6,6 s |
| **optimum (CP-SAT, prouvé)** | **331** | 331 | 900 s |

### Comparaison avec l'optimum

| | Distinctes | Copies en trop | Personas en doublon | Individus dans un doublon | Multiplicité max |
|---|---|---|---|---|---|
| meilleur run sans doubles échanges | 302 | 198 | 29 | 227 | 48 |
| meilleur run, défaut à 6 s | 311 | 189 | 30 | 219 | 37 |
| optimum CP-SAT | **331** | **169** | 10 | 179 | 49 |

- CP-SAT (8 threads, partant d'une population greedy à 302) trouve 330 en 300 s, puis prouve l'optimum 331 en 900 s.
  La population obtenue est vérifiée exacte par `CheckPopulation`.
- L'optimum concentre ses doublons sur 10 personas, dont deux gros paquets (×49 et ×46). Les 321 autres individus y
  sont en exemplaire unique. Le greedy répartit ses doublons sur environ 30 personas : les copies évitables sont dans
  les petits paquets (×2 à ×7).

### Référence de réalisme

**Faut-il réduire les quasi-clones pour être plus réaliste ? Non.** On les calibre sur une population typique, on ne
les minimise pas.

- **Les marges fixent la distance moyenne.** La proportion de paires qui diffèrent sur un attribut ne dépend que de sa
  marge. La distance moyenne entre deux individus vaut donc **4,57 sur 12** pour toute population exacte. Quatre
  attributs sont presque constants (A8 : 98 % sur une valeur, A0 : 92 %, A6 : 91 %, A10 : 90 %). Retirer des paires
  proches en crée forcément plus loin.
- **Une population réelle contient des clones.** Si les marges viennent de 500 vraies personnes, cette population est
  exacte : elle a donc au moins 169 copies en trop.
- **Échantillon typique.** `tools/reference_ipf.py` calcule P* par IPF sur les 3¹² personas, puis l'espérance des
  mesures pour 500 individus tirés selon P*.

| | Distinctes | Paires d = 0 | d = 1 | **d ≤ 1** | d = 2 | Voisins à d ≤ 1 par individu |
|---|---|---|---|---|---|---|
| **échantillon typique de P*** | 232 | 1935 | 7151 | **9086** | 13420 | 36,3 |
| `--dentropy off --deg4 0` | 292,7 | 2200 | 7167 | 9367 | | 37,5 |
| `--dentropy loo --deg4 0` | 292,8 | 1980 | 7066 | 9046 | 13073 | 36,2 |
| **défaut (`--deg4 0.2`)** | 302,5 | 1849 | 6866 | **8714** | 13271 | 34,9 |
| `--dmu 0 --deg4 0` | 276,9 | 1631 | 6893 | 8524 | 13396 | 34,1 |

- Sans doubles échanges, le niveau de ressemblance est celui d'un échantillon typique (9046 contre 9086 paires).
- Le défaut est 4 % en dessous. C'est une conséquence des distinctes gagnées, pas un objectif ; c'est accepté.
- `--dmu 0` vise les ressemblances et perd 16 distinctes : il ne sert que si l'on cherche autre chose que le réalisme
  (confidentialité, diversité imposée).

## Options de `GreedyMain`

Défauts entre parenthèses.

| Groupe | Options |
|---|---|
| construction | `--greedy` (1), `--beta` (2.0, séparation des jumeaux), `--gamma` (1.0, concentration du greedy 2), `--entropy off\|full\|loo` (full), `--alpha` (0.3, poids de ln(mr)), `--lambda` (1.0, poids des H(X_{−b}) en mode loo) |
| faisabilité | `--fc on\|off` (on, look-ahead), `--retries` (3), `--backtrack` (40), `--depth` (2), `--restarts` (10), `--repair` (500), `--perfect` (1000), `--cross` (0) |
| dé-clonage | `--declone` (1500 ms, 0 = désactivé), `--dentropy off\|full\|loo\|near` (loo), `--deg4` (0.2, probabilité d'un double échange), `--dfmax` (1, \|F\| maximal), `--hguard on\|off` (off, refuse tout mouvement qui baisse l'entropie), `--neutral` (1.0, probabilité d'accepter un gain nul) |
| dé-clonage, variantes testées | `--dmu` (−1 = lexicographique ; μ ≥ 0 : μ × distinctes + N × entropie), `--dradius` (1), `--cycle` (0), `--temp` (0, recuit), `--concentrate on\|off` (off), `--small` (0) |
| dé-clonage après chaque colonne | `--coldeclone` (0 ms = désactivé), `--colfrom` (1), `--colneutral` (1.0), `--colspread on\|off` (off), `--colbound on\|off` (off) |
| runs | `--runs`, `--seed`, `--starts` (1 construction par run), `--out` (population.csv), `--quiet` |

## Journal des expériences

Chaque essai commence par sa conclusion. Sauf mention contraire : 20 runs, seeds 0–19, tous exacts. « Défaut » désigne
la configuration de l'époque.

### Construction : historique du greedy 1

**Conclusion : 100 % de runs exacts grâce au look-ahead, au retour arrière sur 2 colonnes et aux relances.**

| Version | Exacts | Temps / run | Distinctes |
|---|---|---|---|
| ColumnGreedy initial | 2/20 | 0,1–1,8 s | ~240 |
| + retour arrière ciblé | 62,5 % | 0,8 s | 244 |
| + look-ahead par capacités | 85,5 % | 1,5 s | 244 |
| + ordre « étapes parfaites », 40 retours arrière | 293/300 (97,7 %) | 0,76 s (max 5,9 s) | 246 |
| + retour arrière sur 2 colonnes, sans relance | 200/200 | 0,47 s (max 5,7 s) | 246,5 |
| + relances | 600/600, dont 596 du premier coup | 0,70 s (max 8,4 s) | 246,3 |
| beta = 2 | 40/40 | 1,1 s | 251,1 |
| + dé-clonage 1,5 s | 40/40 | 2,3 s (max 7,0 s) | 293,0 (max 302) |

### Construction : greedy 2 (`--greedy 2`)

**Conclusion : meilleur sur la construction seule (+6,9), moins bon après dé-clonage (−2,3).**

Une branche est une persona partielle (attributs connus, a = v). Son potentiel est le nombre de personas complètes
actives qui la prolongent : la branche ne peut pas contenir plus de personas distinctes que cela.
- Place libre (potentiel > jumeaux déjà envoyés) : on sépare les jumeaux comme le greedy 1.
- Branche pleine : l'individu sera un doublon, on le met avec le plus de jumeaux (mr × (1 + tw)^gamma).

| gamma | 0,5 | **1** | 2 |
|---|---|---|---|
| Distinctes (sans dé-clonage) | 256,7 | **258,8** | 255,6 |

Une première version, qui préférait toujours les branches à fort potentiel, faisait 243,7 : au début, toutes les
branches ont un potentiel de plusieurs centaines. Hypothèse non vérifiée pour la perte après dé-clonage : les
doublons concentrés laissent moins d'échanges utiles.

### Construction : entropie (`--entropy`)

**Conclusion : `full` avec α = 0,3 donne +5 distinctes reproductibles ; le gain sur H(X) est faible (+0,01 nat).**

Sans dé-clonage. `loo` ajoute au score λ Σ_b δ(n_{−b,v}), le compte de la cellule de X_{−b}, pour pénaliser les
quasi-clones.

| Variante | Distinctes | H(X) | Σ H(X_b \| X_−b) | Paires d = 1 | Σ I paires libres |
|---|---|---|---|---|---|
| off (beta = 2) | 251,9 | 4,943 | 2,478 | 7131 | 0,372 |
| full, α = 0 | 240,4 | 4,863 | 2,320 | 7245 | 0,375 |
| full, α = 0,1 | 248,3 | 4,954 | 2,497 | 6929 | 0,353 |
| **full, α = 0,3** | **256,9** | **4,959** | 2,505 | 7056 | 0,357 |
| full, α = 0,5 / 1 / 3 | 258,1 / 256,7 / 245,6 | 4,948 / 4,924 / 4,867 | | | |
| loo, α = 0, λ = 1 | 240,1 | 4,865 | 2,272 | 7167 | 0,355 |
| loo, α = 0,3, λ = 0,03 / 0,1 / 1 | 251,9 / 247,8 / 242,0 | 4,950 / 4,929 / 4,881 | 2,484 / 2,423 / 2,285 | ~7045 | 0,35–0,36 |
| loo, α = 1, λ = 0,1 / 0,3 / 1 / 3 | 257,1 / 253,3 / 244,6 / 242,0 | 4,935 / 4,938 / 4,915 / 4,872 | | | |

- Le gain glouton seul (α = 0) est myope : il étale chaque groupe sans regarder les marges des cellules.
- Contrôle sur les seeds 100–139 : off H = 4,947 ± 0,029 ; full α = 0,3 : 4,950 ± 0,024 (+4,7 distinctes).
- Après dé-clonage : off 292,5 distinctes (H = 4,994) ; full α = 0,3 : 294,0 (H = 5,002).
- `loo` baisse Σ H(X_b | X_−b), mais surtout parce que H(X) baisse. Σ_b H(X_−b) est meilleur avec `full` (57,0)
  qu'avec `loo` (56,6–56,8).
- Les paires à distance 1 (~7000) sont dominées par les gros paquets imposés : un paquet ×48 voisin d'un ×49 en donne
  2352 à lui seul.

### Dé-clonage : premiers réglages

**Conclusion : le dé-clonage par échanges plafonne en environ 1 s, quelle que soit la stratégie d'acceptation.**

- **Budget** : beta = 1 donne 284,5 en 1 s et 285,3 en 8 s.
- **Mouvements neutres** : indispensables (277 sans eux, contre 285).
- **Recuit** (6 runs, 4 s) : T₀ = 0,3 / 0,6 / 1 / 2 donne 293,0 / 293,5 / 293,7 / 293,0. Aucun gain.
- **Multi-départ** (6 runs) : 1 × 8 s → 293,2 ; 4 × 2 s → 295,3 ; 8 × 1 s → 297,0 (15 s par run).
- **Viser la structure de l'optimum** (8 runs) : la structure se rapproche (20 personas en doublon au lieu de 35), mais
  pas le nombre de distinctes.

| Variante | Distinctes moy. | Meilleure | Personas en doublon (meilleure) |
|---|---|---|---|
| aucune | 293,0 | 296 | 35 |
| `--concentrate on` | 293,4 | 297 | 27 |
| `--small 2` | 293,9 | 297 | 23 |
| `--concentrate on --small 2` | 293,8 | 297 | 20 |

### Dé-clonage après chaque colonne (`--coldeclone`)

**Conclusion : +22 distinctes sur la construction, mais aucun gain après le dé-clonage final. Désactivé.**

Après chaque colonne exacte, on dé-clone la population partielle, restreinte aux attributs connus C : S ⊆ C, et i, j
coïncident sur B(S) ∩ C. Chaque table garde les comptes de sa projection sur C, donc la suite reste faisable. Avec
`--colbound on`, l'objectif est Σ_g min(m_g, pot(g)), avec pot(g) le potentiel du greedy 2.

| Variante | Runs | Distinctes moy. | Meilleure | Temps / run |
|---|---|---|---|---|
| greedy 1 seul | 20 | 251,9 | | 0,9 s |
| greedy 1 + `--coldeclone 100` | 20 | **274,5** | | 1,3 s |
| greedy 1 + `--coldeclone 100 --colspread on` | 20 | 274,0 | | 1,5 s |
| greedy 1 + `--coldeclone 100 --colbound on` | 20 | 258,3 | 268 | 1,7 s |
| greedy 2 + `--coldeclone 100` | 20 | **277,1** | 289 | 2,3 s |
| greedy 2 + `--coldeclone 100 --colbound on` | 20 | 264,4 | 277 | 2,1 s |
| **avec dé-clonage final 1,5 s** | | | | |
| sans (seeds 100–139) | 40 | 292,2 | 302 | 2,5 s |
| `--coldeclone 100` (seeds 100–139) | 40 | 292,2 | 300 | 3,2 s |
| `--coldeclone 100 --colbound on` (seeds 100–139) | 40 | 292,5 | 299 | 3,8 s |
| `--coldeclone 30 / 100 / 300` | 20 | 290,1 / 292,8 / 291,9 | | 2,8 / 3,1 / 4,4 s |
| `--coldeclone 100 --colfrom 5 / 8` | 20 | 292,8 / 291,9 | 299 / 297 | 3,0 s |
| `--coldeclone 100 --colneutral 0.3` | 20 | 291,0 | | 3,0 s |
| `--coldeclone 100 --colspread on` | 20 | 293,0 | | 3,1 s |
| greedy 2 + `--coldeclone 100` (± `--colbound on`) | 20 | 288,9 / 292,4 | 294 / 299 | 3,5 s |

Le dé-clonage final atteint le même plateau (~292) que la construction donne 251 ou 275 distinctes. Les copies
évitables ne viennent pas de la construction, mais de la rigidité du voisinage d'échanges.

### Dé-clonage : critère d'entropie (`--dentropy`)

**Conclusion : `--entropy full --dentropy loo` retire 317 paires à distance ≤ 1 (−3,4 %) sans perdre de distinctes.**

À distinctes égales, le score ajoute un gain d'entropie :
- `full` maximise H(X) ;
- `loo` maximise Σ_b H(X_{−b}) ;
- `near` minimise Σ_b H(X_b | X_{−b}), qui ne compte que la distance exactement 1.

| Construction | `--dentropy` | Distinctes | H(X) | Σ H(X_b \| X_−b) | Paires d = 1 | Σ I paires libres |
|---|---|---|---|---|---|---|
| off | off | 292,5 | 4,993 | 2,607 | 7179 | 0,392 |
| off | full | 292,4 | 5,024 | 2,683 | 7156 | 0,380 |
| off | loo | 291,4 | 5,023 | 2,667 | 7136 | 0,374 |
| off | near | 291,4 | 4,944 | **2,424** | 6983 | 0,388 |
| **full (α = 0,3)** | **full** | **293,4** | **5,036** | 2,716 | 7100 | **0,370** |
| loo (α = 0,3, λ = 0,1) | near | 290,5 | 4,941 | 2,426 | **6903** | 0,383 |

En comptant aussi les clones :

| Construction / `--dentropy` | Distinctes | Paires d = 0 | d = 1 | d ≤ 1 |
|---|---|---|---|---|
| off / off | 292,7 | 2200 | 7167 | 9367 |
| off / near | 291,9 | 2408 | 6982 | 9390 |
| off / full | 292,4 | 2050 | 7166 | 9215 |
| off / loo | 291,4 | 2037 | 7121 | 9157 |
| loo (λ = 0,1) / loo | 290,1 | 2005 | 7103 | 9108 |
| **full / loo** | **292,5** | **1975** | **7075** | **9050** |

- À distinctes égales, H(X) et les quasi-clones s'opposent : monter H étale les doublons en personas voisines.
- `near` retire 185 paires à distance 1 mais ajoute 208 clones : il transforme des quasi-clones en clones. Ce n'est pas
  le bon objectif.
- Les paires de clones viennent surtout des deux gros paquets imposés (×46 à ×49, plus de 1000 paires chacun).

### Dé-clonage : objectif pondéré, rayon 2, cycles

**Conclusion : seul l'objectif pondéré réduit les ressemblances, au prix de distinctes. Le rayon 2 et les cycles
n'apportent rien.**

- **Objectif pondéré** (`--dmu μ`) : on maximise μ × distinctes + N × Σ_B H(X_{−B}). Avec μ petit, on accepte de perdre
  des distinctes pour retirer des ressemblances.
- **Rayon 2** (`--dradius 2`) : on utilise Σ_{|B|=2} H(X_{−B}) sur les 66 paires d'attributs, maximal si et seulement
  si aucune paire n'est à distance ≤ 2.
- **Cycles** (`--cycle p`) : i ← j ← k ← i sur S, les trois individus égaux sur B(S). Un cycle est la composée de deux
  échanges valides : il n'élargit pas le voisinage.

| `Declone` (construction `--entropy full`) | Distinctes | D1 | Paires d = 0 | d = 1 | **d ≤ 1** | d = 2 | D2 |
|---|---|---|---|---|---|---|---|
| lexicographique, r = 1 | **292,8** | 16,87 | 1980 | 7066 | 9046 | 13073 | 109,37 |
| `--dmu 30` | 287,1 | 16,76 | 1881 | 6976 | 8856 | 13193 | 108,89 |
| `--dmu 10` | 285,2 | 16,67 | 1796 | 6956 | 8753 | 13293 | 108,51 |
| `--dmu 3` | 280,2 | 16,60 | 1693 | 6932 | 8624 | 13357 | 108,17 |
| `--dmu 0` | 276,9 | **16,57** | **1631** | 6893 | **8524** | 13396 | **108,03** |
| `--dradius 2` | 290,2 | 16,90 | 1971 | 7060 | 9031 | 13066 | 109,50 |
| `--dradius 2 --dmu 0` | 276,1 | 16,60 | 1634 | **6846** | **8480** | 13419 | 108,05 |
| `--cycle 0.3` / `0.6` | 290,5 / 287,2 | 16,87 / 16,89 | | | 9005 / 8961 | | 109,40 / 109,47 |

Budget 6 s (10 runs) : lexicographique 293,5 distinctes, d ≤ 1 = 9051 ; `--dradius 2` 293,0 / 9028 ;
`--cycle 0.3` 292,8 / 9023 ; `--dmu 0` 276,7 / 8494 ; `--dmu 0 --dradius 2` 275,8 / 8431.

- De μ = −1 à μ = 0 : paires à distance ≤ 1 −6 %, paires de clones −18 %, mais −16 distinctes.
- Rayon 2 : la somme des distances sur toutes les paires est fixée par les marges. Retirer des paires à distance ≤ 1 en
  crée surtout à distance 2. L'objectif reste dominé par d ≤ 1 (une paire de clones compte dans les 66 projections, une
  paire à distance 2 dans une seule).

### Dé-clonage : doubles échanges compensés (`--deg4`)

**Conclusion : +10 distinctes à budget égal, +14 à 6 s, et l'entropie augmente aussi. Retenu par défaut.**

**Diagnostic préalable** (`Deg4Diag`, |F| = 1, sur deux populations finales de `--deg4 0`). L'outil parcourt toutes les
paires de personas et tous les S. L'écart de chaque paire est haché linéairement (hachage de Zobrist additif), ce qui
permet de tester δ₁ + δ₂ = 0 en O(1). Chaque mouvement trouvé est revérifié table par table : aucun faux positif.

| Population | Échanges valides | dont gain de distinctes | Mouvements de degré 4 | dont gain de distinctes | dont entropie en hausse à distinctes égales |
|---|---|---|---|---|---|
| seed 0 (293 distinctes) | 17 835 | **0** | 309 443 | **137** (max +2) | 1870 |
| seed 1 (296 distinctes) | 19 137 | **0** | 371 707 | **211** | 2403 |

La population est un optimum local pour les échanges, mais pas pour le degré 4.

**Réglages :**

| `--deg4` | Distinctes | H(X) | D1 | Paires d = 0 | d = 1 | **d ≤ 1** | d = 2 | Temps / run |
|---|---|---|---|---|---|---|---|---|
| 0 | 292,6 | 5,033 | 16,87 | 1988 | 7080 | 9068 | 13063 | 2,4 s |
| 0,05 | 302,2 | 5,094 | 16,26 | 1851 | 6876 | 8727 | 13294 | 2,2 s |
| 0,1 | 302,1 | 5,092 | 16,28 | 1860 | 6849 | 8709 | 13313 | 2,2 s |
| **0,2** | **302,5** | 5,097 | 16,22 | 1849 | 6866 | 8714 | 13271 | 2,2 s |
| 0,3 | 300,9 | 5,096 | 16,23 | 1821 | 6848 | 8669 | 13300 | 2,2 s |
| 0,5 | 301,5 | **5,105** | **16,15** | 1800 | 6787 | 8587 | 13360 | 2,3 s |
| 0,1 / 0,2 avec `--dfmax 2` | 301,7 / 302,6 | 5,087 / 5,098 | 16,33 / 16,23 | | | 8804 / 8739 | | 2,2 s |
| 0,1 avec `--hguard on` | 293,3 | 5,129 | 15,93 | 1570 | 6745 | 8315 | 13530 | 2,2 s |
| 0 avec `--hguard on` | 280,7 | 5,059 | 16,61 | 1706 | 6934 | 8640 | 13362 | 2,1 s |

Budget 6 s (10 runs) : `--deg4 0` donne 293,5 distinctes (meilleure 297, H = 5,036) ; **`--deg4 0.2` donne 307,9
(meilleure 311, H = 5,121)**, d ≤ 1 = 8615.

- Le plateau vers 292 venait bien du voisinage. Avec les doubles échanges, le dé-clonage progresse encore entre 1,5 s
  et 6 s.
- Le score lexicographique suffit à préserver l'entropie : H(X) +0,06 nat, D1 16,87 → 16,22.
- Toute probabilité entre 0,05 et 0,2 donne le même résultat. |F| = 2 n'apporte rien.
- `--hguard on` (aucun mouvement qui baisse Σ H(X_{−b})) coûte cher en distinctes : les chemins qui en gagnent passent
  souvent par une petite baisse d'entropie.

### Analyses de l'optimum

**Conclusion : les arguments locaux ne bornent rien, et les petits groupes sont rigides.**

- **`Bounds`** propage lo_x ≤ n_x ≤ hi_x cellule par cellule sur les 11 453 personas actives. Il ne conclut rien
  (lo = 0 partout, borne = 500) : chaque cellule contient trop de personas possibles.
- **`Lns`** reconstruit avec `ColumnGreedy` un sous-groupe de m individus proches d'un doublon, sur ses propres comptes.
  Aucune amélioration en 8 s depuis 302. Pour m = 30, 79 % des reconstructions redonnent exactement le même groupe
  (m = 60 : 34/58 ; m = 120 : 109/352). Avec les tables d'arité 3 et 4, les comptes d'un petit groupe le déterminent
  presque entièrement.

## Pistes restantes (en restant dans un greedy)

- Doubles échanges sur deux ensembles S différents pour les deux paires.
- Chaînes d'éjection : un écart compensé en plusieurs paires (degré 6 et plus).
- Énumération exacte des mouvements sur 3 ou 4 individus.
- Combiner le greedy 2 (meilleure population brute) avec un dé-clonage adapté à des doublons concentrés.
- Choisir l'ordre des attributs en fonction de la diversité finale, pas seulement de la faisabilité.

Déjà testé sans gain : dé-clonage après chaque colonne, cycles à trois individus, recuit, rayon 2.
