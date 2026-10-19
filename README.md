# Population synthétique : greedy colonne par colonne

Reconstruit N individus (K attributs catégoriels) qui respectent exactement des tables de marginales,
en limitant le nombre de clones (individus identiques). Tout (K, N, domaines, tables, arités) est lu depuis le CSV
`cell_id,table_id,arity,scope,cell,target`.

La méthode reste un **greedy** (construction + post-traitement local qui garde l'exactitude). CP-SAT ne sert que de
référence, pour connaître l'optimum.

## Fichiers

| Fichier | Rôle |
|---|---|
| `src/ColumnGreedy.java` | greedy 1, colonne par colonne (par défaut) |
| `src/ColumnGreedy2.java` | greedy 2 : places uniques d'abord, doublons concentrés (`--greedy 2`) |
| `src/Declone.java` | dé-clonage : échanges partiels qui gardent les marges exactes (population complète ou partielle) |
| `src/GreedyMain.java` | runs, statistiques, écriture de `population.csv` |
| `src/Evaluation.java` | vérification indépendante (seul juge des résultats) et métriques |
| `src/CheckPopulation.java` | revérifie un fichier population |
| `src/Bounds.java` | analyse : propagation de bornes sur les personas |
| `src/Lns.java` | analyse : reconstruction de sous-groupes (non branchée) |
| `src/Deg4Diag.java` | analyse : compte les doubles échanges compensés (degré 4) améliorants sur une population |
| `src/Main.java`, `GreedyBuilder.java`, `LocalSearch.java`, `Population.java` | ancienne référence : greedy individuel + recuit |
| `tools/reference_cpsat.py` | optimum exact par CP-SAT (référence) |
| `tools/reference_ipf.py` | référence de réalisme : loi d'entropie maximale P* (IPF) et ressemblances d'un échantillon typique |

## Commandes (PowerShell, JDK 21)

```powershell

javac -encoding UTF-8 -d out (Get-ChildItem src\*.java).FullName
java -cp out GreedyMain constraint_cells.csv --runs 10                         # écrit population.csv (meilleur run)
java -cp out GreedyMain constraint_cells.csv --runs 20 --quiet                 # statistiques
java -cp out GreedyMain constraint_cells.csv --runs 20 --quiet --greedy 2      # second greedy
java -cp out GreedyMain constraint_cells.csv --runs 20 --quiet --declone 0     # greedy seul, sans dé-clonage
java -cp out CheckPopulation constraint_cells.csv population.csv               # revérifie un fichier écrit
java -cp out Bounds constraint_cells.csv population.csv                        # bornes par propagation
java -cp out Deg4Diag constraint_cells.csv 1 population.csv                    # mouvements de degré 4 améliorants (|F| ≤ 1)
java -cp out GreedyMain constraint_cells.csv --runs 20 --quiet --deg4 0        # sans doubles échanges (ancien défaut)
java -cp out Main constraint_cells.csv --runs 2 --time 8000                    # ancienne référence

# optimum exact (Python 3.12 : pip install ortools pandas) ; '-' à la place de population.csv = sans solution de départ
python tools\reference_cpsat.py constraint_cells.csv 900 population.csv cpsat_population.csv

# référence de réalisme (numpy) : statistiques d'un échantillon de N individus tiré selon P*
python tools\reference_ipf.py constraint_cells.csv
```

Pour des accents corrects dans la console : `[Console]::OutputEncoding = [Text.Encoding]::UTF8` et `java -Dstdout.encoding=UTF-8 ...`.

Options de `GreedyMain` (défauts entre parenthèses) :
- construction : `--greedy` (1), `--beta` (2.0, séparation des jumeaux), `--gamma` (1.0, concentration du greedy 2),
  `--entropy off|full|loo` (full, préférence par gain d'entropie de Shannon), `--lambda` (1.0, poids des H(X_{-b})),
  `--alpha` (0.3, poids de ln(mr) en mode entropie),
  `--fc on|off` (look-ahead, on), `--retries` (3), `--backtrack` (40), `--depth` (2), `--restarts` (10), `--repair` (500),
  `--perfect` (1000), `--cross` (0) ;
- dé-clonage : `--declone` (1500 ms, 0 = désactivé), `--dentropy off|full|loo|near` (loo, critère d'entropie),
  `--dradius` (1, projections X_{−B} avec |B| = r), `--dmu` (−1 = lexicographique : distinctes puis entropie ; μ ≥ 0 :
  objectif pondéré μ × distinctes + N × entropie), `--cycle` (0, probabilité de tenter un cycle à trois individus), `--deg4` (0.2, probabilité de tenter un double
  échange compensé), `--dfmax` (1, |F| maximal), `--hguard on|off` (off, refuse tout mouvement qui baisse le critère
  d'entropie),
  `--neutral` (1.0), `--temp` (0, recuit), `--concentrate on|off` (off),
  `--small` (0) ;
- dé-clonage après chaque colonne : `--coldeclone` (0 ms par colonne = désactivé), `--colfrom` (1, première position),
  `--colneutral` (1.0), `--colspread on|off` (off), `--colbound on|off` (off, objectif Σ min(m, potentiel)) ;
- runs : `--starts` (1 construction + dé-clonage par run, on garde la meilleure), `--runs`, `--seed`, `--out`, `--quiet`.

## Algorithme

### Greedy 1 (`ColumnGreedy`)

1. **Ordre des attributs** : glouton, on privilégie les étapes « parfaites » (parmi les projections
   scope(t) ∩ connus des tables de a, l'une contient toutes les autres : colonne toujours faisable),
   puis l'attribut le plus relié aux colonnes déjà fixées.
2. **Colonne a** (pour toute la population) : chaque table contenant a impose des comptes exacts sur sa
   projection (connus + a). Une colonne exacte n'est jamais cassée ensuite.
3. Dans une colonne : individus les plus contraints d'abord. Une valeur doit avoir du résiduel dans toutes
   ses cellules et respecter le **look-ahead entier** : pour chaque attribut futur c, un groupe d'individus
   (mêmes valeurs connues sur les attributs liés à c) contient au plus cap(g) = Σ_w min_t cible_t(g, w)
   individus (borne u_x du cours sur les personas partielles). Parmi les valeurs possibles, on sépare les
   jumeaux (individus identiques sur les colonnes déjà fixées) : préférence mr / (1 + tw)^beta, avec mr le plus
   petit résiduel et tw le nombre de jumeaux déjà envoyés sur la valeur.
4. Impasse : micro-réparation de la colonne, nouveaux tirages, puis **retour arrière** sur 1 puis 2 colonnes :
   on refait une ou deux colonnes déjà fixées (attributs des tables violées), a étant provisoirement inconnue,
   puis a.
5. **Relance** : si la population reste inexacte, on relance toute la construction (au plus 10 fois).

### Greedy 2 (`ColumnGreedy2`, `--greedy 2`)

Identique au greedy 1, sauf la préférence entre valeurs faisables. Branche = persona partielle (attributs connus,
a = v). Potentiel = nombre de personas complètes actives (u_x > 0) qui la prolongent : la branche ne peut pas contenir
plus de personas distinctes.
- place libre (potentiel > jumeaux déjà envoyés) : on sépare les jumeaux comme le greedy 1 ;
- branche pleine : l'individu sera un doublon, on le met avec le plus de jumeaux (mr × (1 + tw)^gamma), pour
  concentrer les doublons comme à l'optimum.

### Dé-clonage après chaque colonne (`--coldeclone`)

Après chaque colonne exacte (sauf la dernière), on applique le même dé-clonage à la population **partielle**,
restreint aux attributs connus C : S ⊆ C, et i, j doivent coïncider sur B(S) ∩ C seulement. Pour toute table t,
scope(t) ∩ C est inclus dans S ∪ (B(S) ∩ C) ou disjoint de S : la projection de t sur les attributs connus garde ses
comptes, donc la population partielle reste exacte et les colonnes suivantes restent faisables comme avant. Les
frontières étant plus petites, il y a beaucoup plus d'échanges possibles qu'à la fin. Objectif : personas partielles
distinctes, ou (`--colbound on`) Σ_g min(m_g, pot(g)), avec pot(g) le potentiel du greedy 2 (une persona partielle
de multiplicité m_g ne donnera jamais plus de min(m_g, pot(g)) distinctes).

### Dé-clonage (`Declone`, sur une population exacte)

On échange les valeurs de deux individus i, j sur un ensemble d'attributs S, seulement si i et j coïncident sur la
frontière B(S) (attributs hors S qui partagent une table avec S). Toute table est alors incluse dans S ∪ B(S) ou
disjointe de S : les marges restent exactes à chaque mouvement. On part d'un individu en doublon, on prend le
meilleur partenaire j (gain en distinctes), on accepte les gains positifs et les gains nuls, et on garde la meilleure
population. Variantes : recuit (`--temp`), concentration des doublons à gain égal (`--concentrate`), départ dans les
petits paquets (`--small`), plusieurs constructions par run (`--starts`).

### Double échange compensé (`--deg4`, mouvements de degré 4)

**Le degré 2 est épuisé.** Si deux individus i, j deviennent i′, j′ avec les mêmes projections sur toutes les tables,
alors sur chaque table i′ prend ses valeurs entièrement chez i ou chez j : les attributs pris chez j forment un S
qu'aucune table ne relie à un attribut où i et j diffèrent, c'est exactement l'échange sur S à frontière égale. Le
voisinage de `Declone` contient donc déjà tous les mouvements à deux individus. Le graphe des tables a des cycles
qu'aucune table ne couvre (A3–A9–A11, A0–A2–A4…) : le modèle n'est pas décomposable et sa base de Markov contient des
mouvements de degré ≥ 4.

**Mouvement.** (i, j) échangent S alors qu'ils diffèrent sur F ⊆ B(S), 1 ≤ |F| ≤ `--dfmax` : seul, ce mouvement casse
les tables Tf qui touchent F et S. Avec U = S ∩ ∪Tf et C = ∪Tf \ (S ∪ F), une seconde paire (k, l) échangée sur le même
S annule exactement l'écart si
- k = (valeurs de i sur F, de j sur U, de i sur C), l = (valeurs de j sur F, de i sur U, de i sur C) ;
- k et l coïncident sur B(S) \ F.

Pour t ∈ Tf, les cellules ajoutées par (i, j) sont celles que retire (k, l), et inversement ; toute autre table est
dans S ∪ (B(S) \ F) ou disjointe de S pour chaque paire. Les quatre individus sont forcément distincts. Recherche : i
en doublon, S tiré parmi tous les ensembles (y compris ceux où S ∪ B(S) couvre tout), jusqu'à 4 partenaires j
cassants, puis toutes les paires (k, l) trouvées en un passage sur la population. Score : le même que l'échange
(distinctes, puis Σ_b H(X_{−b})), évalué sur les 4 personas.

**Diagnostic** (`Deg4Diag`, sur des populations finales du défaut, |F| = 1) : toutes les paires de personas et tous
les S, écart haché linéairement (Zobrist additif, δ₁ + δ₂ = 0 en O(1)), chaque mouvement revérifié table par table.

| Population | Échanges valides | dont gain de distinctes | Mouvements de degré 4 exacts | dont gain de distinctes | dont entropie en hausse à distinctes égales |
|---|---|---|---|---|---|
| seed 0 (293 distinctes) | 17 835 | **0** | 309 443 | **137** (max +2) | 1870 |
| seed 1 (296 distinctes) | 19 137 | **0** | 371 707 | **211** | 2403 |

Aucun faux positif du hachage. La population est un optimum local pour les échanges, mais pas pour le degré 4.

### Entropie de Shannon dès la construction (`--entropy`)

Notations : pour un ensemble d'attributs A, H(X_A) = ln N − (1/N) Σ_c n_c ln n_c (nats), sur la loi empirique de la
population. Ajouter un individu à une cellule de n individus augmente Σ n ln n de δ(n) = (n+1) ln(n+1) − n ln n.
- Si A est inclus dans le scope d'une table, H(X_A) est imposée par les marges ; seules les projections non couvertes
  sont optimisables. Si P̂ respecte exactement les tables et P* est la loi d'entropie maximale sous ces tables (IPF),
  KL(P̂ ‖ P*) = H(P*) − H(P̂) : maximiser H(P̂) rapproche de P*.
- Règle de chaîne : H(X) = Σ_k H(X_{a_k} | X_{a_1..a_{k−1}}) ; la construction colonne par colonne ajoute un terme par
  colonne.

Seule la préférence entre valeurs faisables change (faisabilité, look-ahead, réparation, retours arrière inchangés : la
population reste exacte). Score de v pour l'individu i : 5·10⁴ − gain + α ln(mr), avec
- `full` : gain = δ(n_{g,v}), n_{g,v} = jumeaux (même signature sur les connus) déjà envoyés sur v. Chaque colonne
  maximise H(X_a | X_C) de façon gloutonne ;
- `loo` : gain = δ(n_{g,v}) + λ Σ_{b connu} δ(n_{−b,v}), n_{−b,v} = compte de la cellule de X_{(C ∪ a) \ b}. Objectif
  Σ_b H(X_{−b}) = K·H(X) − Σ_b H(X_b | X_{−b}) ; H(X_b | X_{−b}) = 0 ssi aucune paire d'individus ne diffère seulement
  sur b : pénalise les quasi-clones (distance 1). Clés masquées précalculées : O(K) par (i, v).

Dans le dé-clonage final (`--dentropy`), même idée en critère secondaire, à distinctes égales (comme `--concentrate`) :
score = Δdistinctes × 2³² + 10⁶ × gain, avec gain = −Δ Σ n ln n sur les cellules de X (`full`, maximise H(X)), de
tous les X_{−b} (`loo`, maximise Σ_b H(X_{−b})), ou K·Δ_X − Δ_{loo} (`near`, minimise Σ_b H(X_b | X_{−b})). Comptes
tenus à jour incrémentalement (les clés masquées se déduisent du code de la persona) ; la meilleure population est
choisie selon (distinctes, critère secondaire). Les mouvements restent les mêmes : exactitude garantie.

Mesures ajoutées à `Evaluation` : déficit D = K·ln N − Σ_b H(X_{−b}), Σ_b H(X_b | X_{−b}), paires à distance 0
(clones) et 1, Σ I(X_a ; X_b) sur les 19 paires d'attributs couvertes par aucune table.

**Objectif d'abord envisagé : réduire les ressemblances (clones et quasi-clones).** Abandonné ensuite : les
ressemblances du défaut sont déjà au niveau d'une population réaliste (voir « Référence de réalisme »). H(X_{−b}) = ln N ssi tous les individus restent
distincts quand on ignore b ; donc D = 0 ssi aucune paire d'individus n'est à distance ≤ 1. C'est le critère `loo`.
Σ_b H(X_b | X_{−b}) (`near`) ne compte que la distance exactement 1 : il se réduit en transformant des quasi-clones
en clones, et n'est pas le bon objectif. Défaut : `--entropy full --dentropy loo` (plus `--deg4 0.2`, voir plus bas).

`Evaluation.evaluate` recalcule tout depuis zéro : c'est le seul juge des résultats affichés.

## Résultats (constraint_cells.csv, N = 500, K = 12)

Run exact : population de taille N, lignes valides, erreur L1 = 0 sur toutes les cellules. Tous les runs des
versions actuelles sont exacts.

### Synthèse

| | Distinctes moy. | Meilleure | Temps / run |
|---|---|---|---|
| greedy 1 seul (beta = 2) | 251,9 | 266 | 0,7 s |
| greedy 2 seul | 258,8 | 270 | 1,1 s |
| greedy 1 + dé-clonage 1,5 s (ancien défaut) | 292,5 | 297 | 2,6 s |
| + `--entropy full --dentropy loo` | 292,8 | 299 | 2,2 s |
| **+ `--deg4 0.2` (défaut, doubles échanges compensés)** | **302,5** | 307 | 2,2 s |
| défaut, dé-clonage 6 s | 307,9 | 311 | 6,6 s |
| greedy 2 + dé-clonage 1,5 s | 290,2 | 296 | 2,9 s |
| greedy 1 + dé-clonage, 8 constructions par run | 297,0 | 300 | 15 s |
| **optimum (CP-SAT, prouvé)** | **331** | 331 | 900 s |

L'entropie dans le dé-clonage garde les distinctes et réduit les ressemblances : déficit K·ln N − Σ_b H(X_{−b})
17,26 → 16,86, paires à distance ≤ 1 : 9367 → 9050 (voir « Entropie dans le dé-clonage final »). Les doubles échanges
compensés ajoutent 10 distinctes et +0,06 nat (voir « Double échange compensé : réglages testés »).

Le greedy atteint environ 91 % de l'optimum en moyenne (302,5 / 331 ; meilleur run : 311, soit 94 %).

### Historique du greedy 1

| Version | Exacts | Temps / run | Distinctes (exacts) |
|---|---|---|---|
| ColumnGreedy initial | 2/20 | 0,1–1,8 s | ~240 |
| + retour arrière ciblé | 62,5 % | 0,8 s | 244 |
| + look-ahead par capacités | 85,5 % | 1,5 s | 244 |
| + ordre « étapes parfaites », 40 retours arrière | 293/300 (97,7 %) | 0,76 s (max 5,9 s) | 246 |
| + retour arrière sur 2 colonnes, sans relance | 200/200 | 0,47 s (max 5,7 s) | 246,5 |
| + relances | 600/600, dont 596 du premier coup | 0,70 s (max 8,4 s) | 246,3 |
| beta = 2 | 40/40 | 1,1 s | 251,1 |
| + dé-clonage 1,5 s | 40/40 | 2,3 s (max 7,0 s) | 293,0 (max 302) |

### Dé-clonage : réglages testés

- **Budget** : il plafonne en environ 1 s (beta = 1 : 284,5 en 1 s, 285,3 en 8 s). C'est un optimum local des
  échanges à deux.
- **Mouvements neutres** : indispensables (sans eux : 277 au lieu de 285).
- **Recuit** (6 runs, 4 s) : T₀ = 0,3 / 0,6 / 1 / 2 donne 293,0 / 293,5 / 293,7 / 293,0. Aucun gain.
- **Multi-départ** (6 runs, à budget comparable) : 1 × 8 s → 293,2 ; 4 × 2 s → 295,3 ; 8 × 1 s → 297,0 (15 s par run).
  Gain de +2 à +4 pour un temps 4 à 8 fois plus long.
- **Viser la structure de l'optimum** (8 runs) :

| Variante | Distinctes moy. | Meilleure | Personas en doublon (meilleure) |
|---|---|---|---|
| aucune (défaut) | 293,0 | 296 | 35 |
| `--concentrate on` | 293,4 | 297 | 27 |
| `--small 2` | 293,9 | 297 | 23 |
| `--concentrate on --small 2` | 293,8 | 297 | 20 |

La structure se rapproche de l'optimum (20 personas en doublon au lieu de 35 ; l'optimum en a 10), mais le nombre
de distinctes ne bouge pas : il manque des exemplaires uniques (277 contre 321 à l'optimum).

**Limite du voisinage** : seuls 85 ensembles S sur 2047 laissent un attribut hors de S ∪ B(S), car le graphe des
tables est dense. C'est ce qui bloque le dé-clonage, plus que la stratégie d'acceptation.

### Greedy 2 : réglages testés (20 runs, sans dé-clonage)

| gamma | Distinctes moy. | Meilleure |
|---|---|---|
| 0,5 | 256,7 | 268 |
| **1** | **258,8** | 270 |
| 2 | 255,6 | 262 |

Une première version, qui préférait toujours les branches à fort potentiel sans séparer les jumeaux, faisait moins bien
(243,7) : au début, toutes les branches ont un potentiel de plusieurs centaines. Le greedy 2 bat le greedy 1 sur la
construction seule (+6,9), mais pas après dé-clonage (−2,3). Hypothèse non vérifiée : les doublons concentrés laissent
moins d'échanges utiles au dé-clonage.

### Dé-clonage après chaque colonne : réglages testés

Tous les runs sont exacts. Budget par colonne, 10 colonnes dé-clonées (positions 1 à 10).

| Variante | Runs | Distinctes moy. | Meilleure | Temps / run |
|---|---|---|---|---|
| greedy 1 seul | 20 | 251,9 | | 0,9 s |
| greedy 1 + `--coldeclone 100` | 20 | **274,5** | | 1,3 s |
| greedy 1 + `--coldeclone 100 --colspread on` | 20 | 274,0 | | 1,5 s |
| greedy 1 + `--coldeclone 100 --colbound on` | 20 | 258,3 | 268 | 1,7 s |
| greedy 2 + `--coldeclone 100` | 20 | **277,1** | 289 | 2,3 s |
| greedy 2 + `--coldeclone 100 --colbound on` | 20 | 264,4 | 277 | 2,1 s |
| **avec dé-clonage final 1,5 s** | | | | |
| défaut (seeds 100–139) | 40 | 292,2 | 302 | 2,5 s |
| `--coldeclone 100` (seeds 100–139) | 40 | 292,2 | 300 | 3,2 s |
| `--coldeclone 100 --colbound on` (seeds 100–139) | 40 | 292,5 | 299 | 3,8 s |
| `--coldeclone 30 / 100 / 300` (seeds 0–19) | 20 | 290,1 / 292,8 / 291,9 | | 2,8 / 3,1 / 4,4 s |
| `--coldeclone 100 --colfrom 5 / 8` | 20 | 292,8 / 291,9 | 299 / 297 | 3,0 s |
| `--coldeclone 100 --colneutral 0.3` | 20 | 291,0 | | 3,0 s |
| `--coldeclone 100 --colspread on` | 20 | 293,0 | | 3,1 s |
| greedy 2 + `--coldeclone 100` (± `--colbound on`) | 20 | 288,9 / 292,4 | 294 / 299 | 3,5 s |

- **Sur la construction**, l'effet est fort : +22,6 distinctes pour le greedy 1 (251,9 → 274,5) et 277,1 pour le
  greedy 2, soit la meilleure population brute obtenue jusqu'ici.
- **Après le dé-clonage final**, le gain disparaît : 292,2 dans les deux cas sur 40 runs. Le dé-clonage final
  atteint le même plateau (~292) quelle que soit la population de départ (251 ou 275 distinctes).
- L'objectif borné `--colbound` donne une construction moins bonne (258,3) et le même plateau final : forcer
  m_g ≤ pot(g) tôt ne crée pas plus d'exemplaires uniques à la fin.
- Conclusion : les copies évitables restantes (~39 par rapport à l'optimum) ne viennent pas de la construction,
  mais de la rigidité du voisinage d'échanges complet (85 ensembles S utiles sur 2047). C'est désactivé par défaut.

### Entropie à la construction : réglages testés

Sans dé-clonage (`--declone 0`), 20 runs (seeds 0–19), tous exacts. H en nats (ln 500 = 6,215).

| Variante | Distinctes | H(X) | Σ H(X_b \| X_−b) | Paires à distance 1 | Σ I paires libres |
|---|---|---|---|---|---|
| off (beta = 2, défaut) | 251,9 | 4,943 | 2,478 | 7131 | 0,372 |
| full, α = 0 | 240,4 | 4,863 | 2,320 | 7245 | 0,375 |
| full, α = 0,1 | 248,3 | 4,954 | 2,497 | 6929 | 0,353 |
| **full, α = 0,3** | **256,9** | **4,959** | 2,505 | 7056 | 0,357 |
| full, α = 0,5 / 1 / 3 | 258,1 / 256,7 / 245,6 | 4,948 / 4,924 / 4,867 | | | |
| loo, α = 0, λ = 1 | 240,1 | 4,865 | 2,272 | 7167 | 0,355 |
| loo, α = 0,3, λ = 0,03 / 0,1 / 1 | 251,9 / 247,8 / 242,0 | 4,950 / 4,929 / 4,881 | 2,484 / 2,423 / 2,285 | ~7045 | 0,35–0,36 |
| loo, α = 1, λ = 0,1 / 0,3 / 1 / 3 | 257,1 / 253,3 / 244,6 / 242,0 | 4,935 / 4,938 / 4,915 / 4,872 | | | |

Avec dé-clonage final 1,5 s (20 runs) : off 292,5 distinctes, H = 4,994 ; full α = 0,3 : 294,0, H = 5,002 ;
loo α = 0,3 λ = 0,1 : 290,8, H = 4,992.

Contrôle sur les seeds 100–139 (40 runs, sans dé-clonage) : off H = 4,947 ± 0,029 (écart-type) ; full α = 0,3 :
4,950 ± 0,024 (+4,7 distinctes) ; full α = 0,1 : 4,955 ± 0,023.

- **Le gain glouton de Shannon seul (α = 0) est moins bon que l'heuristique actuelle** (4,863 contre 4,943) : il est
  myope (il étale chaque groupe sans regarder la marge des cellules, ce qui force des paquets aux colonnes suivantes).
  Le terme α ln(mr) est indispensable ; le meilleur réglage est α ≈ 0,1–0,3.
- Le gain sur H(X) reste faible (+0,01 nat, à la limite du bruit sur 40 runs) ; full α = 0,3 donne en revanche
  +5 distinctes de façon reproductible, et le léger avantage se conserve après dé-clonage.
- `loo` diminue bien Σ H(X_b | X_−b), mais surtout parce que H(X) baisse : l'objectif propre Σ_b H(X_−b) = K·H − Σ H(b|−b)
  est meilleur avec `full` (57,0) qu'avec `loo` (56,6–56,8) ou off (56,8). Optimiser λ Σ δ colonne par colonne ne suffit
  pas à le maximiser.
- Les paires à distance 1 (~7000) sont dominées par les gros paquets imposés (un paquet ×48 voisin d'un ×49 donne à
  lui seul 2352 paires) : cette mesure bouge peu.

### Entropie dans le dé-clonage final : réglages testés

20 runs (seeds 0–19), dé-clonage 1,5 s, tous exacts.

| Construction | `--dentropy` | Distinctes | H(X) | Σ H(X_b \| X_−b) | Paires à distance 1 | Σ I paires libres |
|---|---|---|---|---|---|---|
| off | off (défaut) | 292,5 | 4,993 | 2,607 | 7179 | 0,392 |
| off | full | 292,4 | 5,024 | 2,683 | 7156 | 0,380 |
| off | loo | 291,4 | 5,023 | 2,667 | 7136 | 0,374 |
| off | near | 291,4 | 4,944 | **2,424** | 6983 | 0,388 |
| **full (α = 0,3)** | **full** | **293,4** | **5,036** | 2,716 | 7100 | **0,370** |
| loo (α = 0,3, λ = 0,1) | near | 290,5 | 4,941 | 2,426 | **6903** | 0,383 |

- **`--dentropy full` augmente H(X) de +0,03 nat sans perdre de distinctes** (écart-type par run ≈ 0,025 : gain net
  sur 20 runs). Combiné à la construction `--entropy full` : **H = 5,036 (+0,043) et 293,4 distinctes**, le meilleur
  couple obtenu.
- **H(X) et les quasi-clones s'opposent à distinctes égales.** Monter H, c'est égaliser les multiplicités : on étale les
  doublons en personas voisines, et Σ H(X_b | X_−b) monte (2,61 → 2,68). `near` fait l'inverse : quasi-clones −7 %
  (2,61 → 2,42) mais H baisse (4,94), car il concentre les doublons sur peu de personas.
- `loo` (Σ H(X_−b)) se comporte comme `full` : à distinctes fixées, il est dominé par le terme K·H(X).

En comptant aussi les clones (paires à distance 0), nouvelle série de 20 runs :

| Construction / `--dentropy` | Distinctes | Paires d = 0 | Paires d = 1 | Total d ≤ 1 |
|---|---|---|---|---|
| off / off (ancien défaut) | 292,7 | 2200 | 7167 | 9367 |
| off / near | 291,9 | 2408 | 6982 | 9390 |
| off / full | 292,4 | 2050 | 7166 | 9215 |
| off / loo | 291,4 | 2037 | 7121 | 9157 |
| loo (λ = 0,1) / loo | 290,1 | 2005 | 7103 | 9108 |
| **full / loo (défaut)** | **292,5** | **1975** | **7075** | **9050** |

`near` baisse les paires à distance 1 (−185) mais ajoute autant de clones (+208) : aucune ressemblance en moins.
`full / loo` retire 317 paires à distance ≤ 1 (−3,4 %) sans perdre de distinctes. Les paires de clones viennent
surtout des deux gros paquets imposés (×46 à ×49, plus de 1000 paires chacun).

### Objectif pondéré, rayon 2, cycles à trois individus

`Declone` évalue maintenant un mouvement modifiant m individus (m = 2 ou 3) ; les comptes de Σ n ln n sont tenus
dans `dec`/`inc` et, pour `loo`, sur toutes les projections X_{−B}, |B| = `--dradius`.
- **Objectif pondéré** (`--dmu μ`) : F = μ × distinctes + N × Σ_B H(X_{−B}). Les clones sont déjà pénalisés par
  Σ_B H(X_{−B}) ; avec μ petit, on accepte de perdre des distinctes pour retirer des ressemblances.
- **Rayon 2** (`--dradius 2`) : Σ_{|B|=2} H(X_{−B}) sur les 66 paires d'attributs ; maximal ssi aucune paire
  d'individus à distance ≤ 2. Mesuré par le déficit D2 = C(K,2)·ln N − Σ_{|B|=2} H(X_{−B}) et les paires à distance 2.
- **Cycles** (`--cycle p`) : i ← j ← k ← i sur S, les trois individus égaux sur B(S) (donc exact). Un cycle est la
  composée de deux échanges valides : il ne rend rien accessible de nouveau, il franchit seulement un état
  intermédiaire moins bon.

20 runs (seeds 0–19), dé-clonage 1,5 s, construction `--entropy full`, tous exacts :

| `Declone` | Distinctes | D1 | Paires d = 0 | d = 1 | **d ≤ 1** | d = 2 | D2 |
|---|---|---|---|---|---|---|---|
| lexicographique, r = 1 (défaut) | **292,8** | 16,87 | 1980 | 7066 | 9046 | 13073 | 109,37 |
| `--dmu 30` | 287,1 | 16,76 | 1881 | 6976 | 8856 | 13193 | 108,89 |
| `--dmu 10` | 285,2 | 16,67 | 1796 | 6956 | 8753 | 13293 | 108,51 |
| `--dmu 3` | 280,2 | 16,60 | 1693 | 6932 | 8624 | 13357 | 108,17 |
| `--dmu 0` | 276,9 | **16,57** | **1631** | 6893 | **8524** | 13396 | **108,03** |
| `--dradius 2` | 290,2 | 16,90 | 1971 | 7060 | 9031 | 13066 | 109,50 |
| `--dradius 2 --dmu 0` | 276,1 | 16,60 | 1634 | **6846** | **8480** | 13419 | 108,05 |
| `--cycle 0.3` / `0.6` | 290,5 / 287,2 | 16,87 / 16,89 | | | 9005 / 8961 | | 109,40 / 109,47 |

Budget 6 s (10 runs) : défaut 293,5 distinctes, d ≤ 1 = 9051 ; `--dradius 2` 293,0 / 9028 ; `--cycle 0.3`
292,8 / 9023 ; `--dmu 0` 276,7 / 8494 ; `--dmu 0 --dradius 2` 275,8 / 8431.

- **L'objectif pondéré est le seul levier efficace**, et il donne un front : de μ = −1 (lexicographique) à μ = 0,
  les paires à distance ≤ 1 baissent de 6 % (9046 → 8524), les paires de clones de 18 % (1980 → 1631), au prix de
  16 distinctes (292,8 → 276,9). Les doublons sont répartis autrement (5 runs) : multiplicité max 36 en moyenne au lieu de 39
  (moins de paires de clones), mais 275 individus en doublon au lieu de 242.
- **Rayon 2 : aucun gain sur la distance 2.** La somme des distances de Hamming sur toutes les paires est fixée par
  les marges à une variable : retirer des paires à distance ≤ 1 en crée ailleurs, et elles arrivent surtout à
  distance 2 (13073 → 13396 avec μ = 0). Dans Σ_{|B|=2} H(X_{−B}), une paire de clones compte dans les 66 projections,
  une paire à distance 1 dans 11, une paire à distance 2 dans une seule : l'objectif reste dominé par d ≤ 1 et ne
  fait que pousser des paires de 2 vers 3. Avec μ = 0 il retire un peu plus de paires à distance 1 (−1 %).
- **Cycles : aucun gain**, à budget égal (1,5 s ou 6 s). Ils coûtent plus cher par itération et, étant des composées
  d'échanges sur le même S, n'élargissent pas le voisinage. La limite reste la condition de frontière B(S)
  (85 ensembles S utiles sur 2047).

Défaut inchangé (lexicographique, r = 1, sans cycle) : il est au niveau de ressemblance réaliste (voir « Référence de
réalisme »). `--dmu 0` (ou `--dmu 0 --dradius 2`) ne sert que si l'on vise moins de ressemblances que le réalisme.

### Référence de réalisme (`tools/reference_ipf.py`)

**Question** : faut-il réduire les quasi-clones (individus qui ne diffèrent que sur un attribut) pour être plus
réaliste ? **Réponse : non, le défaut est déjà au niveau réaliste.**

- **Les marges fixent la distance moyenne.** Pour chaque attribut, la proportion de paires qui diffèrent ne dépend que de
  sa marge à une variable ; la distance de Hamming moyenne entre deux individus vaut donc **4,57 sur 12** pour toute
  population exacte. Quatre attributs sont presque constants (A8 : 98 % sur une valeur, A0 : 92 %, A6 : 91 %,
  A10 : 90 %). Retirer des paires proches en crée forcément plus loin.
- **La population réelle contient des clones.** Si les marges viennent d'une population réelle de 500 personnes, celle-ci
  est une population exacte : d'après CP-SAT (au plus 331 distinctes), elle compte au moins 169 copies en trop.
- **Échantillon typique.** P* = loi d'entropie maximale sous les 45 tables (IPF sur les 3¹² personas), c'est-à-dire la
  loi la moins structurée compatible avec les marges. Espérances pour N = 500 individus tirés selon P* :

| | Distinctes | Paires d = 0 | d = 1 | d ≤ 1 | d = 2 | Voisins à d ≤ 1 par individu |
|---|---|---|---|---|---|---|
| **échantillon typique de P*** | 232 | 1935 | 7151 | **9086** | 13420 | 36,3 |
| défaut (`--entropy full --dentropy loo`) | 292,8 | 1980 | 7066 | **9046** | 13073 | 36,2 |
| ancien défaut (`--entropy off --dentropy off`) | 292,7 | 2200 | 7167 | 9367 | | 37,5 |
| `--dmu 0` | 276,9 | 1631 | 6893 | 8524 | 13396 | 34,1 |

- Le défaut a **le même niveau de ressemblance qu'un échantillon typique** (9046 contre 9086 paires à distance ≤ 1),
  avec plus de distinctes (293 contre 232). Les quasi-clones viennent des données, pas de l'algorithme ; une vraie
  population, plus corrélée que P*, en aurait en général plutôt davantage.
- `--dmu 0` descend 6 % sous ce niveau, au prix de 16 distinctes : **moins réaliste**. À n'utiliser que si l'on vise
  autre chose que le réalisme (confidentialité, diversité imposée).
- **Choix retenu : le défaut.** Les ressemblances se lisent par rapport à la référence P* : on les calibre, on ne les
  minimise pas.

### Optimum et clones imposés

- **CP-SAT** (`tools/reference_cpsat.py`, 8 threads, population greedy à 302 en solution de départ) : 330 (borne 331)
  en 300 s, puis **OPTIMAL = 331 distinctes** en 900 s. La population obtenue est vérifiée exacte par `CheckPopulation`.
  Donc **169 copies en trop sont imposées par les marges**, et il reste au plus 29 copies évitables par rapport au
  meilleur greedy.

| | Distinctes | Copies en trop | Personas en doublon | Individus dans un doublon | Multiplicité max |
|---|---|---|---|---|---|
| greedy 1 + dé-clonage (meilleur run) | 302 | 198 | 29 | 227 | 48 |
| optimum CP-SAT | **331** | **169** | 10 | 179 | 49 |

- **Les gros paquets ne sont pas le problème** : l'optimum garde ×49 et ×46 sur les mêmes personas que le greedy. Il
  concentre tous ses doublons sur 10 personas et met les 321 autres individus en exemplaire unique. Le greedy répartit
  les siens sur 29 personas. Les copies évitables sont dans les petits paquets (×2 à ×7).
- **`Bounds`** (propagation lo_x ≤ n_x ≤ hi_x cellule par cellule sur les 11 453 personas actives) : ne conclut rien
  (lo = 0 partout, borne = 500). Chaque cellule contient trop de personas possibles pour cet argument local.
- **`Lns`** (on reconstruit avec `ColumnGreedy` un sous-groupe de m individus proches d'un doublon, sur ses propres comptes) :
  0 amélioration en 8 s depuis 302. Pour m = 30, 79 % des reconstructions redonnent exactement le même groupe (m = 60 :
  34/58 ; m = 120 : 109/352). Avec les tables d'arité 3 et 4, les comptes d'un petit groupe le déterminent presque
  entièrement : la population est **localement rigide** pour ces voisinages.

### Double échange compensé : réglages testés

20 runs (seeds 0–19), défaut (`--entropy full --dentropy loo`), dé-clonage 1,5 s, tous exacts.

| `--deg4` | Distinctes | H(X) | D1 | Paires d = 0 | d = 1 | **d ≤ 1** | d = 2 | Temps / run |
|---|---|---|---|---|---|---|---|---|
| 0 (ancien défaut) | 292,6 | 5,033 | 16,87 | 1988 | 7080 | 9068 | 13063 | 2,4 s |
| 0,05 | 302,2 | 5,094 | 16,26 | 1851 | 6876 | 8727 | 13294 | 2,2 s |
| 0,1 | 302,1 | 5,092 | 16,28 | 1860 | 6849 | 8709 | 13313 | 2,2 s |
| **0,2 (défaut)** | **302,5** | 5,097 | 16,22 | 1849 | 6866 | 8714 | 13271 | 2,2 s |
| 0,3 | 300,9 | 5,096 | 16,23 | 1821 | 6848 | 8669 | 13300 | 2,2 s |
| 0,5 | 301,5 | **5,105** | **16,15** | 1800 | 6787 | 8587 | 13360 | 2,3 s |
| 0,1 / 0,2 avec `--dfmax 2` | 301,7 / 302,6 | 5,087 / 5,098 | 16,33 / 16,23 | | | 8804 / 8739 | | 2,2 s |
| 0,1 avec `--hguard on` | 293,3 | 5,129 | 15,93 | 1570 | 6745 | 8315 | 13530 | 2,2 s |
| 0 avec `--hguard on` | 280,7 | 5,059 | 16,61 | 1706 | 6934 | 8640 | 13362 | 2,1 s |

Budget 6 s (10 runs) : `--deg4 0` 293,5 distinctes (meilleure 297), H = 5,036 ; **`--deg4 0.2` 307,9 (meilleure
311), H = 5,121**, d ≤ 1 = 8615.

- **+10 distinctes à budget égal** (292,6 → 302,5), **+14 à 6 s** (293,5 → 307,9) : le plateau ~292 était bien dû au
  voisinage. Le meilleur run passe de 302 à 311 (optimum 331 : 94 %).
- **L'entropie progresse aussi** : H(X) +0,06 nat, déficit D1 16,87 → 16,22. Rien n'est perdu sur le critère secondaire :
  le score lexicographique suffit, et les gains de distinctes égalisent les multiplicités.
- Dès p = 0,05 l'effet est là ; p entre 0,05 et 0,2 est équivalent. |F| = 2 n'apporte rien.
- **`--hguard on` (aucun mouvement qui baisse Σ H(X_{−b}))** coûte cher en distinctes (293 avec degré 4, 281 sans) : les
  chemins qui gagnent des distinctes passent souvent par une petite baisse d'entropie. Il maximise H (5,129) mais n'est
  pas retenu.
- **Réalisme** : les paires à distance ≤ 1 passent sous la référence P* (8714 contre 9086, −4 %). C'est une conséquence
  des distinctes en plus, pas d'un objectif de ressemblance. À surveiller au vu de « Référence de réalisme ».
- **Nouveau défaut : `--entropy full --dentropy loo --deg4 0.2`** (`--deg4 0` pour revenir à l'ancien).

## Pistes restantes (en restant dans un greedy)

- ~~Dé-cloner après chaque colonne~~ : testé (`--coldeclone`), +22 distinctes sur la construction mais aucun gain
  après dé-clonage final (voir plus haut). Le plateau ~292 vient du voisinage final, pas du point de départ.
- ~~Cycles à trois individus~~ : testé (`--cycle`), aucun gain (composées d'échanges sur le même S).
- ~~Élargir le voisinage par des mouvements de degré 4~~ : fait (`--deg4`), +10 à +14 distinctes. Suite possible : deux
  ensembles S différents pour les deux paires, chaînes d'éjection (degré 6+), énumération exacte sur 3–4 individus.
- Combiner le greedy 2 (meilleure population brute) avec un dé-clonage adapté à des doublons concentrés.
- Mieux choisir l'ordre des attributs en fonction de la diversité finale, pas seulement de la faisabilité.
