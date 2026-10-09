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
| `src/Declone.java` | dé-clonage : échanges partiels qui gardent les marges exactes |
| `src/GreedyMain.java` | runs, statistiques, écriture de `population.csv` |
| `src/Evaluation.java` | vérification indépendante (seul juge des résultats) et métriques |
| `src/CheckPopulation.java` | revérifie un fichier population |
| `src/Bounds.java` | analyse : propagation de bornes sur les personas |
| `src/Lns.java` | analyse : reconstruction de sous-groupes (non branchée) |
| `src/Main.java`, `GreedyBuilder.java`, `LocalSearch.java`, `Population.java` | ancienne référence : greedy individuel + recuit |
| `tools/reference_cpsat.py` | optimum exact par CP-SAT (référence) |

## Commandes (PowerShell, JDK 21)

```powershell

javac -encoding UTF-8 -d out (Get-ChildItem src\*.java).FullName
java -cp out GreedyMain constraint_cells.csv --runs 10                         # écrit population.csv (meilleur run)
java -cp out GreedyMain constraint_cells.csv --runs 20 --quiet                 # statistiques
java -cp out GreedyMain constraint_cells.csv --runs 20 --quiet --greedy 2      # second greedy
java -cp out GreedyMain constraint_cells.csv --runs 20 --quiet --declone 0     # greedy seul, sans dé-clonage
java -cp out CheckPopulation constraint_cells.csv population.csv               # revérifie un fichier écrit
java -cp out Bounds constraint_cells.csv population.csv                        # bornes par propagation
java -cp out Main constraint_cells.csv --runs 2 --time 8000                    # ancienne référence

# optimum exact (Python 3.12 : pip install ortools pandas) ; '-' à la place de population.csv = sans solution de départ
python tools\reference_cpsat.py constraint_cells.csv 900 population.csv cpsat_population.csv
```

Pour des accents corrects dans la console : `[Console]::OutputEncoding = [Text.Encoding]::UTF8` et `java -Dstdout.encoding=UTF-8 ...`.

Options de `GreedyMain` (défauts entre parenthèses) :
- construction : `--greedy` (1), `--beta` (2.0, séparation des jumeaux), `--gamma` (1.0, concentration du greedy 2),
  `--fc on|off` (look-ahead, on), `--retries` (3), `--backtrack` (40), `--depth` (2), `--restarts` (10), `--repair` (500),
  `--perfect` (1000), `--cross` (0) ;
- dé-clonage : `--declone` (1500 ms, 0 = désactivé), `--neutral` (1.0), `--temp` (0, recuit), `--concentrate on|off` (off),
  `--small` (0) ;
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

### Dé-clonage (`Declone`, sur une population exacte)

On échange les valeurs de deux individus i, j sur un ensemble d'attributs S, seulement si i et j coïncident sur la
frontière B(S) (attributs hors S qui partagent une table avec S). Toute table est alors incluse dans S ∪ B(S) ou
disjointe de S : les marges restent exactes à chaque mouvement. On part d'un individu en doublon, on prend le
meilleur partenaire j (gain en distinctes), on accepte les gains positifs et les gains nuls, et on garde la meilleure
population. Variantes : recuit (`--temp`), concentration des doublons à gain égal (`--concentrate`), départ dans les
petits paquets (`--small`), plusieurs constructions par run (`--starts`).

`Evaluation.evaluate` recalcule tout depuis zéro : c'est le seul juge des résultats affichés.

## Résultats (constraint_cells.csv, N = 500, K = 12)

Run exact : population de taille N, lignes valides, erreur L1 = 0 sur toutes les cellules. Tous les runs des
versions actuelles sont exacts.

### Synthèse

| | Distinctes moy. | Meilleure | Temps / run |
|---|---|---|---|
| greedy 1 seul (beta = 2) | 251,9 | 266 | 0,7 s |
| greedy 2 seul | 258,8 | 270 | 1,1 s |
| **greedy 1 + dé-clonage 1,5 s (défaut)** | **292,5** | 297 | 2,6 s |
| greedy 2 + dé-clonage 1,5 s | 290,2 | 296 | 2,9 s |
| greedy 1 + dé-clonage, 8 constructions par run | 297,0 | 300 | 15 s |
| **optimum (CP-SAT, prouvé)** | **331** | 331 | 900 s |

Le greedy atteint environ 88 % de l'optimum en moyenne (292,5 / 331 ; meilleur run : 302, soit 91 %).

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

