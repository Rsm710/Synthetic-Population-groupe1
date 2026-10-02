# Population synthétique : greedy colonne par colonne

Reconstruit N individus (K attributs catégoriels) qui respectent exactement des tables de marginales.
Tout (K, N, domaines, tables, arités) est lu depuis le CSV `cell_id,table_id,arity,scope,cell,target`.

## Commandes (PowerShell, JDK 21)

```powershell
# si java n'est pas dans le PATH :
$env:Path = "C:\Program Files\Eclipse Adoptium\jdk-21.0.12.101-hotspot\bin;" + $env:Path

javac -encoding UTF-8 -d out (Get-ChildItem src\*.java).FullName
java -cp out GreedyMain constraint_cells.csv --runs 10            # écrit population.csv (meilleur run)
java -cp out GreedyMain constraint_cells.csv --runs 300 --quiet   # taux de runs exacts
java -cp out CheckPopulation constraint_cells.csv population.csv  # revérifie un fichier écrit
java -cp out Main constraint_cells.csv --runs 2 --time 8000       # référence : greedy individuel + recuit
```

Options de `GreedyMain` : `--runs`, `--seed`, `--fc on|off` (look-ahead), `--retries` (3), `--backtrack` (40), `--depth` (2), `--restarts` (10),
`--repair` (500), `--perfect` (1000), `--cross` (0), `--beta` (1.0), `--out`, `--quiet`.
Pour des accents corrects dans la console : `[Console]::OutputEncoding = [Text.Encoding]::UTF8` et `java -Dstdout.encoding=UTF-8 ...`.

## Algorithme (`ColumnGreedy`)

1. **Ordre des attributs** : glouton, on privilégie les étapes « parfaites » (parmi les projections
   scope(t) ∩ connus des tables de a, l'une contient toutes les autres : colonne toujours faisable),
   puis l'attribut le plus relié aux colonnes déjà fixées.
2. **Colonne a** (pour toute la population) : chaque table contenant a impose des comptes exacts sur sa
   projection (connus + a). Une colonne exacte n'est jamais cassée ensuite.
3. Dans une colonne : individus les plus contraints d'abord. Une valeur doit avoir du résiduel dans toutes
   ses cellules et respecter le **look-ahead entier** : pour chaque attribut futur c, un groupe d'individus
   (mêmes valeurs connues sur les attributs liés à c) contient au plus cap(g) = Σ_w min_t cible_t(g, w)
   individus (borne u_x du cours sur les personas partielles).
4. Impasse : micro-réparation de la colonne, nouveaux tirages, puis **retour arrière** sur 1 puis 2 colonnes :
   on refait une ou deux colonnes déjà fixées (attributs des tables violées), a étant provisoirement inconnue,
   puis a.
5. **Relance** : si la population reste inexacte, on relance toute la construction (au plus 10 fois).

`Evaluation.evaluate` recalcule tout depuis zéro : c'est le seul juge des résultats affichés.

## Résultats (constraint_cells.csv, N = 500, K = 12)

Run exact : population de taille N, lignes valides, erreur L1 = 0 sur toutes les cellules (vérifié par `Evaluation.evaluate`).

| Version | Exacts | Temps / run | Distinctes (exacts) |
|---|---|---|---|
| ColumnGreedy initial | 2/20 | 0,1–1,8 s | ~240 |
| + retour arrière ciblé | 62,5 % | 0,8 s | 244 |
| + look-ahead par capacités | 85,5 % | 1,5 s | 244 |
| + ordre « étapes parfaites », 40 retours arrière | 293/300 (97,7 %) | 0,76 s (max 5,9 s) | 246 |
| + retour arrière sur 2 colonnes, sans relance | 200/200 | 0,47 s (max 5,7 s) | 246,5 |
| + relances (actuel) | 600/600, dont 596 du premier coup | 0,70 s (max 8,4 s) | 246,3 |
