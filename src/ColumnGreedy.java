import java.util.*;

/**
 * Greedy colonne par colonne.
 *
 * On fixe les attributs un par un, pour toute la population à la fois. Quand on traite
 * l'attribut a, chaque table t contenant a impose des comptes EXACTS sur la projection
 * de son scope sur les attributs déjà fixés + a (marginalisation de t sur le reste).
 * Ces contraintes ne portent que sur la colonne a : modifier la colonne a ne casse
 * jamais une colonne déjà validée.
 *
 * Pour une colonne :
 *  1. Glouton : on traite d'abord les individus les plus contraints (le moins de valeurs
 *     possibles) ; une valeur est possible si toutes ses cellules projetées ont du résiduel.
 *     Parmi les valeurs possibles, on écarte celles qui dépassent la capacité d'un
 *     attribut futur (capacités u_x des personas partielles), puis on favorise celle qui
 *     sépare l'individu de ses "jumeaux"
 *     (individus identiques sur les attributs déjà fixés) -> diversité.
 *  2. Micro-réparation (si impasse) : descente locale bornée sur la seule colonne a
 *     (changer une valeur, ou échanger deux valeurs), puis nouveaux tirages de la colonne.
 *  3. Si la colonne reste inexacte : retour arrière d'un niveau (on refait une colonne b déjà
 *     fixée, liée aux tables violées, puis a).
 */
public class ColumnGreedy {
    final Instance in;
    final Random rnd;
    public double beta = 1.0;         // poids de la diversité dans le choix de valeur
    public int repairIters = 500;     // budget de la micro-réparation par colonne
    public int[] order;               // ordre des attributs (null = heuristique)
    public boolean lookAhead = true;  // capacités des attributs futurs (borne u_x sur personas partielles)

    // diagnostics du dernier build
    public int columnsRepaired, columnsFailed;
    public final List<String> failLog = new ArrayList<>();   // "attribut@position:résiduel"
    public int[] lastOrder;

    final long[] radix;

    public ColumnGreedy(Instance in, Random rnd) {
        this.in = in; this.rnd = rnd;
        radix = new long[in.K];
        long r = 1;
        for (int a = in.K - 1; a >= 0; a--) { radix[a] = r; r *= in.dom[a]; }
    }

    /** Ordre heuristique : d'abord l'attribut le plus présent, puis celui qui partage le plus de tables avec les déjà choisis. */
    public int[] defaultOrder() {
        int K = in.K;
        int[] ord = new int[K];
        boolean[] done = new boolean[K];
        for (int s = 0; s < K; s++) {
            int best = -1; double bs = -Double.MAX_VALUE;
            for (int a = 0; a < K; a++) {
                if (done[a]) continue;
                double sc = 0;
                for (int t : in.attrTables[a]) {
                    int known = 0;
                    for (int b : in.scope[t]) if (done[b]) known++;
                    sc += 1 + 10.0 * known;    // relier aux colonnes déjà fixées
                }
                long cross = crossCells(done, a);
                sc += (cross == 0 ? perfectBonus : 0) - crossPenalty * cross;
                sc += rnd.nextDouble() * 1e-3;
                if (sc > bs) { bs = sc; best = a; }
            }
            ord[s] = best; done[best] = true;
        }
        return ord;
    }

    public double perfectBonus = 1000; // bonus d'ordre pour une étape parfaite (0 cellule croisée)
    public double crossPenalty = 0;    // pénalité par cellule croisée dans le choix de l'ordre

    /**
     * Cellules croisées de l'étape a : on prend les projections scope(t) ∩ connus des tables contenant a,
     * on garde les ensembles maximaux, et on somme le nombre de cellules de tous sauf le plus grand.
     * 0 = étape "parfaite" (une projection contient toutes les autres) : la colonne est alors un
     * transport exact groupe par groupe, toujours faisable.
     */
    long crossCells(boolean[] known, int a) {
        int[] sets = new int[in.attrTables[a].length];
        int m = 0;
        for (int t : in.attrTables[a]) {
            int s = 0;
            for (int b : in.scope[t]) if (known[b]) s |= 1 << b;
            if (s != 0) sets[m++] = s;
        }
        long sum = 0, max = 0;
        for (int i = 0; i < m; i++) {
            boolean maximal = true;
            for (int j = 0; j < m && maximal; j++)
                if (j != i && (sets[i] & sets[j]) == sets[i] && (sets[i] != sets[j] || j < i)) maximal = false;
            if (!maximal) continue;
            long cells = 1;
            for (int b = 0; b < in.K; b++) if ((sets[i] >> b & 1) != 0) cells *= in.dom[b];
            sum += cells; max = Math.max(max, cells);
        }
        return sum - max;
    }

    public int[][] build() {
        int N = in.N, K = in.K;
        int[] ord = order != null ? order : defaultOrder();
        int[][] pop = new int[N][K];
        boolean[] known = new boolean[K];
        columnsRepaired = 0; columnsFailed = 0; backtracks = 0; deepBacktracks = 0; failLog.clear(); lastOrder = ord;

        for (int pos = 0; pos < K; pos++) {
            int a = ord[pos];
            known[a] = true;
            long left = fillColumnRetry(pop, known, a);
            // Retour arrière : refaire 1, puis 2, ... jusqu'à backtrackDepth colonnes déjà fixées, puis a.
            // Avec a provisoirement inconnue, refaire une colonne b ne touche que les tables de b, qui
            // sont réimposées exactement : aucune autre colonne n'est cassée. Les colonnes refaites
            // sont choisies parmi les attributs des tables encore violées par a.
            if (left > 0 && pos > 0) {
                List<Integer> viol = lastViolated;
                for (int depth = 1; depth <= backtrackDepth && left > 0; depth++) {
                    for (int r = 0; r < backtrackTries && left > 0; r++) {
                        backtracks++;
                        int[] bs = pickColumns(viol, known, a, depth);
                        int[] saveA = column(pop, a);
                        int[][] saveB = new int[bs.length][];
                        for (int q = 0; q < bs.length; q++) saveB[q] = column(pop, bs[q]);
                        known[a] = false;
                        boolean ok = true;                                   // chaque colonne refaite doit rester exacte
                        for (int q = 0; q < bs.length && ok; q++) ok = fillColumnRetry(pop, known, bs[q]) == 0;
                        known[a] = true;
                        long l = ok ? fillColumnRetry(pop, known, a) : Long.MAX_VALUE;
                        if (l < left) { left = l; viol = lastViolated; if (depth > 1) deepBacktracks++; }
                        else {
                            setColumn(pop, a, saveA);
                            for (int q = 0; q < bs.length; q++) setColumn(pop, bs[q], saveB[q]);
                        }
                    }
                }
            }
            if (left > 0) { columnsFailed++; failLog.add(in.attrNames[a] + "@" + pos + ":" + left); }
        }
        return pop;
    }

    public int columnRetries = 3;     // nouveaux tirages d'une colonne qui échoue
    public int backtrackTries = 40;   // tentatives de retour arrière par profondeur
    public int backtrackDepth = 2;    // nombre max de colonnes refaites ensemble
    public int backtracks, deepBacktracks;

    /**
     * depth colonnes distinctes déjà fixées, tirées parmi les attributs des tables violées ;
     * si elles ne suffisent pas, parmi les attributs partageant une table avec a, puis tous les connus.
     */
    private int[] pickColumns(List<Integer> viol, boolean[] known, int a, int depth) {
        List<Integer> pool = new ArrayList<>();
        for (int t : viol) for (int b : in.scope[t]) if (b != a && known[b] && !pool.contains(b)) pool.add(b);
        if (pool.size() < depth)
            for (int t : in.attrTables[a]) for (int b : in.scope[t]) if (b != a && known[b] && !pool.contains(b)) pool.add(b);
        if (pool.size() < depth)
            for (int b = 0; b < in.K; b++) if (b != a && known[b] && !pool.contains(b)) pool.add(b);
        Collections.shuffle(pool, rnd);
        return pool.subList(0, Math.min(depth, pool.size())).stream().mapToInt(Integer::intValue).toArray();
    }
    private List<Integer> lastViolated = new ArrayList<>();   // tables violées par le dernier fillColumn

    private int[] column(int[][] pop, int a) {
        int[] c = new int[pop.length];
        for (int i = 0; i < pop.length; i++) c[i] = pop[i][a];
        return c;
    }

    private void setColumn(int[][] pop, int a, int[] c) {
        for (int i = 0; i < pop.length; i++) pop[i][a] = c[i];
    }

    /** fillColumn, relancée avec d'autres tirages tant qu'elle n'est pas exacte ; garde la meilleure. */
    private long fillColumnRetry(int[][] pop, boolean[] known, int a) {
        long best = fillColumn(pop, known, a);
        int[] bestCol = best > 0 ? column(pop, a) : null;
        List<Integer> bestViol = lastViolated;
        for (int r = 0; r < columnRetries && best > 0; r++) {
            long l = fillColumn(pop, known, a);
            if (l < best) { best = l; bestCol = column(pop, a); bestViol = lastViolated; }
        }
        if (bestCol != null) setColumn(pop, a, bestCol);
        lastViolated = bestViol;
        return best;
    }

    /**
     * Remplit la colonne a (known[a] = true, les autres colonnes connues sont figées).
     * Retourne l'erreur résiduelle sur les marges projetées (0 = colonne exacte).
     */
    private long fillColumn(int[][] pop, boolean[] known, int a) {
        int N = in.N;
        int[] tabs = in.attrTables[a];
        // marginales projetées : index de cellule avec les attributs inconnus mis à 0
        int[][] res = new int[tabs.length][];
        for (int k = 0; k < tabs.length; k++) {
            int t = tabs[k];
            res[k] = new int[in.target[t].length];
            for (int c = 0; c < in.target[t].length; c++) res[k][project(t, c, known)] += in.target[t][c];
        }
        int[][] base = new int[N][tabs.length];  // index projeté de chaque individu avec a = 0
        int[] strideA = new int[tabs.length];
        for (int k = 0; k < tabs.length; k++) {
            int t = tabs[k];
            for (int j = 0; j < in.scope[t].length; j++) if (in.scope[t][j] == a) strideA[k] = in.stride[t][j];
        }
        long[] baseKey = new long[N];            // code des attributs connus avec a = 0 (aussi signature des jumeaux)
        for (int i = 0; i < N; i++) {
            pop[i][a] = 0;
            for (int k = 0; k < tabs.length; k++) base[i][k] = projectInd(tabs[k], pop[i], known);
            baseKey[i] = knownKey(pop[i], known);
        }

        // look-ahead entier (borne u_x sur les personas partielles) : pour chaque attribut futur c,
        // un groupe g (valeurs connues sur les attributs liés à c) contient au plus
        // cap(g) = somme_w min_t cible_t(g, w) individus. Chaque c devient une table virtuelle de capacités.
        int d = in.dom[a];
        Ahead ah = lookAhead ? lookAhead(pop, known, a) : null;
        prepareColumn(known, a);

        // --- 1. glouton : individus les plus contraints d'abord ---
        boolean[] done = new boolean[N];
        HashMap<Long, int[]> twins = new HashMap<>();   // signature -> nb déjà affectés par valeur
        int[] val = new int[N];
        boolean failed = false;
        for (int step = 0; step < N; step++) {
            int bi = -1, bf = Integer.MAX_VALUE, ties = 0;
            for (int i = 0; i < N; i++) {
                if (done[i]) continue;
                int f = 0;
                for (int v = 0; v < d; v++)
                    if (feasible(res, base[i], strideA, v) && (ah == null || ah.ok(i, v))) f++;
                if (f < bf) { bf = f; bi = i; ties = 1; }
                else if (f == bf && rnd.nextInt(++ties) == 0) bi = i;
            }
            int i = bi;
            int[] tw = twins.computeIfAbsent(baseKey[i], s -> new int[d]);
            int bv = -1; double bs = -Double.MAX_VALUE;
            for (int v = 0; v < d; v++) {
                int mr = minRes(res, base[i], strideA, v);
                double s = (mr > 0 ? 1e6 : 0) + (ah == null || ah.ok(i, v) ? 1e5 : 0)
                        + diversity(mr, tw, v, baseKey[i] + v * radix[a]);
                if (s > bs) { bs = s; bv = v; }
            }
            if (minRes(res, base[i], strideA, bv) <= 0) failed = true;
            val[i] = bv; done[i] = true; tw[bv]++;
            if (ah != null) ah.add(i, bv);
            for (int k = 0; k < tabs.length; k++) res[k][base[i][k] + bv * strideA[k]]--;
        }


        // --- 2. micro-réparation sur la colonne a ---
        long left = 0;
        if (failed) {
            columnsRepaired++;
            left = repair(res, base, strideA, val, d);
        }
        lastViolated = new ArrayList<>();
        for (int k = 0; k < tabs.length; k++)
            for (int x : res[k]) if (x != 0) { lastViolated.add(tabs[k]); break; }
        for (int i = 0; i < N; i++) pop[i][a] = val[i];
        return left;
    }


    /** Appelé au début de chaque remplissage de la colonne a (known[a] = true). */
    protected void prepareColumn(boolean[] known, int a) {}

    /**
     * Préférence (entre 0 et 1e5) pour la valeur v d'un individu, parmi les valeurs de même faisabilité.
     * mr = plus petit résiduel des cellules de v, tw[v] = jumeaux déjà envoyés sur v, branch = code des
     * attributs connus avec a = v. Par défaut : séparer les jumeaux (mr / (1 + tw[v])^beta).
     */
    protected double diversity(int mr, int[] tw, int v, long branch) {
        return mr / Math.pow(1 + tw[v], beta) * (1 + 0.1 * rnd.nextDouble());
    }

    /** Tables virtuelles de capacité, une par attribut futur dont les groupes dépendent de a. */
    private static final class Ahead {
        int[][][] gid;   // gid[c][i][v] : groupe de l'individu i si a = v
        int[][] cap, cnt;
        boolean ok(int i, int v) {
            for (int c = 0; c < gid.length; c++) { int g = gid[c][i][v]; if (cnt[c][g] >= cap[c][g]) return false; }
            return true;
        }
        void add(int i, int v) { for (int c = 0; c < gid.length; c++) cnt[c][gid[c][i][v]]++; }
    }

    private Ahead lookAhead(int[][] pop, boolean[] known, int a) {
        int N = in.N, d = in.dom[a];
        List<int[][]> gids = new ArrayList<>();
        List<int[]> caps = new ArrayList<>();
        boolean[] kc = known.clone();
        for (int c = 0; c < in.K; c++) {
            if (known[c]) continue;
            // U = attributs connus (a compris) qui partagent une table avec c
            boolean[] U = new boolean[in.K];
            for (int t : in.attrTables[c]) for (int b : in.scope[t]) if (known[b]) U[b] = true;
            if (!U[a]) continue;                       // les groupes de c ne dépendent pas de a
            kc[c] = true;
            int[] tabs = in.attrTables[c];
            int[][] pt = new int[tabs.length][];       // cibles projetées sur scope(t) ∩ (connus + c)
            for (int k = 0; k < tabs.length; k++) {
                int t = tabs[k];
                pt[k] = new int[in.target[t].length];
                for (int cell = 0; cell < in.target[t].length; cell++) pt[k][project(t, cell, kc)] += in.target[t][cell];
            }
            HashMap<Long, Integer> ids = new HashMap<>();
            List<Integer> capList = new ArrayList<>();
            int[][] gid = new int[N][d];
            for (int i = 0; i < N; i++) {
                int[] x = pop[i];
                int saveA = x[a], saveC = x[c];
                for (int v = 0; v < d; v++) {
                    x[a] = v;
                    long key = 0;
                    for (int b = 0; b < in.K; b++) if (U[b]) key += x[b] * radix[b];
                    Integer id = ids.get(key);
                    if (id == null) {
                        int cp = 0;
                        for (int w = 0; w < in.dom[c]; w++) {
                            x[c] = w;
                            int m = Integer.MAX_VALUE;
                            for (int k = 0; k < tabs.length; k++) m = Math.min(m, pt[k][projectInd(tabs[k], x, kc)]);
                            cp += m;
                        }
                        id = capList.size();
                        ids.put(key, id);
                        capList.add(cp);
                    }
                    gid[i][v] = id;
                }
                x[a] = saveA; x[c] = saveC;
            }
            kc[c] = false;
            gids.add(gid);
            caps.add(capList.stream().mapToInt(Integer::intValue).toArray());
        }
        Ahead ah = new Ahead();
        ah.gid = gids.toArray(new int[0][][]);
        ah.cap = caps.toArray(new int[0][]);
        ah.cnt = new int[ah.cap.length][];
        for (int c = 0; c < ah.cap.length; c++) ah.cnt[c] = new int[ah.cap[c].length];
        return ah;
    }


    private boolean feasible(int[][] res, int[] base, int[] strideA, int v) {
        for (int k = 0; k < res.length; k++) if (res[k][base[k] + v * strideA[k]] <= 0) return false;
        return true;
    }

    private int minRes(int[][] res, int[] base, int[] strideA, int v) {
        int m = Integer.MAX_VALUE;
        for (int k = 0; k < res.length; k++) m = Math.min(m, res[k][base[k] + v * strideA[k]]);
        return m;
    }

    /** Descente locale sur la colonne : erreur = somme |res|. */
    private long repair(int[][] res, int[][] base, int[] strideA, int[] val, int d) {
        int N = in.N;
        long err = 0;
        for (int[] r : res) for (int x : r) err += Math.abs(x);
        for (int it = 0; it < repairIters && err > 0; it++) {
            // individu dans une cellule en excès (res < 0)
            int i = rnd.nextInt(N), tries = 0;
            while (tries++ < 200 && !inOverfull(res, base[i], strideA, val[i])) i = rnd.nextInt(N);
            // meilleur mouvement : changer val[i], ou échanger avec j
            long bestDelta = Long.MAX_VALUE; int bj = -1, bv = -1;
            for (int v = 0; v < d; v++) {
                if (v == val[i]) continue;
                long dl = moveDelta(res, base[i], strideA, val[i], v);
                if (dl < bestDelta) { bestDelta = dl; bv = v; bj = -1; }
            }
            for (int j = 0; j < N; j++) {
                if (val[j] == val[i]) continue;
                long dl = swapDelta(res, base, strideA, val, i, j);
                if (dl < bestDelta || (dl == bestDelta && rnd.nextInt(4) == 0)) { bestDelta = dl; bj = j; bv = -1; }
            }
            if (bestDelta > 0 || (bestDelta == 0 && rnd.nextInt(3) != 0)) continue;   // plateau : déplacement latéral parfois
            if (bj < 0) { applyMove(res, base[i], strideA, val[i], bv); val[i] = bv; }
            else {
                int vi = val[i], vj = val[bj];
                applyMove(res, base[i], strideA, vi, vj); applyMove(res, base[bj], strideA, vj, vi);
                val[i] = vj; val[bj] = vi;
            }
            err += bestDelta;
        }
        return err;
    }

    private boolean inOverfull(int[][] res, int[] base, int[] strideA, int v) {
        for (int k = 0; k < res.length; k++) if (res[k][base[k] + v * strideA[k]] < 0) return true;
        return false;
    }

    /** Variation de somme |res| si un individu passe de v0 à v1 (res = cible - compte). */
    private long moveDelta(int[][] res, int[] base, int[] strideA, int v0, int v1) {
        long dl = 0;
        for (int k = 0; k < res.length; k++) {
            int c0 = base[k] + v0 * strideA[k], c1 = base[k] + v1 * strideA[k];
            if (c0 == c1) continue;
            dl += Math.abs(res[k][c0] + 1) - Math.abs(res[k][c0]);
            dl += Math.abs(res[k][c1] - 1) - Math.abs(res[k][c1]);
        }
        return dl;
    }

    private long swapDelta(int[][] res, int[][] base, int[] strideA, int[] val, int i, int j) {
        int vi = val[i], vj = val[j];
        long d1 = moveDelta(res, base[i], strideA, vi, vj);
        applyMove(res, base[i], strideA, vi, vj);
        long d2 = moveDelta(res, base[j], strideA, vj, vi);
        applyMove(res, base[i], strideA, vj, vi);
        return d1 + d2;
    }

    private void applyMove(int[][] res, int[] base, int[] strideA, int v0, int v1) {
        for (int k = 0; k < res.length; k++) {
            res[k][base[k] + v0 * strideA[k]]++;
            res[k][base[k] + v1 * strideA[k]]--;
        }
    }

    /** Code d'une persona restreinte aux attributs connus (les autres comptent pour 0). */
    protected long knownKey(int[] x, boolean[] known) {
        long k = 0;
        for (int b = 0; b < in.K; b++) if (known[b]) k += x[b] * radix[b];
        return k;
    }

    /** Index de cellule c de la table t avec les chiffres des attributs non connus mis à 0. */
    private int project(int t, int c, boolean[] known) {
        int idx = 0;
        for (int j = 0; j < in.scope[t].length; j++) {
            int b = in.scope[t][j];
            if (known[b]) idx += ((c / in.stride[t][j]) % in.dom[b]) * in.stride[t][j];
        }
        return idx;
    }

    /** Index projeté d'un individu (attributs connus, a compris, avec pop[i][a] courant). */
    private int projectInd(int t, int[] x, boolean[] known) {
        int idx = 0;
        for (int j = 0; j < in.scope[t].length; j++) if (known[in.scope[t][j]]) idx += x[in.scope[t][j]] * in.stride[t][j];
        return idx;
    }
}
