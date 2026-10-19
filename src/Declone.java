import java.util.*;

/**
 * Post-optimisation qui conserve l'exactitude : réduit les clones par échanges partiels.
 *
 * Mouvement : deux individus i et j, un ensemble d'attributs S. On échange leurs valeurs sur S.
 * Soit B(S) la frontière de S (attributs hors S qui partagent une table avec un attribut de S).
 * Si i et j coïncident sur B(S), toute table est soit dans S ∪ B(S), soit hors de S : ses comptes
 * sont inchangés. Une population exacte reste donc exacte à chaque mouvement.
 *
 * Objectif : nombre de personas distinctes. Recuit : gains nuls acceptés avec probabilité neutral, pertes avec
 * probabilité exp(gain / T), T décroissant linéairement de temp à 0 sur le budget (temp = 0 : simple montée) ;
 * on garde la meilleure population rencontrée.
 */
public final class Declone {
    final Instance in;
    final Random rnd;
    public long timeMs = 3000;        // budget de temps
    public double neutral = 0.3;      // probabilité d'accepter un mouvement de gain nul
    public double temp = 0;           // température initiale du recuit (0 = pas de mouvement perdant)
    public boolean concentrate = false; // à distinctes égales, préférer les échanges qui concentrent les doublons (Σ m² croît)
    public boolean spread = false;    // à distinctes égales, préférer les échanges qui étalent les doublons (Σ m² décroît)
    public double small = 0;          // individu de départ dans un paquet de taille m accepté avec proba (2/m)^small
    /**
     * Objectif borné (population partielle) : cap[g] = nombre de personas actives qui prolongent la persona
     * partielle g. Un paquet de m individus en g donnera au plus min(m, cap[g]) distinctes : on maximise
     * Σ_g min(m_g, cap[g]) au lieu du nombre de distinctes (null = distinctes).
     */
    public HashMap<Long, Integer> cap;
    /**
     * Critère secondaire, à distinctes égales (entropie de Shannon, en nats) : 0 = aucun, 1 = H(X),
     * 2 = Σ_b H(X_{-b}) sur les attributs connus, 3 = -Σ_b H(X_b | X_{-b}) = Σ_b H(X_{-b}) - K·H(X) (minimise les
     * quasi-clones : H(X_b | X_{-b}) = 0 ssi aucune paire ne diffère seulement sur b).
     * On garde aussi la meilleure population au sens (distinctes, critère secondaire).
     */
    public int entropy = 0;
    public long iterations, accepted; // diagnostics

    final long[] radix;
    final int[][] sAttrs, bAttrs;     // pour chaque ensemble S retenu : ses attributs, sa frontière
    final int[] knownAttrs;           // attributs connus (masque du constructeur)

    public Declone(Instance in, Random rnd) { this(in, rnd, (1 << in.K) - 1); }

    /**
     * Dé-clonage restreint aux attributs connus (masque known), pour une population partielle dont les autres
     * colonnes valent 0 (construction colonne par colonne). Mouvements : S ⊆ known, i et j égaux sur
     * B(S) ∩ known. Chaque table garde les comptes de sa projection sur les attributs connus : une population
     * partielle exacte le reste.
     */
    public Declone(Instance in, Random rnd, int known) {
        this.in = in; this.rnd = rnd;
        int K = in.K;
        knownAttrs = bits(known);
        radix = new long[K];
        long r = 1;
        for (int a = K - 1; a >= 0; a--) { radix[a] = r; r *= in.dom[a]; }
        int[] nb = new int[K];                        // voisinage de chaque attribut, parmi les connus (masque)
        for (int t = 0; t < in.T; t++) {
            int m = 0;
            for (int b : in.scope[t]) m |= 1 << b;
            for (int b : in.scope[t]) nb[b] |= m & known;
        }
        // Échanger S ou son complément (dans known) revient au même : on prend les S sans le dernier attribut connu.
        // Le mouvement n'a d'effet que si un attribut connu reste hors de S ∪ B(S).
        List<int[]> ss = new ArrayList<>(), bs = new ArrayList<>();
        int full = known, rest = known & ~Integer.highestOneBit(known), nk = Integer.bitCount(known);
        List<Integer> masks = new ArrayList<>();
        if (K <= 20) { for (int s = 1; s <= rest; s++) if ((s & ~rest) == 0) masks.add(s); }
        else for (int q = 0; q < 1 << 16; q++) {      // trop d'attributs : ensembles connexes tirés au hasard
            int s = 0, size = 1 + rnd.nextInt(Math.max(1, nk / 2));
            for (int p = rnd.nextInt(nk), b = 0; b < K; b++) if ((known >> b & 1) != 0 && p-- == 0) { s = 1 << b; break; }
            while (Integer.bitCount(s) < size) {
                int cand = 0;
                for (int b = 0; b < K; b++) if ((s >> b & 1) != 0) cand |= nb[b];
                cand &= ~s;
                if (cand == 0) break;
                int pick = rnd.nextInt(Integer.bitCount(cand));
                for (int b = 0; b < K; b++) if ((cand >> b & 1) != 0 && pick-- == 0) { s |= 1 << b; break; }
            }
            masks.add(s);
        }
        for (int s : masks) {
            int n = 0;
            for (int b = 0; b < K; b++) if ((s >> b & 1) != 0) n |= nb[b];
            int bnd = n & ~s;
            if ((s | bnd) == full) continue;
            ss.add(bits(s)); bs.add(bits(bnd));
        }
        sAttrs = ss.toArray(new int[0][]);
        bAttrs = bs.toArray(new int[0][]);
    }

    private int[] bits(int m) {
        int[] r = new int[Integer.bitCount(m)];
        int k = 0;
        for (int b = 0; b < in.K; b++) if ((m >> b & 1) != 0) r[k++] = b;
        return r;
    }

    private long key(int[] x, int[] attrs) {
        long k = 0;
        for (int b : attrs) k += x[b] * radix[b];
        return k;
    }

    private long key(int[] x) {
        long k = 0;
        for (int b = 0; b < in.K; b++) k += x[b] * radix[b];
        return k;
    }

    /**
     * Rayon des projections du critère loo : on utilise X_{-B} pour tous les B ⊆ connus de taille radius.
     * H(X_{-B}) = ln N ssi les individus restent tous distincts quand on ignore B : Σ_B H(X_{-B}) est maximal
     * ssi aucune paire d'individus n'est à distance ≤ radius.
     */
    public int radius = 1;
    /**
     * Poids des distinctes dans un objectif pondéré F = mu × distinctes + N × (critère d'entropie) (mu ≥ 0) ;
     * mu < 0 : lexicographique (distinctes d'abord, entropie à égalité).
     */
    public double mu = -1;
    public double cycle = 0;          // probabilité de tenter un cycle à trois individus i <- j <- k <- i sur S
    public long cycles;               // diagnostics : cycles acceptés

    static final double LEX = 1e6;    // poids des distinctes en mode lexicographique (> 2 × |critère secondaire| d'un mouvement)

    /** Améliore pop en place (et la remplace par la meilleure trouvée) ; renvoie le nombre de distinctes. */
    public int run(int[][] pop) {
        int N = pop.length;
        if (sAttrs.length == 0) return distinct(pop);
        long[] full = new long[N];
        HashMap<Long, Integer> cnt = new HashMap<>();
        for (int i = 0; i < N; i++) { full[i] = key(pop[i]); cnt.merge(full[i], 1, Integer::sum); }
        int cur = 0;
        for (Map.Entry<Long, Integer> e : cnt.entrySet()) cur += value(e.getKey(), e.getValue());
        if (entropy >= 2) {
            masks = new ArrayList<>(); loo = new ArrayList<>();
            subsets(0, new int[Math.min(radius, knownAttrs.length)], 0);
            for (int[] m : masks) {
                HashMap<Long, Integer> h = new HashMap<>();
                for (int i = 0; i < N; i++) h.merge(masked(full[i], m), 1, Integer::sum);
                loo.add(h);
            }
        }
        W = mu < 0 ? LEX : mu;
        double obj = 0, bestObj = 0;      // objectif cumulé (relatif au départ)
        int best = cur;
        int[][] bestPop = copy(pop);
        long start = System.nanoTime(), end = start + timeMs * 1_000_000L;
        double T = temp;
        int[] cand = new int[N], candK = new int[N];
        iterations = 0; accepted = 0; cycles = 0;
        long[] o = new long[3], n = new long[3];
        while (true) {
            if ((iterations & 255) == 0) {
                long now = System.nanoTime();
                if (now >= end) break;
                T = temp * (end - now) / (double) (end - start);
            }
            iterations++;
            if (cur == N && entropy == 0) break;
            int i = rnd.nextInt(N);                   // individu en doublon, de préférence dans un petit paquet
            for (int tries = 0; tries < 4 * N; tries++) {
                int m = cnt.get(full[i]);
                if (m >= 2 && (small == 0 || rnd.nextDouble() < Math.pow(2.0 / m, small))) break;
                i = rnd.nextInt(N);
            }
            int q = rnd.nextInt(sAttrs.length);
            int[] S = sAttrs[q], B = bAttrs[q];
            long bi = key(pop[i], B), si = key(pop[i], S);
            boolean cyc = cycle > 0 && rnd.nextDouble() < cycle;
            double bestScore = Double.NEGATIVE_INFINITY; int nc = 0, bestD = 0;
            if (!cyc) {
                // meilleur partenaire j : même frontière, différent sur S et hors de S ∪ B
                for (int j = 0; j < N; j++) {
                    if (j == i || key(pop[j], B) != bi) continue;
                    long sj = key(pop[j], S);
                    if (sj == si) continue;
                    o[0] = full[i]; o[1] = full[j];
                    n[0] = full[i] - si + sj; n[1] = full[j] - sj + si;
                    if (n[0] == full[j]) continue;    // i et j identiques hors de S : simple permutation
                    double sc = score(cnt, o, n, 2);
                    if (sc > bestScore) { bestScore = sc; nc = 0; bestD = lastD; }
                    if (sc == bestScore) cand[nc++] = j;
                }
            } else {
                // cycle i <- j <- k <- i (i prend S_j, j prend S_k, k prend S_i) : composé de deux échanges valides,
                // il franchit d'un coup un état intermédiaire moins bon. Quelques j tirés, le meilleur k pour chacun.
                List<Integer> g = new ArrayList<>();
                for (int j = 0; j < N; j++) if (j != i && key(pop[j], B) == bi && key(pop[j], S) != si) g.add(j);
                if (g.size() < 2) continue;
                for (int t = 0; t < Math.min(8, g.size()); t++) {
                    int j = g.get(rnd.nextInt(g.size()));
                    long sj = key(pop[j], S);
                    for (int k : g) {
                        if (k == j) continue;
                        long sk = key(pop[k], S);
                        if (sk == sj) continue;           // trois valeurs de S distinctes, sinon c'est un échange
                        o[0] = full[i]; o[1] = full[j]; o[2] = full[k];
                        n[0] = full[i] - si + sj; n[1] = full[j] - sj + sk; n[2] = full[k] - sk + si;
                        double sc = score(cnt, o, n, 3);
                        if (sc > bestScore) { bestScore = sc; nc = 0; bestD = lastD; }
                        if (sc == bestScore && nc < N) { cand[nc] = j; candK[nc] = k; nc++; }
                    }
                }
            }
            if (nc == 0) continue;
            // acceptation : gain > 0 toujours ; gain nul (lexicographique : distinctes inchangées et critère
            // secondaire ≤ 0) avec probabilité neutral ; perte avec probabilité exp(perte / T), en unités de distinctes.
            boolean neutralMove = mu < 0 ? (bestD == 0 && bestScore <= 0) : Math.abs(bestScore) < 1e-9;
            if (neutralMove) { if (rnd.nextDouble() >= neutral) continue; }
            else if (bestScore < 0 && (T <= 0 || rnd.nextDouble() >= Math.exp(bestScore / W / T))) continue;
            int c = rnd.nextInt(nc), j = cand[c];
            if (!cyc) {
                long sj = key(pop[j], S);
                o[0] = full[i]; o[1] = full[j]; n[0] = full[i] - si + sj; n[1] = full[j] - sj + si;
                apply(cnt, o, n, 2);
                for (int b : S) { int v = pop[i][b]; pop[i][b] = pop[j][b]; pop[j][b] = v; }
                full[i] = n[0]; full[j] = n[1];
            } else {
                int k = candK[c];
                long sj = key(pop[j], S), sk = key(pop[k], S);
                o[0] = full[i]; o[1] = full[j]; o[2] = full[k];
                n[0] = full[i] - si + sj; n[1] = full[j] - sj + sk; n[2] = full[k] - sk + si;
                apply(cnt, o, n, 3);
                for (int b : S) { int v = pop[i][b]; pop[i][b] = pop[j][b]; pop[j][b] = pop[k][b]; pop[k][b] = v; }
                full[i] = n[0]; full[j] = n[1]; full[k] = n[2];
                cycles++;
            }
            cur += bestD; obj += bestScore; accepted++;
            boolean better = mu < 0 && entropy == 0 && !concentrate && !spread ? cur > best : obj > bestObj + 1e-9;
            if (better) { best = cur; bestObj = obj; bestPop = copy(pop); }
        }
        for (int i = 0; i < N; i++) pop[i] = bestPop[i];
        return distinct(pop);
    }

    private double W;                 // poids des distinctes dans le score
    private int lastD;                // variation des distinctes du dernier score
    private long sq;                  // variation de Σ m² accumulée par dec / inc
    private double nlnF, nlnL;        // variations de Σ n ln n sur les cellules de X, et de tous les X_{-B}, accumulées par apply
    private List<int[]> masks;        // entropy ≥ 2 : ensembles B (taille radius) des projections X_{-B}
    private List<HashMap<Long, Integer>> loo;   // comptes des cellules de chaque X_{-B}
    private double[] deltaTab;

    private void subsets(int from, int[] cur, int size) {
        if (size == cur.length) { masks.add(cur.clone()); return; }
        for (int q = from; q < knownAttrs.length; q++) { cur[size] = knownAttrs[q]; subsets(q + 1, cur, size + 1); }
    }

    /** δ(n) = (n+1) ln(n+1) - n ln n : hausse de Σ n ln n quand une cellule de n individus en reçoit un de plus. */
    private double delta(int n) {
        if (deltaTab == null || n >= deltaTab.length) {
            deltaTab = new double[Math.max(n + 1, in.N + 2)];
            for (int k = 0; k < deltaTab.length; k++) deltaTab[k] = (k + 1) * Math.log(k + 1) - (k == 0 ? 0 : k * Math.log(k));
        }
        return deltaTab[n];
    }

    /** Code de la persona k privée des attributs de B (chiffres mis à 0). */
    private long masked(long k, int[] B) {
        for (int b : B) k -= (k / radix[b] % in.dom[b]) * radix[b];
        return k;
    }

    private void bump(HashMap<Long, Integer> h, long k, int d) {
        int c = h.getOrDefault(k, 0);
        nlnL += d > 0 ? delta(c) : -delta(c - 1);
        if (c + d == 0) h.remove(k); else h.put(k, c + d);
    }

    /**
     * Score si les personas o[0..m) deviennent n[0..m) : W × variation des distinctes, plus (si concentrate) la variation
     * de Σ m², ou moins (si spread), plus (si entropy) N × la variation d'entropie (= -Δ Σ n ln n) ; tables remises en l'état.
     */
    private double score(HashMap<Long, Integer> cnt, long[] o, long[] n, int m) {
        sq = 0; nlnF = 0; nlnL = 0;
        int d = apply(cnt, o, n, m);
        long s2 = sq;
        double e = entropy == 1 ? nlnF : entropy == 2 ? nlnL : nlnL - masks.size() * nlnF;
        apply(cnt, n, o, m);
        lastD = d;
        return d * W + (concentrate ? s2 : spread ? -s2 : 0) + (entropy > 0 ? -e : 0);
    }

    private int apply(HashMap<Long, Integer> cnt, long[] o, long[] n, int m) {
        int d = 0;
        for (int t = 0; t < m; t++) d += dec(cnt, o[t]);
        for (int t = 0; t < m; t++) d += inc(cnt, n[t]);
        if (entropy >= 2)
            for (int q = 0; q < masks.size(); q++) {
                HashMap<Long, Integer> h = loo.get(q); int[] M = masks.get(q);
                for (int t = 0; t < m; t++) bump(h, masked(o[t], M), -1);
                for (int t = 0; t < m; t++) bump(h, masked(n[t], M), +1);
            }
        return d;
    }

    /** Valeur d'une persona de multiplicité m : 1 si présente, ou min(m, cap) si cap est donné. */
    private int value(long k, int m) {
        return cap == null ? (m > 0 ? 1 : 0) : Math.min(m, cap.getOrDefault(k, 0));
    }

    private int dec(HashMap<Long, Integer> cnt, long k) {
        int c = cnt.get(k);
        sq -= 2L * c - 1;
        nlnF -= delta(c - 1);
        if (c == 1) cnt.remove(k); else cnt.put(k, c - 1);
        return value(k, c - 1) - value(k, c);
    }

    private int inc(HashMap<Long, Integer> cnt, long k) {
        Integer c0 = cnt.get(k);
        int c = c0 == null ? 0 : c0;
        sq += 2L * c + 1;
        nlnF += delta(c);
        cnt.put(k, c + 1);
        return value(k, c + 1) - value(k, c);
    }

    private int distinct(int[][] pop) {
        HashSet<Long> h = new HashSet<>();
        for (int[] x : pop) h.add(key(x));
        return h.size();
    }

    private static int[][] copy(int[][] pop) {
        int[][] c = new int[pop.length][];
        for (int i = 0; i < pop.length; i++) c[i] = pop[i].clone();
        return c;
    }
}
