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
    public double small = 0;          // individu de départ dans un paquet de taille m accepté avec proba (2/m)^small
    public long iterations, accepted; // diagnostics

    final long[] radix;
    final int[][] sAttrs, bAttrs;     // pour chaque ensemble S retenu : ses attributs, sa frontière

    public Declone(Instance in, Random rnd) {
        this.in = in; this.rnd = rnd;
        int K = in.K;
        radix = new long[K];
        long r = 1;
        for (int a = K - 1; a >= 0; a--) { radix[a] = r; r *= in.dom[a]; }
        int[] nb = new int[K];                        // voisinage de chaque attribut (masque)
        for (int t = 0; t < in.T; t++) {
            int m = 0;
            for (int b : in.scope[t]) m |= 1 << b;
            for (int b : in.scope[t]) nb[b] |= m;
        }
        // Échanger S ou son complément revient au même : on prend les S sans le dernier attribut.
        // Le mouvement n'a d'effet que si un attribut reste hors de S ∪ B(S).
        List<int[]> ss = new ArrayList<>(), bs = new ArrayList<>();
        int full = (1 << K) - 1;
        List<Integer> masks = new ArrayList<>();
        if (K <= 20) for (int s = 1; s < 1 << (K - 1); s++) masks.add(s);
        else for (int q = 0; q < 1 << 16; q++) {      // trop d'attributs : ensembles connexes tirés au hasard
            int s = 1 << rnd.nextInt(K), size = 1 + rnd.nextInt(K / 2);
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

    /** Améliore pop en place (et la remplace par la meilleure trouvée) ; renvoie le nombre de distinctes. */
    public int run(int[][] pop) {
        int N = pop.length;
        if (sAttrs.length == 0) return distinct(pop);
        long[] full = new long[N];
        HashMap<Long, Integer> cnt = new HashMap<>();
        for (int i = 0; i < N; i++) { full[i] = key(pop[i]); cnt.merge(full[i], 1, Integer::sum); }
        int cur = cnt.size(), best = cur;
        int[][] bestPop = copy(pop);
        long start = System.nanoTime(), end = start + timeMs * 1_000_000L;
        double T = temp;
        int[] cand = new int[N];
        long[] candI = new long[N], candJ = new long[N];
        iterations = 0; accepted = 0;
        while (true) {
            if ((iterations & 255) == 0) {
                long now = System.nanoTime();
                if (now >= end) break;
                T = temp * (end - now) / (double) (end - start);
            }
            iterations++;
            if (cur == N) break;
            int i = rnd.nextInt(N);                   // individu en doublon, de préférence dans un petit paquet
            for (int tries = 0; tries < 4 * N; tries++) {
                int m = cnt.get(full[i]);
                if (m >= 2 && (small == 0 || rnd.nextDouble() < Math.pow(2.0 / m, small))) break;
                i = rnd.nextInt(N);
            }
            int q = rnd.nextInt(sAttrs.length);
            int[] S = sAttrs[q], B = bAttrs[q];
            long bi = key(pop[i], B), si = key(pop[i], S);
            // meilleur partenaire j : même frontière, différent sur S et hors de S ∪ B
            long bestScore = Long.MIN_VALUE; int nc = 0;
            for (int j = 0; j < N; j++) {
                if (j == i || key(pop[j], B) != bi) continue;
                long sj = key(pop[j], S);
                if (sj == si) continue;
                long ni = full[i] - si + sj, nj = full[j] - sj + si;
                if (ni == full[j]) continue;          // i et j identiques hors de S : simple permutation
                long sc = score(cnt, full[i], full[j], ni, nj);
                if (sc > bestScore) { bestScore = sc; nc = 0; }
                if (sc == bestScore) { cand[nc] = j; candI[nc] = ni; candJ[nc] = nj; nc++; }
            }
            if (nc == 0) continue;
            int bestDelta = (int) Math.floorDiv(bestScore + SCALE / 2, SCALE);
            long dsq = bestScore - bestDelta * SCALE;
            if (bestDelta == 0 && dsq <= 0 && rnd.nextDouble() >= neutral) continue;
            if (bestDelta < 0 && (T <= 0 || rnd.nextDouble() >= Math.exp(bestDelta / T))) continue;
            int c = rnd.nextInt(nc), j = cand[c];
            apply(cnt, full[i], full[j], candI[c], candJ[c]);
            for (int b : S) { int v = pop[i][b]; pop[i][b] = pop[j][b]; pop[j][b] = v; }
            full[i] = candI[c]; full[j] = candJ[c];
            cur += bestDelta; accepted++;
            if (cur > best) { best = cur; bestPop = copy(pop); }
        }
        for (int i = 0; i < N; i++) pop[i] = bestPop[i];
        return best;
    }

    static final long SCALE = 1L << 32;
    private long sq;                  // variation de Σ m² accumulée par dec / inc

    /**
     * Score si (oi, oj) devient (ni, nj) : variation des distinctes × SCALE, plus (si concentrate)
     * la variation de Σ m² ; la table est remise en l'état.
     */
    private long score(HashMap<Long, Integer> cnt, long oi, long oj, long ni, long nj) {
        sq = 0;
        int d = apply(cnt, oi, oj, ni, nj);
        long s2 = sq;
        apply(cnt, ni, nj, oi, oj);
        return d * SCALE + (concentrate ? s2 : 0);
    }

    private int apply(HashMap<Long, Integer> cnt, long oi, long oj, long ni, long nj) {
        return dec(cnt, oi) + dec(cnt, oj) + inc(cnt, ni) + inc(cnt, nj);
    }

    private int dec(HashMap<Long, Integer> cnt, long k) {
        int c = cnt.get(k);
        sq -= 2L * c - 1;
        if (c == 1) { cnt.remove(k); return -1; }
        cnt.put(k, c - 1);
        return 0;
    }

    private int inc(HashMap<Long, Integer> cnt, long k) {
        Integer c = cnt.get(k);
        sq += c == null ? 1 : 2L * c + 1;
        cnt.put(k, c == null ? 1 : c + 1);
        return c == null ? 1 : 0;
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
