import java.nio.file.*;
import java.util.*;

/**
 * Diagnostic : existe-t-il, sur une population finale, des mouvements de degré 4 « double échange compensé » qui
 * améliorent le score lexicographique (distinctes, puis Σ_b H(X_{-b})) alors qu'aucun échange valide ne le fait ?
 *
 * Double échange : (i, j) échangent S alors qu'ils diffèrent sur F ⊆ B(S) (1 ≤ |F| ≤ fmax) : seul, il casse les tables
 * qui contiennent un attribut de F et touchent S. On cherche une seconde paire (k, l), échangée sur le même S, dont
 * l'écart sur les tables est exactement l'opposé. Écart d'une paire = Σ_t [e(i'_t) + e(j'_t) - e(i_t) - e(j_t)],
 * haché linéairement (Zobrist additif) : δ1 + δ2 = 0 se cherche en O(1) ; chaque mouvement trouvé est revérifié
 * table par table.
 *
 * Usage : java -cp out Deg4Diag <contraintes.csv> <fmax> <population.csv>...
 */
public final class Deg4Diag {
    static Instance in;
    static long[] radix;
    static int K, N;

    public static void main(String[] args) throws Exception {
        in = Instance.load(Paths.get(args[0]));
        int fmax = Integer.parseInt(args[1]);
        K = in.K; N = in.N;
        radix = new long[K];
        long r = 1;
        for (int a = K - 1; a >= 0; a--) { radix[a] = r; r *= in.dom[a]; }
        for (int f = 2; f < args.length; f++) {
            long t0 = System.nanoTime();
            diag(CheckPopulation.read(in, Paths.get(args[f])), fmax);
            System.out.printf(Locale.ROOT, "  (%s, %.1f s)%n", args[f], (System.nanoTime() - t0) / 1e9);
        }
    }

    static long key(int[] x) { long k = 0; for (int b = 0; b < K; b++) k += x[b] * radix[b]; return k; }

    static long masked(long k, int b) { return k - (k / radix[b] % in.dom[b]) * radix[b]; }

    static HashMap<Long, Integer> cnt;
    static List<HashMap<Long, Integer>> loo;

    static double dl(int n) { return n == 0 ? 0 : n * Math.log(n); }

    /** Variation (distinctes, -Σ_b Δ Σ n ln n sur X_{-b}) si les personas o deviennent n ; null si o indisponible. */
    static double[] eval(long[] o, long[] n) {
        HashMap<Long, Integer> need = new HashMap<>();
        for (long k : o) need.merge(k, 1, Integer::sum);
        for (Map.Entry<Long, Integer> e : need.entrySet()) if (cnt.getOrDefault(e.getKey(), 0) < e.getValue()) return null;
        HashMap<Long, Integer> d = new HashMap<>();
        for (long k : o) d.merge(k, -1, Integer::sum);
        for (long k : n) d.merge(k, 1, Integer::sum);
        int dd = 0;
        for (Map.Entry<Long, Integer> e : d.entrySet()) {
            int c = cnt.getOrDefault(e.getKey(), 0), c2 = c + e.getValue();
            dd += (c2 > 0 ? 1 : 0) - (c > 0 ? 1 : 0);
        }
        double nl = 0;
        for (int b = 0; b < K; b++) {
            HashMap<Long, Integer> db = new HashMap<>();
            for (Map.Entry<Long, Integer> e : d.entrySet()) if (e.getValue() != 0) db.merge(masked(e.getKey(), b), e.getValue(), Integer::sum);
            for (Map.Entry<Long, Integer> e : db.entrySet()) {
                int c = loo.get(b).getOrDefault(e.getKey(), 0);
                nl += dl(c + e.getValue()) - dl(c);
            }
        }
        return new double[] { dd, -nl };
    }

    static int[] decode(long k) {
        int[] x = new int[K];
        for (int b = 0; b < K; b++) x[b] = (int) (k / radix[b] % in.dom[b]);
        return x;
    }

    /** Les multiensembles de cellules de o et n coïncident sur toutes les tables. */
    static boolean exact(long[] o, long[] n) {
        for (int t = 0; t < in.T; t++) {
            HashMap<Integer, Integer> h = new HashMap<>();
            for (long k : o) h.merge(in.cellOf(t, decode(k)), 1, Integer::sum);
            for (long k : n) h.merge(in.cellOf(t, decode(k)), -1, Integer::sum);
            for (int v : h.values()) if (v != 0) return false;
        }
        return true;
    }

    static long swapS(long a, long b, int S) {          // a avec les valeurs de b sur S
        long r = a;
        for (int x = 0; x < K; x++) if ((S >> x & 1) != 0) r += (b / radix[x] % in.dom[x] - a / radix[x] % in.dom[x]) * radix[x];
        return r;
    }

    static void diag(int[][] pop, int fmax) {
        long[] full = new long[N];
        cnt = new HashMap<>();
        for (int i = 0; i < N; i++) { full[i] = key(pop[i]); cnt.merge(full[i], 1, Integer::sum); }
        loo = new ArrayList<>();
        for (int b = 0; b < K; b++) {
            HashMap<Long, Integer> h = new HashMap<>();
            for (int i = 0; i < N; i++) h.merge(masked(full[i], b), 1, Integer::sum);
            loo.add(h);
        }
        // personas présentes (une paire de personas suffit : l'effet d'un mouvement ne dépend que des personas)
        long[] ps = cnt.keySet().stream().mapToLong(Long::longValue).sorted().toArray();
        int P = ps.length;
        int[][] px = new int[P][];
        for (int p = 0; p < P; p++) px[p] = decode(ps[p]);
        int[] nb = new int[K];
        for (int t = 0; t < in.T; t++) { int m = 0; for (int b : in.scope[t]) m |= 1 << b; for (int b : in.scope[t]) nb[b] |= m; }
        int[] tmask = new int[in.T];
        for (int t = 0; t < in.T; t++) for (int b : in.scope[t]) tmask[t] |= 1 << b;
        Random rz = new Random(12345);
        long[][] z = new long[in.T][];
        for (int t = 0; t < in.T; t++) { z[t] = new long[in.target[t].length]; for (int c = 0; c < z[t].length; c++) z[t][c] = rz.nextLong(); }
        int[][] diff = new int[P][P];
        for (int p = 0; p < P; p++) for (int q = 0; q < P; q++) { int m = 0; for (int b = 0; b < K; b++) if (px[p][b] != px[q][b]) m |= 1 << b; diff[p][q] = m; }

        long sw2 = 0, sw2Up = 0, sw2Lex = 0;               // échanges valides (référence)
        long pairs4 = 0, moves4 = 0, bad = 0, unavailable = 0, up4 = 0, lex4 = 0, neutral4 = 0;
        double bestLoo = 0; int bestD = Integer.MIN_VALUE;
        HashSet<String> seen = new HashSet<>();            // mouvements (multiensembles anciens -> nouveaux) déjà comptés
        int rest = (1 << (K - 1)) - 1;
        for (int S = 1; S <= rest; S++) {
            int nbS = 0;
            for (int b = 0; b < K; b++) if ((S >> b & 1) != 0) nbS |= nb[b];
            int B = nbS & ~S;
            HashMap<Long, List<int[]>> byH = new HashMap<>();
            for (int p = 0; p < P; p++) for (int q = p + 1; q < P; q++) {
                int d = diff[p][q];
                if ((d & S) == 0 || (d & ~S) == 0) continue;
                int F = d & B;
                if (F == 0) {                               // échange valide
                    long[] o = { ps[p], ps[q] }, n = { swapS(ps[p], ps[q], S), swapS(ps[q], ps[p], S) };
                    double[] e = eval(o, n);
                    if (e == null) continue;
                    sw2++;
                    if (e[0] > 0) sw2Up++;
                    else if (e[0] == 0 && e[1] > 1e-9) sw2Lex++;
                    continue;
                }
                if (Integer.bitCount(F) > fmax) continue;
                long h = 0;
                long ni = swapS(ps[p], ps[q], S), nj = swapS(ps[q], ps[p], S);
                int[] xi = px[p], xj = px[q], yi = decode(ni), yj = decode(nj);
                for (int t = 0; t < in.T; t++) {
                    if ((tmask[t] & F) == 0 || (tmask[t] & S) == 0) continue;
                    h += z[t][in.cellOf(t, yi)] + z[t][in.cellOf(t, yj)] - z[t][in.cellOf(t, xi)] - z[t][in.cellOf(t, xj)];
                }
                if (h == 0) continue;                       // écart nul : échange valide sur un S plus petit, déjà couvert
                pairs4++;
                byH.computeIfAbsent(h, k -> new ArrayList<>()).add(new int[] { p, q });
            }
            for (Map.Entry<Long, List<int[]>> e : byH.entrySet()) {
                if (e.getKey() < 0 && byH.containsKey(-e.getKey())) continue;   // chaque couple (h, -h) une seule fois
                List<int[]> other = byH.get(-e.getKey());
                if (other == null) continue;
                for (int[] a : e.getValue()) for (int[] c : other) {
                    long[] o = { ps[a[0]], ps[a[1]], ps[c[0]], ps[c[1]] };
                    long[] n = { swapS(o[0], o[1], S), swapS(o[1], o[0], S), swapS(o[2], o[3], S), swapS(o[3], o[2], S) };
                    long[] so = o.clone(), sn = n.clone(); Arrays.sort(so); Arrays.sort(sn);
                    if (Arrays.equals(so, sn)) continue;    // simple permutation
                    if (!seen.add(Arrays.toString(so) + Arrays.toString(sn))) continue;
                    if (!exact(o, n)) { bad++; continue; }
                    double[] v = eval(o, n);
                    if (v == null) { unavailable++; continue; }
                    moves4++;
                    if (v[0] > 0) up4++;
                    else if (v[0] == 0 && v[1] > 1e-9) lex4++;
                    else if (v[0] == 0) neutral4++;
                    if (v[0] > bestD || (v[0] == bestD && v[1] > bestLoo)) { bestD = (int) v[0]; bestLoo = v[1]; }
                }
            }
        }
        System.out.printf(Locale.ROOT, "personas=%d | échanges valides : %d (gain distinctes : %d, gain entropie à distinctes égales : %d)%n",
                P, sw2, sw2Up, sw2Lex);
        System.out.printf(Locale.ROOT, "  degré 4 (|F| ≤ %d) : paires cassantes=%d, mouvements exacts distincts=%d (faux positifs=%d, personas indisponibles=%d)%n",
                fmax, pairs4, moves4, bad, unavailable);
        System.out.printf(Locale.ROOT, "  gain distinctes > 0 : %d | distinctes égales et entropie loo en hausse : %d | neutres : %d | meilleur : Δdistinctes=%d, ΔN·ΣH(X_-b)=%.3f%n",
                up4, lex4, neutral4, bestD, bestLoo);
    }
}
