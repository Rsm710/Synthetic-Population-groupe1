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
 * Avec --cross, la seconde paire peut échanger un autre ensemble S' : l'écart d'une paire ne porte que sur les tables
 * qui touchent S et F (les autres sont dans S ∪ (B(S) \ F) ou disjointes de S), c'est donc l'écart complet, et les
 * écarts sont indexés tous S confondus.
 *
 * Usage : java -cp out Deg4Diag <contraintes.csv> <fmax> [--cross] <population.csv>...
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
        boolean cross = args.length > 2 && args[2].equals("--cross");
        for (int f = cross ? 3 : 2; f < args.length; f++) {
            long t0 = System.nanoTime();
            diag(CheckPopulation.read(in, Paths.get(args[f])), fmax, cross);
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
        // variations de comptes, regroupées par clé (au plus o.length + n.length clés)
        int m = o.length + n.length, u = 0;
        long[] key = new long[m];
        int[] dv = new int[m];
        for (int r = 0; r < m; r++) {
            long k = r < o.length ? o[r] : n[r - o.length];
            int s = r < o.length ? -1 : 1, x = 0;
            while (x < u && key[x] != k) x++;
            if (x == u) { key[u] = k; dv[u++] = 0; }
            dv[x] += s;
        }
        int dd = 0;
        for (int x = 0; x < u; x++) {
            int c = cnt.getOrDefault(key[x], 0), c2 = c + dv[x];
            if (c2 < 0) return null;                        // persona ancienne indisponible
            dd += (c2 > 0 ? 1 : 0) - (c > 0 ? 1 : 0);
        }
        double nl = 0;
        long[] mk = new long[u];
        int[] mv = new int[u];
        for (int b = 0; b < K; b++) {
            int w = 0;
            for (int x = 0; x < u; x++) {
                if (dv[x] == 0) continue;
                long k = masked(key[x], b);
                int y = 0;
                while (y < w && mk[y] != k) y++;
                if (y == w) { mk[w] = k; mv[w++] = 0; }
                mv[y] += dv[x];
            }
            for (int y = 0; y < w; y++) {
                if (mv[y] == 0) continue;
                int c = loo.get(b).getOrDefault(mk[y], 0);
                nl += dl(c + mv[y]) - dl(c);
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
        int[][] xo = new int[o.length][], xn = new int[n.length][];
        for (int r = 0; r < o.length; r++) xo[r] = decode(o[r]);
        for (int r = 0; r < n.length; r++) xn[r] = decode(n[r]);
        int[] co = new int[o.length], cn = new int[n.length];
        for (int t = 0; t < in.T; t++) {
            for (int r = 0; r < o.length; r++) co[r] = in.cellOf(t, xo[r]);
            for (int r = 0; r < n.length; r++) cn[r] = in.cellOf(t, xn[r]);
            Arrays.sort(co); Arrays.sort(cn);
            if (!Arrays.equals(co, cn)) return false;
        }
        return true;
    }

    static long swapS(long a, long b, int S) {          // a avec les valeurs de b sur S
        long r = a;
        for (int x = 0; x < K; x++) if ((S >> x & 1) != 0) r += (b / radix[x] % in.dom[x] - a / radix[x] % in.dom[x]) * radix[x];
        return r;
    }

    static int nbOf(int[] nb, int S) { int m = 0; for (int b = 0; b < K; b++) if ((S >> b & 1) != 0) m |= nb[b]; return m; }

    static void diag(int[][] pop, int fmax, boolean cross) {
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

        // Mouvement de paire canonique (p, q, A) : échanger S ne dépend que de A = S ∩ d (d = attributs où p, q diffèrent),
        // et A, d \ A donnent le même mouvement : on prend A ∋ plus petit attribut de d. F = d ∩ B(S) ne fait que grandir
        // quand on ajoute à S des attributs hors de d : il est minimal pour S = A (ou d \ A).
        long sw2 = 0, sw2Up = 0, sw2Lex = 0;               // échanges valides (référence)
        long pairs4 = 0, moves4 = 0, bad = 0, unavailable = 0, up4 = 0, lex4 = 0, neutral4 = 0;
        double bestLoo = 0; int bestD = Integer.MIN_VALUE;
        HashSet<Long> seen = new HashSet<>();              // hachage 64 bits des mouvements (anciens -> nouveaux) déjà comptés
        long[] gh = new long[1 << 16], ginfo = new long[1 << 16];   // écarts (h, A<<20 | p<<10 | q)
        int gn = 0;
        for (int p = 0; p < P; p++) for (int q = p + 1; q < P; q++) {
            int d = diff[p][q], low = d & -d, rest = d & ~low;
            for (int sub = rest; ; sub = (sub - 1) & rest) {   // A = low | sub, A ≠ d
                int A = low | sub;
                if (A != d) {
                    int A2 = d & ~A;
                    int F = Math.min(Integer.bitCount(d & nbOf(nb, A) & ~A), Integer.bitCount(d & nbOf(nb, A2) & ~A2));
                    long ni = swapS(ps[p], ps[q], A), nj = swapS(ps[q], ps[p], A);
                    if (F == 0) {                           // échange valide
                        double[] e = eval(new long[] { ps[p], ps[q] }, new long[] { ni, nj });
                        if (e != null) {
                            sw2++;
                            if (e[0] > 0) sw2Up++;
                            else if (e[0] == 0 && e[1] > 1e-9) sw2Lex++;
                        }
                    } else if (F <= fmax) {
                        // seules les tables qui touchent A et d \ A changent : h est l'écart complet
                        long h = 0;
                        int[] xi = px[p], xj = px[q], yi = decode(ni), yj = decode(nj);
                        for (int t = 0; t < in.T; t++) {
                            if ((tmask[t] & A) == 0 || (tmask[t] & A2) == 0) continue;
                            h += z[t][in.cellOf(t, yi)] + z[t][in.cellOf(t, yj)] - z[t][in.cellOf(t, xi)] - z[t][in.cellOf(t, xj)];
                        }
                        if (h != 0) {
                            pairs4++;
                            if (gn == gh.length) { gh = Arrays.copyOf(gh, 2 * gn); ginfo = Arrays.copyOf(ginfo, 2 * gn); }
                            gh[gn] = h; ginfo[gn++] = (long) A << 20 | p << 10 | q;
                        }
                    }
                }
                if (sub == 0) break;
            }
        }
        // écarts h dont l'opposé existe, puis regroupement des seules paires concernées
        long[] sorted = Arrays.copyOf(gh, gn);
        Arrays.sort(sorted);
        HashMap<Long, List<long[]>> byH = new HashMap<>();
        for (int g = 0; g < gn; g++) {
            long h = gh[g];
            if (h == Long.MIN_VALUE || Arrays.binarySearch(sorted, -h) < 0) continue;
            byH.computeIfAbsent(h, k -> new ArrayList<>()).add(new long[] { ginfo[g] >> 20, ginfo[g] >> 10 & 1023, ginfo[g] & 1023 });
        }
        long couples = 0, sameUp = 0, sameLex = 0, sameMoves = 0;
        for (Map.Entry<Long, List<long[]>> e : byH.entrySet()) if (e.getKey() > 0) couples += (long) e.getValue().size() * byH.get(-e.getKey()).size();
        System.out.printf(Locale.ROOT, "  %d écarts ayant un opposé, %d couples de paires à tester%n", byH.size(), couples);
        for (Map.Entry<Long, List<long[]>> e : byH.entrySet()) {
            if (e.getKey() < 0) continue;
            for (long[] a : e.getValue()) for (long[] c : byH.get(-e.getKey())) {
                int A1 = (int) a[0], A2 = (int) c[0];
                int p1 = (int) a[1], q1 = (int) a[2], p2 = (int) c[1], q2 = (int) c[2];
                // un même S existe si A1 (ou d1 \ A1) et A2 (ou d2 \ A2) coïncident sur d1 ∩ d2
                int d1 = diff[p1][q1], d2 = diff[p2][q2], m = d1 & d2;
                boolean same = (A1 & m) == (A2 & m) || (A1 & m) == (~A2 & m);
                if (!same && !cross) continue;
                long[] o = { ps[p1], ps[q1], ps[p2], ps[q2] };
                long[] n = { swapS(o[0], o[1], A1), swapS(o[1], o[0], A1), swapS(o[2], o[3], A2), swapS(o[3], o[2], A2) };
                long[] so = o.clone(), sn = n.clone(); Arrays.sort(so); Arrays.sort(sn);
                if (Arrays.equals(so, sn)) continue;        // simple permutation
                long hk = 0;
                for (long k : so) hk = hk * 0x9E3779B97F4A7C15L + k;
                for (long k : sn) hk = hk * 0x9E3779B97F4A7C15L + k;
                if (!seen.add(hk)) continue;
                if (!exact(o, n)) { bad++; continue; }
                double[] v = eval(o, n);
                if (v == null) { unavailable++; continue; }
                moves4++;
                if (same) sameMoves++;
                if (v[0] > 0) { up4++; if (same) sameUp++; }
                else if (v[0] == 0 && v[1] > 1e-9) { lex4++; if (same) sameLex++; }
                else if (v[0] == 0) neutral4++;
                if (v[0] > bestD || (v[0] == bestD && v[1] > bestLoo)) { bestD = (int) v[0]; bestLoo = v[1]; }
            }
        }
        System.out.printf(Locale.ROOT, "personas=%d | échanges valides : %d (gain distinctes : %d, gain entropie à distinctes égales : %d)%n",
                P, sw2, sw2Up, sw2Lex);
        System.out.printf(Locale.ROOT, "  degré 4 (|F| ≤ %d) : paires cassantes=%d, mouvements exacts distincts=%d (faux positifs=%d, personas indisponibles=%d)%n",
                fmax, pairs4, moves4, bad, unavailable);
        System.out.printf(Locale.ROOT, "  gain distinctes > 0 : %d | distinctes égales et entropie loo en hausse : %d | neutres : %d | meilleur : Δdistinctes=%d, ΔN·ΣH(X_-b)=%.3f%n",
                up4, lex4, neutral4, bestD, bestLoo);
        if (cross) System.out.printf(Locale.ROOT, "  dont même S : %d mouvements, gain distinctes %d, entropie en hausse %d | S' ≠ S : gain distinctes %d, entropie en hausse %d%n",
                sameMoves, sameUp, sameLex, up4 - sameUp, lex4 - sameLex);
    }
}
