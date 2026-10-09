import java.nio.file.*;
import java.util.*;

/**
 * Bornes rigoureuses par propagation sur les personas (énumération de l'espace, si sa taille le permet).
 *
 * Variables : n_x = nombre de copies de la persona x, avec lo_x <= n_x <= hi_x (au départ 0 et u_x = min des cibles).
 * Chaque cellule c impose  somme_{x dans c} n_x = cible_c, d'où, jusqu'au point fixe :
 *     lo_x >= cible_c - somme_{y dans c, y != x} hi_y
 *     hi_x <= cible_c - somme_{y dans c, y != x} lo_y
 * Copies imposées : une persona avec lo_x >= 2 a au moins lo_x - 1 copies en trop, dans toute solution exacte.
 * Borne sup. sur les distinctes :
 *     distinctes <= min_t somme_c min(#{x dans c : hi_x > 0}, cible_c - somme_{x dans c} max(lo_x - 1, 0))
 *
 * Usage : java -cp out Bounds <contraintes.csv> [population.csv]
 */
public final class Bounds {
    final Instance in;
    public int[][] personas;   // personas actives (hi > 0 au départ)
    public int[] lo, hi;
    int[][] cellOf;            // cellOf[k][t]
    int[][] sumLo, sumHi;      // par table et cellule
    public int rounds;

    public Bounds(Instance in) {
        this.in = in;
        int K = in.K;
        double space = 1;
        for (int d : in.dom) space *= d;
        if (space > 50_000_000L) throw new IllegalArgumentException("espace des personas trop grand : " + (long) space);
        List<int[]> ps = new ArrayList<>();
        List<Integer> us = new ArrayList<>();
        int[] x = new int[K];
        for (long id = 0; id < (long) space; id++) {
            long rem = id;
            for (int a = K - 1; a >= 0; a--) { x[a] = (int) (rem % in.dom[a]); rem /= in.dom[a]; }
            int u = in.N;
            for (int t = 0; t < in.T && u > 0; t++) u = Math.min(u, in.target[t][in.cellOf(t, x)]);
            if (u > 0) { ps.add(x.clone()); us.add(u); }
        }
        personas = ps.toArray(new int[0][]);
        int P = personas.length;
        lo = new int[P]; hi = new int[P];
        cellOf = new int[P][in.T];
        for (int k = 0; k < P; k++) {
            hi[k] = us.get(k);
            for (int t = 0; t < in.T; t++) cellOf[k][t] = in.cellOf(t, personas[k]);
        }
    }

    /** Propagation jusqu'au point fixe ; renvoie false si une contradiction apparaît (instance infaisable). */
    public boolean propagate() {
        int P = personas.length;
        boolean changed = true;
        rounds = 0;
        while (changed) {
            changed = false; rounds++;
            sums();
            for (int k = 0; k < P; k++) {
                for (int t = 0; t < in.T; t++) {
                    int c = cellOf[k][t], tg = in.target[t][c];
                    int l = tg - (sumHi[t][c] - hi[k]), h = tg - (sumLo[t][c] - lo[k]);
                    if (l > lo[k]) { for (int s = 0; s < in.T; s++) sumLo[s][cellOf[k][s]] += l - lo[k]; lo[k] = l; changed = true; }
                    if (h < hi[k]) { for (int s = 0; s < in.T; s++) sumHi[s][cellOf[k][s]] -= hi[k] - h; hi[k] = h; changed = true; }
                    if (lo[k] > hi[k]) return false;
                }
            }
        }
        return true;
    }

    private void sums() {
        sumLo = new int[in.T][]; sumHi = new int[in.T][];
        for (int t = 0; t < in.T; t++) { sumLo[t] = new int[in.target[t].length]; sumHi[t] = new int[in.target[t].length]; }
        for (int k = 0; k < personas.length; k++)
            for (int t = 0; t < in.T; t++) { sumLo[t][cellOf[k][t]] += lo[k]; sumHi[t][cellOf[k][t]] += hi[k]; }
    }

    public int forcedExtra() {
        int s = 0;
        for (int l : lo) s += Math.max(l - 1, 0);
        return s;
    }

    public int distinctUpperBound() {
        int ub = in.N;
        for (int t = 0; t < in.T; t++) {
            int[] act = new int[in.target[t].length], forced = new int[in.target[t].length];
            for (int k = 0; k < personas.length; k++) {
                int c = cellOf[k][t];
                if (hi[k] > 0) act[c]++;
                forced[c] += Math.max(lo[k] - 1, 0);
            }
            int s = 0;
            for (int c = 0; c < act.length; c++) s += Math.min(act[c], in.target[t][c] - forced[c]);
            ub = Math.min(ub, s);
        }
        return ub;
    }

    public static void main(String[] args) throws Exception {
        Instance in = Instance.load(Paths.get(args[0]));
        Bounds b = new Bounds(in);
        int P = b.personas.length;
        System.out.printf("personas actives au départ : %d%n", P);
        if (!b.propagate()) { System.out.println("contradiction : instance infaisable"); return; }
        int alive = 0, fixed = 0, withLo = 0;
        for (int k = 0; k < P; k++) { if (b.hi[k] > 0) alive++; if (b.lo[k] == b.hi[k]) fixed++; if (b.lo[k] > 0) withLo++; }
        System.out.printf(Locale.ROOT, "après propagation (%d tours) : %d personas encore possibles, %d présentes de force (lo > 0), %d fixées (lo = hi)%n",
                b.rounds, alive, withLo, fixed);
        System.out.printf("copies en trop imposées : >= %d  =>  distinctes <= %d (borne par cellules : %d)%n",
                b.forcedExtra(), in.N - b.forcedExtra(), b.distinctUpperBound());

        Integer[] idx = new Integer[P];
        for (int k = 0; k < P; k++) idx[k] = k;
        Arrays.sort(idx, (x, y) -> b.lo[y] - b.lo[x]);
        System.out.println("personas les plus imposées (lo..hi) :");
        for (int q = 0; q < Math.min(10, P) && b.lo[idx[q]] > 1; q++)
            System.out.printf("  %s  %d..%d%n", Arrays.toString(b.personas[idx[q]]), b.lo[idx[q]], b.hi[idx[q]]);

        if (args.length > 1) {   // comparer une population aux bornes
            int[][] pop = CheckPopulation.read(in, Paths.get(args[1]));
            Map<String, Integer> h = new HashMap<>();
            for (int[] x : pop) h.merge(Arrays.toString(x), 1, Integer::sum);
            Map<String, Integer> pos = new HashMap<>();
            for (int k = 0; k < P; k++) pos.put(Arrays.toString(b.personas[k]), k);
            List<Map.Entry<String, Integer>> l = new ArrayList<>(h.entrySet());
            l.sort((x, y) -> y.getValue() - x.getValue());
            int avoidable = 0;
            for (Map.Entry<String, Integer> e : l) {
                int k = pos.get(e.getKey());
                avoidable += e.getValue() - Math.max(b.lo[k], 1);
            }
            System.out.printf("%s : %d distinctes, %d copies en trop dont au plus %d évitables d'après lo%n",
                    args[1], h.size(), pop.length - h.size(), avoidable);
            System.out.println("plus gros paquets (copies, lo..hi) :");
            for (int q = 0; q < Math.min(10, l.size()); q++) {
                int k = pos.get(l.get(q).getKey());
                System.out.printf("  %s  x%d  %d..%d%n", l.get(q).getKey(), l.get(q).getValue(), b.lo[k], b.hi[k]);
            }
        }
    }
}
