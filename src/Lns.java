import java.util.*;

/**
 * Recherche à grand voisinage qui conserve l'exactitude.
 *
 * On choisit un sous-ensemble I d'individus : des copies d'une persona en doublon (tirée selon son nombre
 * de copies en trop) et des individus proches d'elle (distance de Hamming). Les comptes de I dans chaque
 * cellule deviennent les cibles d'un sous-problème de taille |I|, reconstruit de zéro par ColumnGreedy.
 * Si la reconstruction est exacte, remplacer I par elle laisse toutes les marges globales inchangées.
 * On accepte si le nombre de distinctes ne baisse pas et on garde la meilleure population.
 */
public final class Lns {
    final Instance in;
    final Random rnd;
    public long timeMs = 3000;      // budget de temps
    public int size = 60;           // taille du sous-ensemble reconstruit
    public double beta = 2.0;       // diversité du greedy des sous-problèmes
    public long iterations, exactSub, improved, same;   // diagnostics (same : reconstruction = même multiensemble)

    public Lns(Instance in, Random rnd) { this.in = in; this.rnd = rnd; }

    public int run(int[][] pop) {
        int N = pop.length, K = in.K;
        HashMap<List<Integer>, List<Integer>> groups = new HashMap<>();
        int cur = distinct(pop, groups), best = cur;
        int[][] bestPop = copy(pop);
        long end = System.nanoTime() + timeMs * 1_000_000L;
        iterations = 0; exactSub = 0; improved = 0; same = 0;
        while (System.nanoTime() < end && cur < N) {
            iterations++;
            // persona en doublon, tirée proportionnellement à ses copies en trop
            int extra = N - cur, pick = rnd.nextInt(extra);
            List<Integer> g = null;
            for (List<Integer> idx : groups.values()) { pick -= idx.size() - 1; if (pick < 0) { g = idx; break; } }
            int[] p = pop[g.get(0)];
            // sous-ensemble : copies de p (au plus la moitié), puis les plus proches par Hamming (ex aequo au hasard)
            boolean[] in_ = new boolean[N];
            List<Integer> sub = new ArrayList<>();
            List<Integer> copies = new ArrayList<>(g);
            Collections.shuffle(copies, rnd);
            for (int i : copies) { if (sub.size() >= size / 2) break; sub.add(i); in_[i] = true; }
            Integer[] order = new Integer[N];
            double[] dist = new double[N];
            for (int i = 0; i < N; i++) {
                order[i] = i;
                int d = 0;
                for (int a = 0; a < K; a++) if (pop[i][a] != p[a]) d++;
                dist[i] = d + rnd.nextDouble();
            }
            Arrays.sort(order, Comparator.comparingDouble(i -> dist[i]));
            for (int i : order) { if (sub.size() >= size) break; if (!in_[i]) { sub.add(i); in_[i] = true; } }

            int m = sub.size();
            int[][] tg = new int[in.T][];
            for (int t = 0; t < in.T; t++) {
                tg[t] = new int[in.target[t].length];
                for (int i : sub) tg[t][in.cellOf(t, pop[i])]++;
            }
            Instance si = in.withTargets(m, tg);
            ColumnGreedy cg = new ColumnGreedy(si, rnd);
            cg.beta = beta;
            int[][] np = cg.build();
            if (Evaluation.evaluate(si, np).error != 0) continue;
            exactSub++;
            if (sameMultiset(np, sub, pop)) { same++; continue; }
            int[][] old = new int[m][];
            for (int q = 0; q < m; q++) { old[q] = pop[sub.get(q)]; pop[sub.get(q)] = np[q]; }
            int nd = distinct(pop, groups);
            if (nd >= cur) {
                if (nd > cur) improved++;
                cur = nd;
                if (cur > best) { best = cur; bestPop = copy(pop); }
            } else {
                for (int q = 0; q < m; q++) pop[sub.get(q)] = old[q];
                distinct(pop, groups);
            }
        }
        for (int i = 0; i < N; i++) pop[i] = bestPop[i];
        return best;
    }

    private static boolean sameMultiset(int[][] np, List<Integer> sub, int[][] pop) {
        HashMap<List<Integer>, Integer> c = new HashMap<>();
        for (int[] x : np) c.merge(Arrays.stream(x).boxed().toList(), 1, Integer::sum);
        for (int i : sub) if (c.merge(Arrays.stream(pop[i]).boxed().toList(), -1, Integer::sum) < 0) return false;
        return true;
    }

    /** Nombre de distinctes ; remplit groups (persona -> indices de ses copies). */
    private static int distinct(int[][] pop, HashMap<List<Integer>, List<Integer>> groups) {
        groups.clear();
        for (int i = 0; i < pop.length; i++) {
            List<Integer> key = new ArrayList<>(pop[i].length);
            for (int v : pop[i]) key.add(v);
            groups.computeIfAbsent(key, k -> new ArrayList<>()).add(i);
        }
        return groups.size();
    }

    private static int[][] copy(int[][] pop) {
        int[][] c = new int[pop.length][];
        for (int i = 0; i < pop.length; i++) c[i] = pop[i].clone();
        return c;
    }
}
