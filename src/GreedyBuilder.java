import java.util.*;

/**
 * Plan A : construction gloutonne randomisée, individu par individu (DFS bornée par individu).
 *
 * Pour chaque individu, on fixe les attributs un par un. Le poids d'une valeur v pour
 * l'attribut a, sachant l'affectation partielle P, est :
 *     w(v) = min sur les tables t contenant a de  R_t(P ∪ {a=v})
 * où R_t est la capacité résiduelle (target - count) sommée sur les cellules de t
 * compatibles avec l'affectation partielle. C'est une version dynamique de la borne u_x.
 * L'attribut choisi à chaque étape est celui qui a le moins de valeurs possibles (fail-first).
 * Au dernier attribut, les valeurs qui recréent une persona déjà présente sont pénalisées.
 */
public final class GreedyBuilder {
    final Instance in;
    final Random rnd;

    public GreedyBuilder(Instance in, Random rnd) { this.in = in; this.rnd = rnd; }

    public int nodeBudget = 20000;     // noeuds max de la recherche par individu
    long deadEnds = 0;

    public Population build() {
        Population p = new Population(in);
        int K = in.K;
        int[] cur = new int[K];
        boolean[] assigned = new boolean[K];
        for (int i = 0; i < in.N; i++) {
            // 1) persona nouvelle et compatible ; 2) doublon compatible ; 3) glouton best-effort
            boolean ok = false;
            for (int pass = 0; pass < 2 && !ok; pass++) {
                Arrays.fill(assigned, false);
                nodes = 0;
                ok = dfs(p, cur, assigned, 0, pass == 0);
            }
            if (!ok) { deadEnds++; bestEffort(p, cur, assigned); }
            System.arraycopy(cur, 0, p.ind[i], 0, K);
            p.add(i);
        }
        return p;
    }

    private int nodes;

    /** DFS randomisée : chaque valeur doit laisser une capacité résiduelle > 0 dans toutes les tables. */
    private boolean dfs(Population p, int[] cur, boolean[] assigned, int depth, boolean requireNew) {
        if (depth == in.K) return !requireNew || !p.persona.containsKey(p.key(cur));
        if (++nodes > nodeBudget) return false;
        // fail-first : attribut avec le moins de valeurs à poids > 0
        int a = -1, bestFeas = Integer.MAX_VALUE, ties = 0;
        long[] wa = null;
        for (int b = 0; b < in.K; b++) {
            if (assigned[b]) continue;
            long[] w = new long[in.dom[b]];
            int feas = 0;
            for (int v = 0; v < in.dom[b]; v++) { cur[b] = v; w[v] = weight(p, cur, assigned, b); if (w[v] > 0) feas++; }
            if (feas == 0) return false;
            if (feas < bestFeas) { bestFeas = feas; a = b; wa = w; ties = 1; }
            else if (feas == bestFeas && rnd.nextInt(++ties) == 0) { a = b; wa = w; }
        }
        // valeurs dans un ordre tiré au hasard proportionnellement aux poids
        int d = in.dom[a];
        double[] pw = new double[d];
        for (int v = 0; v < d; v++) pw[v] = Math.max(0, wa[v]);
        assigned[a] = true;
        for (int k = 0; k < d; k++) {
            double tot = 0;
            for (double x : pw) tot += x;
            if (tot <= 0) break;
            double r = rnd.nextDouble() * tot;
            int v = 0;
            while (v < d - 1 && (r -= pw[v]) >= 0) v++;
            if (pw[v] <= 0) { v = 0; while (pw[v] <= 0) v++; }
            pw[v] = 0;
            cur[a] = v;
            if (dfs(p, cur, assigned, depth + 1, requireNew)) return true;
            if (nodes > nodeBudget) break;
        }
        assigned[a] = false;
        return false;
    }

    /**
     * Impasse globale : plus aucune persona n'a de capacité partout. On cherche (branch & bound
     * borné) la persona qui déborde dans le moins de cellules possible.
     */
    private int[] bbBest;
    private int bbBestCost;

    private void bestEffort(Population p, int[] cur, boolean[] assigned) {
        Arrays.fill(assigned, false);
        bbBest = null; bbBestCost = Integer.MAX_VALUE; nodes = 0;
        bb(p, cur, assigned, 0, 0);
        System.arraycopy(bbBest, 0, cur, 0, in.K);
    }

    private void bb(Population p, int[] cur, boolean[] assigned, int a, int cost) {
        if (cost >= bbBestCost) return;
        if (a == in.K) { bbBestCost = cost; bbBest = cur.clone(); return; }
        if (++nodes > nodeBudget && bbBest != null) return;
        int d = in.dom[a];
        int[] add = new int[d];
        Integer[] order = new Integer[d];
        for (int v = 0; v < d; v++) {
            order[v] = v;
            cur[a] = v;
            // tables dont a est le dernier attribut du scope (ordre 0..K-1) : cellule désormais fixée
            for (int t : in.attrTables[a]) {
                boolean complete = true;
                for (int b : in.scope[t]) if (b > a) { complete = false; break; }
                if (complete) {
                    int c = in.cellOf(t, cur);
                    if (in.target[t][c] - p.count[t][c] <= 0) add[v]++;
                }
            }
        }
        Collections.shuffle(Arrays.asList(order), rnd);
        Arrays.sort(order, Comparator.comparingInt(v -> add[v]));
        assigned[a] = true;
        for (int v : order) { cur[a] = v; bb(p, cur, assigned, a + 1, cost + add[v]); }
        assigned[a] = false;
    }

    /** Poids = min des capacités résiduelles (signées) sur les tables contenant a. */
    private long weight(Population p, int[] cur, boolean[] assigned, int a) {
        long w = Long.MAX_VALUE;
        for (int t : in.attrTables[a]) {
            long r = residualSlice(p, t, cur, assigned, a);
            if (r < w) w = r;
        }
        return w == Long.MAX_VALUE ? 1 : w;
    }

    /** Somme de (target - count) sur les cellules de t compatibles avec l'affectation partielle + a. */
    private long residualSlice(Population p, int t, int[] cur, boolean[] assigned, int a) {
        int[] sc = in.scope[t], st = in.stride[t];
        int[] tg = in.target[t], cn = p.count[t];
        long s = 0;
        outer:
        for (int c = 0; c < tg.length; c++) {
            for (int j = 0; j < sc.length; j++) {
                int b = sc[j];
                if (b == a || assigned[b]) {
                    int v = (c / st[j]) % in.dom[b];
                    if (v != cur[b]) continue outer;
                }
            }
            s += tg[c] - cn[c];
        }
        return s;
    }

}
