import java.util.*;

/** Vérification indépendante d'une population, métriques de diversité et bornes. */
public final class Evaluation {

    public static final class Report {
        public long error; public int violatedCells;
        public int size, invalidRows;   // taille de la population, lignes mal formées (longueur ou valeur hors domaine)
        public int distinct, maxMult, duplicatedIndividuals;
        public int extraCopies() { return size - distinct; }   // individus en trop par rapport à une population sans doublon
        public double entropy, entropyNorm;
        // Entropies de Shannon (nats) sur les projections « sans un attribut » X_{-b} :
        public double entropyLoo;       // Σ_b H(X_{-b})
        public double nearDup;          // Σ_b H(X_b | X_{-b}) = K·H(X) - Σ_b H(X_{-b}) ; 0 = aucune paire à distance 1
        public double looDeficit;       // K·ln N - Σ_b H(X_{-b}) ; 0 ssi aucune paire à distance ≤ 1 (ni clone, ni quasi-clone)
        public long pairsD1;            // paires d'individus à distance de Hamming exactement 1
        public long pairsD0;            // paires de clones (distance 0) : Σ C(m, 2)
        public long pairsD2;            // paires d'individus à distance de Hamming exactement 2
        public double looDeficit2;      // C(K,2)·ln N - Σ_{|B|=2} H(X_{-B}) ; 0 ssi aucune paire à distance ≤ 2
        public double freeMI;           // Σ I(X_a ; X_b) sur les paires (a, b) couvertes par aucune table
        @Override public String toString() {
            return String.format(Locale.ROOT,
                "taille=%d (lignes invalides=%d) | erreur L1=%d (cellules violées=%d) | personas distinctes=%d | individus en doublon=%d | copies en trop=%d | "
                + "multiplicité max=%d | entropie=%.4f nats (%.2f%% du max ln N) | Σ H(X_-b)=%.4f (déficit %.4f) | Σ H(X_b|X_-b)=%.4f | déficit rayon 2=%.4f | paires à distance 0 / 1 / 2=%d / %d / %d | "
                + "Σ I paires libres=%.4f",
                size, invalidRows, error, violatedCells, distinct, duplicatedIndividuals, extraCopies(), maxMult, entropy, 100 * entropyNorm,
                entropyLoo, looDeficit, nearDup, looDeficit2, pairsD0, pairsD1, pairsD2, freeMI);
        }
    }

    /** Recalcule tout à partir de zéro, sans réutiliser l'état incrémental. */
    public static Report evaluate(Instance in, int[][] pop) {
        Report r = new Report();
        r.size = pop.length;
        for (int[] x : pop) {
            boolean ok = x.length == in.K;
            for (int a = 0; ok && a < in.K; a++) if (x[a] < 0 || x[a] >= in.dom[a]) ok = false;
            if (!ok) r.invalidRows++;
        }
        if (r.invalidRows > 0) { r.error = Long.MAX_VALUE; return r; }
        for (int t = 0; t < in.T; t++) {
            int[] c = new int[in.target[t].length];
            for (int[] x : pop) c[in.cellOf(t, x)]++;
            for (int k = 0; k < c.length; k++) {
                int d = Math.abs(c[k] - in.target[t][k]);
                r.error += d;
                if (d != 0) r.violatedCells++;
            }
        }
        Map<String, Integer> h = new HashMap<>();
        for (int[] x : pop) h.merge(Arrays.toString(x), 1, Integer::sum);
        r.distinct = h.size();
        double n = pop.length, H = 0;
        for (int m : h.values()) {
            r.maxMult = Math.max(r.maxMult, m);
            if (m > 1) r.duplicatedIndividuals += m;
            r.pairsD0 += (long) m * (m - 1) / 2;
            H -= (m / n) * Math.log(m / n);
        }
        r.entropy = H;
        r.entropyNorm = H / Math.log(n);
        // X_{-b} : on regroupe les individus qui coïncident hors de b. Dans un groupe de taille m dont les valeurs
        // de b ont les comptes m_v, les paires à distance 1 (différentes seulement sur b) sont C(m,2) - Σ_v C(m_v,2).
        for (int b = 0; b < in.K; b++) {
            Map<String, int[]> g = new HashMap<>();
            int db = in.dom[b];
            for (int[] x : pop) {
                int[] y = x.clone(); y[b] = -1;
                g.computeIfAbsent(Arrays.toString(y), k -> new int[db])[x[b]]++;
            }
            for (int[] cv : g.values()) {
                long m = 0, same = 0;
                for (int c : cv) { m += c; same += (long) c * (c - 1) / 2; }
                r.entropyLoo -= (m / n) * Math.log(m / n);
                r.pairsD1 += m * (m - 1) / 2 - same;
            }
        }
        r.nearDup = in.K * H - r.entropyLoo;
        r.looDeficit = in.K * Math.log(n) - r.entropyLoo;
        // rayon 2 : X_{-B} pour les C(K,2) paires B = {b, c}
        for (int b = 0; b < in.K; b++)
            for (int c = b + 1; c < in.K; c++) {
                Map<String, Integer> g = new HashMap<>();
                for (int[] x : pop) { int[] y = x.clone(); y[b] = -1; y[c] = -1; g.merge(Arrays.toString(y), 1, Integer::sum); }
                double hb = 0;
                for (int m : g.values()) hb -= (m / n) * Math.log(m / n);
                r.looDeficit2 += Math.log(n) - hb;
            }
        for (int i = 0; i < pop.length; i++)
            for (int j = i + 1; j < pop.length; j++) {
                int d = 0;
                for (int a = 0; a < in.K && d <= 2; a++) if (pop[i][a] != pop[j][a]) d++;
                if (d == 2) r.pairsD2++;
            }
        // information mutuelle des paires libres : H(X_a) + H(X_b) - H(X_a, X_b)
        for (int[] ab : freePairs(in)) {
            int a = ab[0], b = ab[1];
            double[][] c = new double[in.dom[a]][in.dom[b]];
            double[] ca = new double[in.dom[a]], cb = new double[in.dom[b]];
            for (int[] x : pop) { c[x[a]][x[b]]++; ca[x[a]]++; cb[x[b]]++; }
            for (int u = 0; u < in.dom[a]; u++)
                for (int v = 0; v < in.dom[b]; v++)
                    if (c[u][v] > 0) r.freeMI += (c[u][v] / n) * Math.log(c[u][v] * n / (ca[u] * cb[v]));
        }
        return r;
    }

    /** Paires d'attributs (a < b) qu'aucune table ne contient ensemble : leur loi jointe n'est pas imposée. */
    public static List<int[]> freePairs(Instance in) {
        boolean[][] cov = new boolean[in.K][in.K];
        for (int t = 0; t < in.T; t++) for (int a : in.scope[t]) for (int b : in.scope[t]) cov[a][b] = true;
        List<int[]> l = new ArrayList<>();
        for (int a = 0; a < in.K; a++) for (int b = a + 1; b < in.K; b++) if (!cov[a][b]) l.add(new int[]{a, b});
        return l;
    }

    /**
     * Bornes du cours (u_x = min des targets des cellules couvrant x), par énumération
     * de l'espace des personas, si sa taille le permet. Donne aussi une borne supérieure
     * sur le nombre de personas distinctes :
     *     distinct <= min_t  somme_c min(target_c, nb de personas actives dans c)
     */
    public static void printBounds(Instance in, long maxSpace) {
        double space = 1;
        for (int d : in.dom) space *= d;
        if (space > maxSpace) { System.out.println("Espace des personas trop grand pour l'énumération ("
                + (long) space + "), bornes ignorées."); return; }
        int K = in.K;
        int[] x = new int[K];
        long active = 0, sumU = 0; int maxU = 0;
        int[][] activeInCell = new int[in.T][];
        for (int t = 0; t < in.T; t++) activeInCell[t] = new int[in.target[t].length];
        int[] cells = new int[in.T];
        for (long id = 0; id < (long) space; id++) {
            long rem = id;
            for (int a = K - 1; a >= 0; a--) { x[a] = (int) (rem % in.dom[a]); rem /= in.dom[a]; }
            int u = Integer.MAX_VALUE;
            for (int t = 0; t < in.T && u > 0; t++) {
                cells[t] = in.cellOf(t, x);
                u = Math.min(u, in.target[t][cells[t]]);
            }
            if (u == Integer.MAX_VALUE) u = in.N;
            if (u > 0) {
                active++; sumU += u; maxU = Math.max(maxU, u);
                for (int t = 0; t < in.T; t++) activeInCell[t][cells[t]]++;
            }
        }
        int ub = in.N;
        for (int t = 0; t < in.T; t++) {
            int s = 0;
            for (int c = 0; c < in.target[t].length; c++) s += Math.min(in.target[t][c], activeInCell[t][c]);
            ub = Math.min(ub, s);
        }
        System.out.printf(Locale.ROOT,
            "Bornes : %d personas, %d actives (%.1f%% éliminées), somme u_x=%d, u_x max=%d, "
            + "borne sup. sur les distinctes=%d%n",
            (long) space, active, 100.0 * (space - active) / space, sumU, maxU, ub);
    }
}
