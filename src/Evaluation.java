import java.util.*;

/** Vérification indépendante d'une population, métriques de diversité et bornes. */
public final class Evaluation {

    public static final class Report {
        public long error; public int violatedCells;
        public int size, invalidRows;   // taille de la population, lignes mal formées (longueur ou valeur hors domaine)
        public int distinct, maxMult, duplicatedIndividuals;
        public int extraCopies() { return size - distinct; }   // individus en trop par rapport à une population sans doublon
        public double entropy, entropyNorm;
        @Override public String toString() {
            return String.format(Locale.ROOT,
                "taille=%d (lignes invalides=%d) | erreur L1=%d (cellules violées=%d) | personas distinctes=%d | individus en doublon=%d | copies en trop=%d | "
                + "multiplicité max=%d | entropie=%.4f nats (%.2f%% du max ln N)",
                size, invalidRows, error, violatedCells, distinct, duplicatedIndividuals, extraCopies(), maxMult, entropy, 100 * entropyNorm);
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
            H -= (m / n) * Math.log(m / n);
        }
        r.entropy = H;
        r.entropyNorm = H / Math.log(n);
        return r;
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
