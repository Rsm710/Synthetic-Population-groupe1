import java.util.*;

/**
 * Variante de ColumnGreedy : « places uniques d'abord, doublons concentrés ».
 *
 * Seule la préférence entre valeurs faisables change (ordre, look-ahead, réparation, retours arrière identiques).
 * Une branche = persona partielle (attributs connus, a compris). Son potentiel pot = nombre de personas complètes
 * actives (u_x > 0) qui la prolongent : la branche ne peut pas contenir plus de pot personas distinctes.
 *  - s'il reste une branche avec des places libres (pot - jumeaux déjà envoyés > 0), on y va, en séparant
 *    les jumeaux comme ColumnGreedy (mr / (1 + tw)^beta) ;
 *  - sinon l'individu sera forcément un doublon : on le met avec le plus de jumeaux déjà présents
 *    (concentration des doublons, comme à l'optimum CP-SAT).
 * Sans énumération possible de l'espace des personas, pot = produit des domaines non encore fixés.
 */
public final class ColumnGreedy2 extends ColumnGreedy {
    public double gamma = 1.0;          // force de la concentration des doublons
    final int[][] active;               // personas actives, si l'espace est énumérable
    private HashMap<Long, Integer> pot; // branche -> potentiel, pour la colonne en cours
    private boolean[] known;

    public ColumnGreedy2(Instance in, Random rnd) {
        super(in, rnd);
        double space = 1;
        for (int d : in.dom) space *= d;
        List<int[]> ps = new ArrayList<>();
        if (space <= 5_000_000) {
            int[] x = new int[in.K];
            for (long id = 0; id < (long) space; id++) {
                long rem = id;
                for (int b = in.K - 1; b >= 0; b--) { x[b] = (int) (rem % in.dom[b]); rem /= in.dom[b]; }
                boolean ok = true;
                for (int t = 0; t < in.T && ok; t++) if (in.target[t][in.cellOf(t, x)] == 0) ok = false;
                if (ok) ps.add(x.clone());
            }
            active = ps.toArray(new int[0][]);
        } else active = null;
    }

    @Override protected void prepareColumn(boolean[] known, int a) {
        this.known = known;
        if (active == null) return;
        pot = new HashMap<>();
        for (int[] x : active) pot.merge(knownKey(x, known), 1, Integer::sum);
    }

    private int potential(long branch) {
        if (active != null) return pot.getOrDefault(branch, 0);
        long p = 1;
        for (int b = 0; b < in.K; b++) if (!known[b]) p = Math.min(p * in.dom[b], in.N);
        return (int) p;
    }

    @Override protected double diversity(int mr, int[] tw, int v, long branch) {
        int free = potential(branch) - tw[v];
        double noise = 1 + 0.1 * rnd.nextDouble();
        if (free > 0) return 5e4 + mr / Math.pow(1 + tw[v], beta) * noise;   // place libre : séparer comme ColumnGreedy
        return Math.min(4e4, mr * Math.pow(1 + tw[v], gamma)) * noise;
    }
}
