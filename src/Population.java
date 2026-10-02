import java.util.*;

/**
 * État courant : N individus, compteurs de chaque cellule, histogramme des personas.
 * Les deux coûts sont maintenus de façon incrémentale :
 *   error   = somme |count - target| sur toutes les cellules (0 = toutes les marges exactes)
 *   dupCost = somme n_x ln n_x sur les personas (0 = aucun doublon ; minimiser <=> maximiser l'entropie)
 */
public final class Population {
    final Instance in;
    final int[][] ind;          // ind[i][a]
    final int[][] count;        // count[t][cellule]
    final int[] tableErr;       // erreur par table
    final HashMap<Long, Integer> persona = new HashMap<>();
    final long[] radix;         // codage des personas en long
    final double[] flog;        // flog[n] = n ln n
    long error = 0;
    double dupCost = 0;
    int distinct = 0;

    public Population(Instance in) {
        this.in = in;
        ind = new int[in.N][in.K];
        count = new int[in.T][];
        for (int t = 0; t < in.T; t++) count[t] = new int[in.target[t].length];
        tableErr = new int[in.T];
        radix = new long[in.K];
        double bits = 0;
        long r = 1;
        for (int a = in.K - 1; a >= 0; a--) { radix[a] = r; r *= in.dom[a]; bits += Math.log(in.dom[a]) / Math.log(2); }
        if (bits > 62) throw new IllegalArgumentException("Espace des personas trop grand pour un codage long");
        flog = new double[in.N + 2];
        for (int n = 1; n < flog.length; n++) flog[n] = n * Math.log(n);
        // état vide : toutes les cellules sont à count = 0
        for (int t = 0; t < in.T; t++)
            for (int c = 0; c < count[t].length; c++) { tableErr[t] += in.target[t][c]; error += in.target[t][c]; }
    }

    public long key(int[] v) {
        long k = 0;
        for (int a = 0; a < in.K; a++) k += v[a] * radix[a];
        return k;
    }

    private void bump(int t, int c, int d) {
        int tg = in.target[t][c], old = count[t][c], nw = old + d;
        int de = Math.abs(nw - tg) - Math.abs(old - tg);
        count[t][c] = nw;
        tableErr[t] += de;
        error += de;
    }

    private void bumpPersona(long k, int d) {
        int old = persona.getOrDefault(k, 0), nw = old + d;
        if (nw == 0) persona.remove(k); else persona.put(k, nw);
        dupCost += flog[nw] - flog[old];
        if (old == 0) distinct++;
        if (nw == 0) distinct--;
    }

    /** Ajoute l'individu i (ses valeurs doivent être déjà écrites dans ind[i]). */
    public void add(int i) {
        for (int t = 0; t < in.T; t++) bump(t, in.cellOf(t, ind[i]), +1);
        bumpPersona(key(ind[i]), +1);
    }

    // ----- mouvements élémentaires, réversibles -----

    /** Change ind[i][a] en v (mise à jour des seules tables contenant a). */
    public void set(int i, int a, int v) {
        int old = ind[i][a];
        if (old == v) return;
        bumpPersona(key(ind[i]), -1);
        for (int t : in.attrTables[a]) bump(t, in.cellOf(t, ind[i]), -1);
        ind[i][a] = v;
        for (int t : in.attrTables[a]) bump(t, in.cellOf(t, ind[i]), +1);
        bumpPersona(key(ind[i]), +1);
    }

    /** Échange les valeurs des attributs attrs[0..m) entre i et j. */
    public void swap(int i, int j, int[] attrs, int m) {
        bumpPersona(key(ind[i]), -1);
        bumpPersona(key(ind[j]), -1);
        // tables touchées : contenant au moins un attribut échangé
        boolean[] touched = touchedTables(attrs, m);
        for (int t = 0; t < in.T; t++) if (touched[t]) {
            bump(t, in.cellOf(t, ind[i]), -1);
            bump(t, in.cellOf(t, ind[j]), -1);
        }
        for (int q = 0; q < m; q++) {
            int a = attrs[q], tmp = ind[i][a];
            ind[i][a] = ind[j][a];
            ind[j][a] = tmp;
        }
        for (int t = 0; t < in.T; t++) if (touched[t]) {
            bump(t, in.cellOf(t, ind[i]), +1);
            bump(t, in.cellOf(t, ind[j]), +1);
        }
        bumpPersona(key(ind[i]), +1);
        bumpPersona(key(ind[j]), +1);
    }

    private boolean[] touchedTables(int[] attrs, int m) {
        boolean[] b = new boolean[in.T];
        for (int q = 0; q < m; q++) for (int t : in.attrTables[attrs[q]]) b[t] = true;
        return b;
    }

    public int personaCount(int i) { return persona.getOrDefault(key(ind[i]), 0); }

    public int[][] snapshot() {
        int[][] s = new int[in.N][];
        for (int i = 0; i < in.N; i++) s[i] = ind[i].clone();
        return s;
    }
}
