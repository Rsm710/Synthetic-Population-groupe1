import java.util.*;

/**
 * Plan B : recuit simulé sur  F = W * error + dupCost.
 *  - error   : écart L1 aux marges (doit finir à 0)
 *  - dupCost : somme n_x ln n_x (0 si aucun doublon)
 * Mouvements :
 *  - set  : changer un attribut d'un individu
 *  - swap : échanger un attribut entre deux individus (préserve la marge unaire)
 *  - swapScope : échanger tout le scope d'une table entre deux individus (préserve cette table)
 * Les individus à modifier sont choisis en priorité dans les cellules en excès (tant que
 * error > 0), puis parmi les doublons (une fois error = 0).
 */
public final class LocalSearch {
    final Instance in;
    final Random rnd;
    public double W = 3.0, T0 = 2.0, T1 = 0.02;
    public long maxIters = 3_000_000;
    public long timeLimitMs = 10_000;

    public long iters;
    int[][] best;
    long bestError;
    double bestDup;

    public LocalSearch(Instance in, Random rnd) { this.in = in; this.rnd = rnd; }

    public int[][] run(Population p) {
        best = p.snapshot(); bestError = p.error; bestDup = p.dupCost;
        int K = in.K;
        int[] attrs = new int[K];
        long start = System.currentTimeMillis();
        double logRatio = Math.log(T1 / T0);
        double frac = 0, temp = T0;
        for (long it = 0; it < maxIters; it++) {
            if ((it & 1023) == 0) {
                long el = System.currentTimeMillis() - start;
                if (el > timeLimitMs) break;
                // température pilotée par la plus avancée des deux horloges (itérations / temps)
                frac = Math.max(it / (double) maxIters, el / (double) timeLimitMs);
                temp = T0 * Math.exp(logRatio * frac);
            }
            iters = it;
            double before = W * p.error + p.dupCost;

            int i = pickIndividual(p);
            double r = rnd.nextDouble();
            if (r < 0.3) {                                  // set
                int a = rnd.nextInt(K);
                if (in.dom[a] < 2) continue;
                int old = p.ind[i][a];
                int v = rnd.nextInt(in.dom[a] - 1); if (v >= old) v++;
                p.set(i, a, v);
                if (!accept(p.error * W + p.dupCost - before, temp)) p.set(i, a, old);
            } else {                                        // swap / swapScope
                int j = rnd.nextInt(in.N);
                if (j == i) continue;
                int m;
                if (r < 0.65) { attrs[0] = rnd.nextInt(K); m = 1; }
                else {
                    int t = rnd.nextInt(in.T);
                    m = in.scope[t].length;
                    System.arraycopy(in.scope[t], 0, attrs, 0, m);
                }
                p.swap(i, j, attrs, m);
                if (!accept(p.error * W + p.dupCost - before, temp)) p.swap(i, j, attrs, m);
            }

            if (p.error < bestError || (p.error == bestError && p.dupCost < bestDup - 1e-9)) {
                bestError = p.error; bestDup = p.dupCost; best = p.snapshot();
            }
        }
        return best;
    }

    private boolean accept(double delta, double temp) {
        return delta <= 0 || rnd.nextDouble() < Math.exp(-delta / temp);
    }

    /** Individu ciblé : dans une cellule en excès si error > 0, sinon un doublon si possible. */
    private int pickIndividual(Population p) {
        if (p.error > 0 && rnd.nextDouble() < 0.7) {
            int nViol = 0, tSel = -1;
            for (int t = 0; t < in.T; t++) if (p.tableErr[t] > 0 && rnd.nextInt(++nViol) == 0) tSel = t;
            if (tSel >= 0) {
                int nOver = 0, cSel = -1;
                for (int c = 0; c < p.count[tSel].length; c++)
                    if (p.count[tSel][c] > in.target[tSel][c] && rnd.nextInt(++nOver) == 0) cSel = c;
                if (cSel >= 0) {
                    int nIn = 0, iSel = -1;
                    for (int i = 0; i < in.N; i++)
                        if (in.cellOf(tSel, p.ind[i]) == cSel && rnd.nextInt(++nIn) == 0) iSel = i;
                    if (iSel >= 0) return iSel;
                }
            }
        } else if (p.error == 0 && p.distinct < in.N && rnd.nextDouble() < 0.8) {
            for (int tries = 0; tries < 50; tries++) {
                int i = rnd.nextInt(in.N);
                if (p.personaCount(i) > 1) return i;
            }
        }
        return rnd.nextInt(in.N);
    }
}
