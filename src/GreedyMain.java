import java.nio.file.*;
import java.util.*;

/**
 * Usage : java -cp out GreedyMain <contraintes.csv> [--runs R] [--seed s]
 *                                 [--fc on|off] [--retries r] [--backtrack b] [--depth d] [--restarts m] [--repair iters] [--perfect bonus] [--cross lambda] [--beta b] [--out population.csv] [--quiet]
 */
public final class GreedyMain {
    public static void main(String[] args) throws Exception {
        Instance in = Instance.load(Paths.get(args[0]));
        int runs = 10; long seed = 0; boolean fc = true, quiet = false;
        double beta = 1.0; int retries = 3, backtrack = 40, repair = 500, depth = 2, restarts = 10; double cross = 0, perfect = 1000;
        String out = "population.csv";
        for (int k = 1; k < args.length; k++) {
            switch (args[k]) {
                case "--runs" -> runs = Integer.parseInt(args[++k]);
                case "--seed" -> seed = Long.parseLong(args[++k]);
                case "--fc" -> fc = args[++k].equals("on");
                case "--beta" -> beta = Double.parseDouble(args[++k]);
                case "--out" -> out = args[++k];
                case "--quiet" -> quiet = true;
                case "--depth" -> depth = Integer.parseInt(args[++k]);
                case "--restarts" -> restarts = Integer.parseInt(args[++k]);
                case "--cross" -> cross = Double.parseDouble(args[++k]);
                case "--perfect" -> perfect = Double.parseDouble(args[++k]);
                case "--repair" -> repair = Integer.parseInt(args[++k]);
                case "--retries" -> retries = Integer.parseInt(args[++k]);
                case "--backtrack" -> backtrack = Integer.parseInt(args[++k]);
                default -> throw new IllegalArgumentException("Option inconnue : " + args[k]);
            }
        }
        System.out.printf("Instance : N=%d, K=%d, %d tables%n", in.N, in.K, in.T);
        // Un run = au plus 1 + restarts constructions complètes (même générateur aléatoire, tirages différents),
        // arrêtées dès qu'une population est exacte selon Evaluation.evaluate.
        int exact = 0, firstTry = 0, sizeOk = 0, sumAttempts = 0, maxAttempts = 0; int[][] best = null; Evaluation.Report bestRep = null;
        double sumMs = 0, maxMs = 0, sumDist = 0, sumH = 0;
        for (int r = 0; r < runs; r++) {
            Random rnd = new Random(seed + r);
            long t1 = System.nanoTime();
            int[][] pop = null; Evaluation.Report rep = null; ColumnGreedy g = null;
            int attempts = 0, deep = 0;
            while (attempts <= restarts) {
                attempts++;
                g = new ColumnGreedy(in, rnd);
                g.beta = beta; g.lookAhead = fc; g.crossPenalty = cross; g.perfectBonus = perfect; g.columnRetries = retries;
                g.backtrackTries = backtrack; g.backtrackDepth = depth; g.repairIters = repair;
                pop = g.build();
                rep = Evaluation.evaluate(in, pop);
                deep += g.deepBacktracks;
                if (rep.error == 0) break;
            }
            double ms = (System.nanoTime() - t1) / 1e6;
            sumMs += ms; maxMs = Math.max(maxMs, ms);
            sumAttempts += attempts; maxAttempts = Math.max(maxAttempts, attempts);
            if (rep.size == in.N && rep.invalidRows == 0) sizeOk++;
            if (rep.error == 0) { exact++; sumDist += rep.distinct; sumH += rep.entropyNorm; if (attempts == 1) firstTry++; }
            if (r == 0) {
                StringBuilder sb = new StringBuilder();
                for (int a : g.lastOrder) sb.append(in.attrNames[a]).append(' ');
                System.out.println("ordre (run 0) : " + sb);
            }
            if (!quiet)
                System.out.printf(Locale.ROOT, "run %2d  %5.0f ms  tentatives=%d  retours à 2+ colonnes=%d  taille=%d  erreur=%3d  distinctes=%3d  en doublon=%3d  copies en trop=%3d  H=%.4f (%.2f%%)%n",
                        r, ms, attempts, deep, rep.size, rep.error, rep.distinct, rep.duplicatedIndividuals, rep.extraCopies(), rep.entropy, 100 * rep.entropyNorm);
            if (bestRep == null || Main.better(rep, bestRep)) { bestRep = rep; best = pop; }
        }
        System.out.printf(Locale.ROOT, "runs exacts : %d/%d (%.1f%%), dont %d du premier coup ; tentatives moy. %.2f (max %d) ; temps moyen %.0f ms (max %.0f)%n",
                exact, runs, 100.0 * exact / runs, firstTry, (double) sumAttempts / runs, maxAttempts, sumMs / runs, maxMs);
        if (exact > 0) System.out.printf(Locale.ROOT, "sur les exacts : distinctes moy. %.1f, H moy. %.2f%%%n", sumDist / exact, 100 * sumH / exact);
        System.out.printf("taille = N = %d : %d/%d runs%n", in.N, sizeOk, runs);
        System.out.println("meilleure : " + bestRep);
        Main.write(in, best, Paths.get(out));
    }
}
