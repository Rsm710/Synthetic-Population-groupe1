import java.nio.file.*;
import java.util.*;

/**
 * Usage : java -cp out GreedyMain <contraintes.csv> [--runs R] [--seed s]
 *                                 [--fc on|off] [--retries r] [--backtrack b] [--depth d] [--restarts m] [--repair iters] [--perfect bonus] [--cross lambda] [--beta b] [--greedy 1|2] [--gamma g] [--entropy off|full|loo] [--lambda l] [--alpha a] [--declone ms] [--dentropy off|full|loo|near] [--dradius r] [--dmu mu] [--cycle p] [--deg4 p] [--dfmax f] [--hguard on|off] [--neutral p] [--temp T] [--concentrate on|off] [--small s] [--coldeclone ms] [--colfrom p] [--colneutral p] [--colspread on|off] [--colbound on|off] [--starts m] [--out population.csv] [--quiet]
 */
public final class GreedyMain {
    public static void main(String[] args) throws Exception {
        Instance in = Instance.load(Paths.get(args[0]));
        int runs = 10; long seed = 0; boolean fc = true, quiet = false;
        int greedy = 1; double gamma = 1.0; int entropy = 1; double lambda = 1.0, alpha = 0.3;
        double beta = 2.0; int retries = 3, backtrack = 40, repair = 500, depth = 2, restarts = 10; double cross = 0, perfect = 1000;
        long declone = 1500; int dentropy = 2, dradius = 1; double dmu = -1, cycle = 0, deg4 = 0.2; int dfmax = 1; boolean hguard = false; double neutral = 1.0, temp = 0, small = 0; int starts = 1; boolean concentrate = false;
        long colDeclone = 0; int colFrom = 1; double colNeutral = 1.0; boolean colSpread = false, colBound = false;
        String out = "population.csv";
        for (int k = 1; k < args.length; k++) {
            switch (args[k]) {
                case "--runs" -> runs = Integer.parseInt(args[++k]);
                case "--seed" -> seed = Long.parseLong(args[++k]);
                case "--fc" -> fc = args[++k].equals("on");
                case "--beta" -> beta = Double.parseDouble(args[++k]);
                case "--greedy" -> greedy = Integer.parseInt(args[++k]);
                case "--gamma" -> gamma = Double.parseDouble(args[++k]);
                case "--entropy" -> entropy = switch (args[++k]) { case "off" -> 0; case "full" -> 1; case "loo" -> 2;
                    default -> throw new IllegalArgumentException("--entropy off|full|loo"); };
                case "--lambda" -> lambda = Double.parseDouble(args[++k]);
                case "--alpha" -> alpha = Double.parseDouble(args[++k]);
                case "--out" -> out = args[++k];
                case "--quiet" -> quiet = true;
                case "--declone" -> declone = Long.parseLong(args[++k]);
                case "--dentropy" -> dentropy = switch (args[++k]) { case "off" -> 0; case "full" -> 1; case "loo" -> 2; case "near" -> 3;
                    default -> throw new IllegalArgumentException("--dentropy off|full|loo|near"); };
                case "--dradius" -> dradius = Integer.parseInt(args[++k]);
                case "--dmu" -> dmu = Double.parseDouble(args[++k]);
                case "--cycle" -> cycle = Double.parseDouble(args[++k]);
                case "--deg4" -> deg4 = Double.parseDouble(args[++k]);
                case "--dfmax" -> dfmax = Integer.parseInt(args[++k]);
                case "--hguard" -> hguard = args[++k].equals("on");
                case "--neutral" -> neutral = Double.parseDouble(args[++k]);
                case "--temp" -> temp = Double.parseDouble(args[++k]);
                case "--concentrate" -> concentrate = args[++k].equals("on");
                case "--small" -> small = Double.parseDouble(args[++k]);
                case "--starts" -> starts = Integer.parseInt(args[++k]);
                case "--coldeclone" -> colDeclone = Long.parseLong(args[++k]);
                case "--colfrom" -> colFrom = Integer.parseInt(args[++k]);
                case "--colneutral" -> colNeutral = Double.parseDouble(args[++k]);
                case "--colspread" -> colSpread = args[++k].equals("on");
                case "--colbound" -> colBound = args[++k].equals("on");
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
        System.out.printf("Instance : N=%d, K=%d, %d tables, greedy %d, entropie %s%n", in.N, in.K, in.T, greedy,
                entropy == 0 ? "off" : (entropy == 1 ? "full" : "loo (lambda=" + lambda + ")") + " alpha=" + alpha);
        // Un run = starts constructions ; chacune = au plus 1 + restarts tentatives complètes (même générateur aléatoire, tirages différents),
        // arrêtées dès qu'une population est exacte selon Evaluation.evaluate.
        int exact = 0, firstTry = 0, sizeOk = 0, sumAttempts = 0, maxAttempts = 0; int[][] best = null; Evaluation.Report bestRep = null;
        double sumMs = 0, maxMs = 0, sumDist = 0, sumH = 0, sumNear = 0, sumMI = 0, sumD1 = 0, sumD0 = 0, sumDef = 0, sumDef2 = 0, sumD2 = 0, sumH0 = 0, sumNear0 = 0, sumD10 = 0;
        for (int r = 0; r < runs; r++) {
            Random rnd = new Random(seed + r);
            long t1 = System.nanoTime();
            // starts constructions indépendantes (greedy puis dé-clonage) ; on garde la meilleure.
            int[][] pop = null; Evaluation.Report rep = null; ColumnGreedy g = null;
            int attempts = 0, deep = 0, before = 0; Evaluation.Report built = null;   // built : population avant dé-clonage
            for (int s = 0; s < starts; s++) {
                int[][] p = null; Evaluation.Report rp = null; int att = 0;
                while (att <= restarts) {
                    att++;
                    if (greedy == 2) { ColumnGreedy2 g2 = new ColumnGreedy2(in, rnd); g2.gamma = gamma; g = g2; }
                    else g = new ColumnGreedy(in, rnd);
                    g.entropy = entropy; g.lambda = lambda; g.alpha = alpha;
                    g.beta = beta; g.lookAhead = fc; g.crossPenalty = cross; g.perfectBonus = perfect; g.columnRetries = retries;
                    g.backtrackTries = backtrack; g.backtrackDepth = depth; g.repairIters = repair;
                    g.columnDeclone = colDeclone; g.columnDecloneFrom = colFrom; g.columnNeutral = colNeutral; g.columnSpread = colSpread; g.columnBound = colBound;
                    p = g.build();
                    rp = Evaluation.evaluate(in, p);
                    deep += g.deepBacktracks;
                    if (rp.error == 0) break;
                }
                attempts = Math.max(attempts, att);
                int b0 = rp.distinct; Evaluation.Report r0 = rp;
                if (rp.error == 0 && declone > 0) {          // post-optimisation : échanges qui gardent l'exactitude
                    Declone dc = new Declone(in, rnd);
                    dc.timeMs = declone; dc.neutral = neutral; dc.temp = temp; dc.concentrate = concentrate; dc.small = small; dc.entropy = dentropy; dc.radius = dradius; dc.mu = dmu; dc.cycle = cycle; dc.deg4 = deg4; dc.fmax = dfmax; dc.hguard = hguard;
                    dc.run(p);
                    rp = Evaluation.evaluate(in, p);
                }
                if (rep == null || Main.better(rp, rep)) { rep = rp; pop = p; before = b0; built = r0; }
            }
            double ms = (System.nanoTime() - t1) / 1e6;
            sumMs += ms; maxMs = Math.max(maxMs, ms);
            sumAttempts += attempts; maxAttempts = Math.max(maxAttempts, attempts);
            if (rep.size == in.N && rep.invalidRows == 0) sizeOk++;
            if (rep.error == 0) {
                exact++; sumDist += rep.distinct; sumH += rep.entropy; sumNear += rep.nearDup; sumMI += rep.freeMI; sumD1 += rep.pairsD1; sumD0 += rep.pairsD0; sumDef += rep.looDeficit; sumDef2 += rep.looDeficit2; sumD2 += rep.pairsD2;
                sumH0 += built.entropy; sumNear0 += built.nearDup; sumD10 += built.pairsD1;
                if (attempts == 1) firstTry++;
            }
            if (r == 0) {
                StringBuilder sb = new StringBuilder();
                for (int a : g.lastOrder) sb.append(in.attrNames[a]).append(' ');
                System.out.println("ordre (run 0) : " + sb);
            }
            if (!quiet)
                System.out.printf(Locale.ROOT, "run %2d  %5.0f ms  tentatives=%d  retours à 2+ colonnes=%d  taille=%d  erreur=%3d  distinctes=%3d -> %3d  en doublon=%3d  copies en trop=%3d  H=%.4f -> %.4f  ΣH(b|-b)=%.3f -> %.3f  d1=%d -> %d  I libres=%.4f%n",
                        r, ms, attempts, deep, rep.size, rep.error, before, rep.distinct, rep.duplicatedIndividuals, rep.extraCopies(),
                        built.entropy, rep.entropy, built.nearDup, rep.nearDup, built.pairsD1, rep.pairsD1, rep.freeMI);
            if (bestRep == null || Main.better(rep, bestRep)) { bestRep = rep; best = pop; }
        }
        System.out.printf(Locale.ROOT, "runs exacts : %d/%d (%.1f%%), dont %d du premier coup ; tentatives moy. %.2f (max %d) ; temps moyen %.0f ms (max %.0f)%n",
                exact, runs, 100.0 * exact / runs, firstTry, (double) sumAttempts / runs, maxAttempts, sumMs / runs, maxMs);
        if (exact > 0) {
            System.out.printf(Locale.ROOT, "sur les exacts : distinctes moy. %.1f, H moy. %.4f nats (%.2f%% de ln N), déficit K ln N - Σ H(X_-b) moy. %.3f, Σ H(X_b|X_-b) moy. %.3f, paires à distance 0 / 1 / 2 moy. %.1f / %.1f / %.1f, déficit rayon 2 moy. %.3f, Σ I paires libres moy. %.4f%n",
                    sumDist / exact, sumH / exact, 100 * sumH / exact / Math.log(in.N), sumDef / exact, sumNear / exact, sumD0 / exact, sumD1 / exact, sumD2 / exact, sumDef2 / exact, sumMI / exact);
            System.out.printf(Locale.ROOT, "  construction seule (avant dé-clonage) : H moy. %.4f, Σ H(X_b|X_-b) moy. %.3f, paires à distance 1 moy. %.1f%n",
                    sumH0 / exact, sumNear0 / exact, sumD10 / exact);
        }
        System.out.printf("taille = N = %d : %d/%d runs%n", in.N, sizeOk, runs);
        System.out.println("meilleure : " + bestRep);
        Main.write(in, best, Paths.get(out));
    }
}
