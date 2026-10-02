import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Usage : java -cp out Main <contraintes.csv> [--runs R] [--time ms_par_run] [--iters n]
 *                                            [--seed s] [--out population.csv]
 */
public final class Main {
    public static void main(String[] args) throws IOException {
        Path csv = Paths.get(args[0]);
        int runs = 5; long timeMs = 5_000, iters = 3_000_000, seed = 42; double W = 3.0;
        String out = "population.csv";
        for (int k = 1; k < args.length; k += 2) {
            switch (args[k]) {
                case "--runs" -> runs = Integer.parseInt(args[k + 1]);
                case "--time" -> timeMs = Long.parseLong(args[k + 1]);
                case "--iters" -> iters = Long.parseLong(args[k + 1]);
                case "--W" -> W = Double.parseDouble(args[k + 1]);
                case "--seed" -> seed = Long.parseLong(args[k + 1]);
                case "--out" -> out = args[k + 1];
                default -> throw new IllegalArgumentException("Option inconnue : " + args[k]);
            }
        }
        Instance in = Instance.load(csv);
        System.out.printf("Instance : N=%d, K=%d, domaines=%s, %d tables%n",
                in.N, in.K, Arrays.toString(in.dom), in.T);
        Evaluation.printBounds(in, 50_000_000L);

        int[][] bestPop = null; Evaluation.Report bestRep = null;
        for (int r = 0; r < runs; r++) {
            Random rnd = new Random(seed + r);
            long t0 = System.nanoTime();
            GreedyBuilder gb = new GreedyBuilder(in, rnd);
            Population p = gb.build();
            long t1 = System.nanoTime();
            Evaluation.Report g = Evaluation.evaluate(in, p.snapshot());
            LocalSearch ls = new LocalSearch(in, rnd);
            ls.maxIters = iters; ls.timeLimitMs = timeMs; ls.W = W;
            int[][] pop = ls.run(p);
            long t2 = System.nanoTime();
            Evaluation.Report rep = Evaluation.evaluate(in, pop);
            System.out.printf(Locale.ROOT, "Run %d%n  greedy (%.0f ms, %d impasses) : %s%n  + LS   (%.0f ms, %d it.) : %s%n",
                    r, (t1 - t0) / 1e6, gb.deadEnds, g, (t2 - t1) / 1e6, ls.iters, rep);
            if (bestRep == null || better(rep, bestRep)) { bestRep = rep; bestPop = pop; }
        }
        System.out.println("Meilleure solution : " + bestRep);
        write(in, bestPop, Paths.get(out));
        System.out.println("Population écrite dans " + out);
    }

    static boolean better(Evaluation.Report a, Evaluation.Report b) {
        if (a.error != b.error) return a.error < b.error;
        if (a.distinct != b.distinct) return a.distinct > b.distinct;
        return a.entropy > b.entropy;
    }

    static void write(Instance in, int[][] pop, Path path) throws IOException {
        try (PrintWriter w = new PrintWriter(Files.newBufferedWriter(path))) {
            StringBuilder h = new StringBuilder("id");
            for (String a : in.attrNames) h.append(',').append(a);
            w.println(h);
            for (int i = 0; i < pop.length; i++) {
                StringBuilder sb = new StringBuilder().append(i);
                for (int a = 0; a < in.K; a++) sb.append(',').append(in.valueNames[a][pop[i][a]]);
                w.println(sb);
            }
        }
    }
}
