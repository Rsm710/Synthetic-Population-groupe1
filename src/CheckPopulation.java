import java.nio.file.*;
import java.util.*;

/**
 * Vérifie un fichier population écrit sur disque, indépendamment du générateur.
 * Usage : java -cp out CheckPopulation <contraintes.csv> <population.csv>
 */
public final class CheckPopulation {
    public static void main(String[] args) throws Exception {
        Instance in = Instance.load(Paths.get(args[0]));
        List<String> lines = Files.readAllLines(Paths.get(args[1]));
        String[] header = lines.get(0).split(",");
        int[] col = new int[in.K];                     // colonne du fichier pour chaque attribut
        Arrays.fill(col, -1);
        for (int j = 0; j < header.length; j++)
            for (int a = 0; a < in.K; a++) if (header[j].trim().equals(in.attrNames[a])) col[a] = j;
        for (int a = 0; a < in.K; a++)
            if (col[a] < 0) throw new IllegalStateException("Attribut absent du fichier : " + in.attrNames[a]);
        List<int[]> pop = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) {
            if (lines.get(i).isBlank()) continue;
            String[] f = lines.get(i).split(",");
            int[] x = new int[in.K];
            for (int a = 0; a < in.K; a++) x[a] = Arrays.asList(in.valueNames[a]).indexOf(f[col[a]].trim());
            pop.add(x);
        }
        Evaluation.Report rep = Evaluation.evaluate(in, pop.toArray(new int[0][]));
        System.out.println(args[1] + " : N attendu=" + in.N + " | " + rep);
    }
}
