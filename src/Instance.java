import java.io.*;
import java.nio.file.*;
import java.util.*;

/**
 * Instance générique : K attributs catégoriels, domaines quelconques,
 * tables de marginales d'arité quelconque. Tout est déduit du CSV
 * (colonnes : cell_id,table_id,arity,scope,cell,target).
 */
public final class Instance {
    public final int K;                 // nombre d'attributs
    public final int N;                 // taille de la population
    public final String[] attrNames;    // nom de chaque attribut
    public final String[][] valueNames; // valueNames[a][v]
    public final int[] dom;             // taille du domaine de chaque attribut

    public final int T;                 // nombre de tables
    public final int[] tableIds;        // id d'origine des tables
    public final int[][] scope;         // scope[t] = attributs de la table t
    public final int[][] stride;        // index de cellule = somme v[scope[j]] * stride[t][j]
    public final int[][] target;        // target[t][cellule]
    public final int[][] attrTables;    // attrTables[a] = tables contenant a

    private Instance(int K, int N, String[] attrNames, String[][] valueNames, int[] dom,
                     int[] tableIds, int[][] scope, int[][] stride, int[][] target) {
        this.K = K; this.N = N; this.attrNames = attrNames; this.valueNames = valueNames; this.dom = dom;
        this.T = scope.length; this.tableIds = tableIds; this.scope = scope; this.stride = stride; this.target = target;
        List<List<Integer>> at = new ArrayList<>();
        for (int a = 0; a < K; a++) at.add(new ArrayList<>());
        for (int t = 0; t < T; t++) for (int a : scope[t]) at.get(a).add(t);
        attrTables = new int[K][];
        for (int a = 0; a < K; a++) attrTables[a] = at.get(a).stream().mapToInt(Integer::intValue).toArray();
    }

    /** Index de la cellule de la table t pour un individu (valeurs complètes). */
    public int cellOf(int t, int[] values) {
        int idx = 0;
        for (int j = 0; j < scope[t].length; j++) idx += values[scope[t][j]] * stride[t][j];
        return idx;
    }

    public static Instance load(Path csv) throws IOException {
        List<String> lines = Files.readAllLines(csv);
        // 1) collecte des noms d'attributs et de valeurs
        Map<String, TreeSet<String>> values = new LinkedHashMap<>();
        Map<Integer, String> scopeOf = new LinkedHashMap<>();
        List<String[]> rows = new ArrayList<>();
        for (int i = 1; i < lines.size(); i++) {
            String line = lines.get(i).trim();
            if (line.isEmpty()) continue;
            String[] f = line.split(",");
            rows.add(f);
            int tid = Integer.parseInt(f[1].trim());
            scopeOf.putIfAbsent(tid, f[3].trim());
            for (String kv : f[4].trim().split(";")) {
                String[] p = kv.split("=");
                values.computeIfAbsent(p[0].trim(), k -> new TreeSet<>(NATURAL)).add(p[1].trim());
            }
        }
        List<String> attrs = new ArrayList<>(values.keySet());
        attrs.sort(NATURAL);
        int K = attrs.size();
        Map<String, Integer> attrIdx = new HashMap<>();
        String[] attrNames = new String[K];
        String[][] valueNames = new String[K][];
        int[] dom = new int[K];
        List<Map<String, Integer>> valIdx = new ArrayList<>();
        for (int a = 0; a < K; a++) {
            attrNames[a] = attrs.get(a);
            attrIdx.put(attrs.get(a), a);
            valueNames[a] = values.get(attrs.get(a)).toArray(new String[0]);
            dom[a] = valueNames[a].length;
            Map<String, Integer> m = new HashMap<>();
            for (int v = 0; v < dom[a]; v++) m.put(valueNames[a][v], v);
            valIdx.add(m);
        }
        // 2) tables
        List<Integer> tids = new ArrayList<>(scopeOf.keySet());
        Collections.sort(tids);
        int T = tids.size();
        Map<Integer, Integer> tIdx = new HashMap<>();
        int[] tableIds = new int[T];
        int[][] scope = new int[T][], stride = new int[T][], target = new int[T][];
        for (int t = 0; t < T; t++) {
            int tid = tids.get(t);
            tableIds[t] = tid;
            tIdx.put(tid, t);
            String[] sc = scopeOf.get(tid).split("\\|");
            scope[t] = new int[sc.length];
            stride[t] = new int[sc.length];
            int size = 1;
            for (int j = sc.length - 1; j >= 0; j--) {
                scope[t][j] = attrIdx.get(sc[j].trim());
                stride[t][j] = size;
                size *= dom[scope[t][j]];
            }
            target[t] = new int[size];
            Arrays.fill(target[t], -1);
        }
        for (String[] f : rows) {
            int t = tIdx.get(Integer.parseInt(f[1].trim()));
            int[] vals = new int[K];
            for (String kv : f[4].trim().split(";")) {
                String[] p = kv.split("=");
                int a = attrIdx.get(p[0].trim());
                vals[a] = valIdx.get(a).get(p[1].trim());
            }
            int idx = 0;
            for (int j = 0; j < scope[t].length; j++) idx += vals[scope[t][j]] * stride[t][j];
            target[t][idx] = Integer.parseInt(f[5].trim());
        }
        // 3) cohérence : cellules manquantes = 0, toutes les tables somment à N
        int N = -1;
        for (int t = 0; t < T; t++) {
            int s = 0;
            for (int c = 0; c < target[t].length; c++) {
                if (target[t][c] < 0) target[t][c] = 0;
                s += target[t][c];
            }
            if (N < 0) N = s;
            else if (s != N)
                System.err.println("Attention : la table " + tableIds[t] + " somme à " + s + " au lieu de " + N);
        }
        return new Instance(K, N, attrNames, valueNames, dom, tableIds, scope, stride, target);
    }

    /** Tri naturel : A2 < A10, v2 < v10. */
    static final Comparator<String> NATURAL = (x, y) -> {
        String px = x.replaceAll("\\d+$", ""), py = y.replaceAll("\\d+$", "");
        int c = px.compareTo(py);
        if (c != 0) return x.compareTo(y);
        String nx = x.substring(px.length()), ny = y.substring(py.length());
        if (nx.isEmpty() || ny.isEmpty()) return x.compareTo(y);
        return Long.compare(Long.parseLong(nx), Long.parseLong(ny));
    };
}
