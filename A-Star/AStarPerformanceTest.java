import java.util.*;
import java.util.List;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.chart.axis.LogarithmicAxis;
import org.jfree.chart.plot.CategoryPlot;
import javax.swing.*;
import java.awt.*;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.plot.XYPlot;


public class AStarPerformanceTest {
    static class Vertex {
        int id;
        int x, y;
        long distance;
        long potential;
        long distwithPotential;
        boolean processed;
        List<Edge> edges;

        public Vertex(int id, int x, int y) {
            this.id = id;
            this.x = x;
            this.y = y;
            this.distance = Long.MAX_VALUE;
            this.potential = 0;
            this.distwithPotential = Long.MAX_VALUE;
            this.processed = false;
            this.edges = new ArrayList<>();
        }
    }

    static class Edge {
        Vertex to;
        int weight;

        public Edge(Vertex to, int weight) {
            this.to = to;
            this.weight = weight;
        }
    }

    static class Result {
        long distance;
        int expandedNodes;
        long time;
        double theoreticalComplexity;

        public Result(long distance, int expandedNodes, long time, double theoreticalComplexity) {
            this.distance = distance;
            this.expandedNodes = expandedNodes;
            this.time = time;
            this.theoreticalComplexity = theoreticalComplexity;
        }
    }

    public static Vertex[] generateSparseGraph(int n, int avgDegree) {
        Vertex[] graph = new Vertex[n];
        Random rand = new Random();

        for (int i = 0; i < n; i++) {
            graph[i] = new Vertex(i, rand.nextInt(1000), rand.nextInt(1000));
        }

        for (int i = 0; i < n; i++) {
            int degree = rand.nextInt(avgDegree * 2) + 1;
            for (int j = 0; j < degree; j++) {
                int to = rand.nextInt(n);
                if (to != i && !hasEdge(graph[i], graph[to])) {
                    int weight = rand.nextInt(100) + 1;
                    graph[i].edges.add(new Edge(graph[to], weight));
                }
            }
        }

        return graph;
    }

    public static Vertex[] generateDenseGraph(int n) {
        Vertex[] graph = new Vertex[n];
        Random rand = new Random();

        for (int i = 0; i < n; i++) {
            graph[i] = new Vertex(i, rand.nextInt(1000), rand.nextInt(1000));
        }

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    int weight = rand.nextInt(100) + 1;
                    graph[i].edges.add(new Edge(graph[j], weight));
                }
            }
        }

        return graph;
    }

    private static boolean hasEdge(Vertex from, Vertex to) {
        for (Edge edge : from.edges) {
            if (edge.to == to) return true;
        }
        return false;
    }

    public static Result aStarSearch(Vertex[] graph, int start, int goal, boolean isSparse) {
        long startTime = System.nanoTime();
        int expandedNodes = 0;

        PriorityQueue<Vertex> openSet = new PriorityQueue<>(
            Comparator.comparingLong(v -> v.distwithPotential)
        );

        graph[start].distance = 0;
        graph[start].potential = heuristic(graph[start], graph[goal]);
        graph[start].distwithPotential = graph[start].distance + graph[start].potential;
        openSet.add(graph[start]);

        while (!openSet.isEmpty()) {
            Vertex current = openSet.poll();
            expandedNodes++;

            if (current.id == goal) {
                long endTime = System.nanoTime();
                long duration = (endTime - startTime) / 1_000_000;
                double theoreticalComplexity = isSparse ? 
                    calculateSparseComplexity(graph.length) : 
                    calculateDenseComplexity(graph.length);
                return new Result(current.distance, expandedNodes, duration, theoreticalComplexity);
            }

            current.processed = true;

            for (Edge edge : current.edges) {
                Vertex neighbor = edge.to;
                if (neighbor.processed) continue;

                long tentativeDistance = current.distance + edge.weight;

                if (tentativeDistance < neighbor.distance) {
                    neighbor.distance = tentativeDistance;
                    neighbor.potential = heuristic(neighbor, graph[goal]);
                    neighbor.distwithPotential = neighbor.distance + neighbor.potential;

                    if (!openSet.contains(neighbor)) {
                        openSet.add(neighbor);
                    } else {
                        openSet.remove(neighbor);
                        openSet.add(neighbor);
                    }
                }
            }
        }

        long endTime = System.nanoTime();
        long duration = (endTime - startTime) / 1_000_000;
        double theoreticalComplexity = isSparse ? 
            calculateSparseComplexity(graph.length) : 
            calculateDenseComplexity(graph.length);
        return new Result(-1, expandedNodes, duration, theoreticalComplexity);
    }

    private static long heuristic(Vertex a, Vertex b) {
        return (long) Math.sqrt(Math.pow(a.x - b.x, 2) + Math.pow(a.y - b.y, 2));
    }

    private static double calculateSparseComplexity(int v) {
        return v * Math.log(v);
    }

    private static double calculateDenseComplexity(int v) {
        return Math.pow(v, 2)* Math.log(v);
    }

    private static void createPerformanceChart(String title, List<Result> results, int[] sizes) {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        double smallPositiveValue = 0.1;

        for (int i = 0; i < sizes.length; i++) {
            Result result = results.get(i);
            int size = sizes[i];

            dataset.addValue(Math.max(result.time, smallPositiveValue), "Actual Time (ms)", String.valueOf(size));
            dataset.addValue(Math.max(result.expandedNodes, smallPositiveValue), "Expanded Nodes", String.valueOf(size));
            dataset.addValue(Math.max(result.distance, smallPositiveValue), "Distance", String.valueOf(size));
            dataset.addValue(Math.max(result.theoreticalComplexity, smallPositiveValue), "Theoretical Complexity", String.valueOf(size));
        }

        JFreeChart chart = ChartFactory.createBarChart(
                title,
                "Graph Size",
                "Value (log scale)",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false);

        CategoryPlot plot = (CategoryPlot) chart.getPlot();
        LogarithmicAxis logAxis = new LogarithmicAxis("Value (log scale)");
        logAxis.setAllowNegativesFlag(true);
        plot.setRangeAxis(logAxis);

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    private static void createComparisonChart(String title, List<Result> results, int[] sizes) {
        XYSeriesCollection dataset = new XYSeriesCollection();
        XYSeries actualTimeSeries = new XYSeries("Actual Time");
        XYSeries theoreticalComplexitySeries = new XYSeries("Theoretical Complexity");

        double smallPositiveValue = 0.1;

        for (int i = 0; i < sizes.length; i++) {
            Result result = results.get(i);
            int size = sizes[i];
            actualTimeSeries.add(size, Math.max(result.time, smallPositiveValue));
            theoreticalComplexitySeries.add(size, Math.max(result.theoreticalComplexity, smallPositiveValue));
        }

        dataset.addSeries(actualTimeSeries);
        dataset.addSeries(theoreticalComplexitySeries);

        JFreeChart chart = ChartFactory.createXYLineChart(
                title,
                "Graph Size",
                "Time / Complexity",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false);

        XYPlot plot = (XYPlot) chart.getPlot();
        LogarithmicAxis xAxis = new LogarithmicAxis("Graph Size (log scale)");
        LogarithmicAxis yAxis = new LogarithmicAxis("Time / Complexity (log scale)");
        plot.setDomainAxis(xAxis);
        plot.setRangeAxis(yAxis);

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        JFrame frame = new JFrame(title);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(chartPanel);
        frame.pack();
        frame.setVisible(true);
    }

    public static void main(String[] args) {
        int[] sizes = {100, 1000, 10000};
        List<Result> sparseResults = new ArrayList<>();
        List<Result> denseResults = new ArrayList<>();

        System.out.println("Testing Sparse Graphs");
        for (int size : sizes) {
            Vertex[] sparseGraph = generateSparseGraph(size, 5);
            Result result = aStarSearch(sparseGraph, 0, size - 1, true);
            sparseResults.add(result);
            System.out.printf("Size: %d, Time: %d ms, Expanded Nodes: %d, Distance: %d, Theoretical Complexity: %.2f%n", 
                              size, result.time, result.expandedNodes, result.distance, result.theoreticalComplexity);
        }

        System.out.println("\nTesting Dense Graphs");
        for (int size : sizes) {
            Vertex[] denseGraph = generateDenseGraph(size);
            Result result = aStarSearch(denseGraph, 0, size - 1, false);
            denseResults.add(result);
            System.out.printf("Size: %d, Time: %d ms, Expanded Nodes: %d, Distance: %d, Theoretical Complexity: %.2f%n", 
                              size, result.time, result.expandedNodes, result.distance, result.theoreticalComplexity);
        }

        createPerformanceChart("A* Performance on Sparse Graphs", sparseResults, sizes);
        createPerformanceChart("A* Performance on Dense Graphs", denseResults, sizes);
        createComparisonChart("A* Time vs Theoretical Complexity (Sparse)", sparseResults, sizes);
        createComparisonChart("A* Time vs Theoretical Complexity (Dense)", denseResults, sizes);
    }
}

