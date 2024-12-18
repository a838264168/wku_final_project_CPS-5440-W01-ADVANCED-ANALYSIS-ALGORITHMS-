import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

import java.awt.Color;
import javax.swing.*;
import java.util.*;

public class DijkstraComparison extends ApplicationFrame {
    private static class Node implements Comparable<Node> {
        int vertex;
        int distance;

        Node(int vertex, int distance) {
            this.vertex = vertex;
            this.distance = distance;
        }

        @Override
        public int compareTo(Node other) {
            return Integer.compare(this.distance, other.distance);
        }
    }

    public static Result dijkstra(int[][] graph, int start, int end) {
        int n = graph.length;
        int[] distances = new int[n];
        boolean[] visited = new boolean[n];
        int nodesExtended = 0;

        Arrays.fill(distances, Integer.MAX_VALUE);
        distances[start] = 0;

        PriorityQueue<Node> pq = new PriorityQueue<>();
        pq.offer(new Node(start, 0));

        while (!pq.isEmpty()) {
            Node current = pq.poll();
            int currentVertex = current.vertex;

            if (currentVertex == end) {
                return new Result(distances[end], nodesExtended);
            }

            if (visited[currentVertex]) continue;
            visited[currentVertex] = true;
            nodesExtended++;

            for (int neighbor = 0; neighbor < n; neighbor++) {
                if (graph[currentVertex][neighbor] > 0) {
                    int newDistance = distances[currentVertex] + graph[currentVertex][neighbor];
                    if (newDistance < distances[neighbor]) {
                        distances[neighbor] = newDistance;
                        pq.offer(new Node(neighbor, newDistance));
                    }
                }
            }
        }

        return new Result(-1, nodesExtended); // Path not found
    }

    private static int[][] generateGraph(int size, boolean isSparse) {
        int[][] graph = new int[size][size];
        Random random = new Random();
        int edgesPerVertex = isSparse ? 2 : size / 2;

        for (int i = 0; i < size; i++) {
            Set<Integer> connectedVertices = new HashSet<>();
            while (connectedVertices.size() < edgesPerVertex) {
                int j = random.nextInt(size);
                if (i != j && !connectedVertices.contains(j)) {
                    int weight = random.nextInt(100) + 1;
                    graph[i][j] = weight;
                    connectedVertices.add(j);
                }
            }
        }

        return graph;
    }

    private static class Result {
        int pathLength;
        int nodesExtended;

        Result(int pathLength, int nodesExtended) {
            this.pathLength = pathLength;
            this.nodesExtended = nodesExtended;
        }
    }

    public DijkstraComparison(String title) {
        super(title);
        JPanel mainPanel = new JPanel();
        mainPanel.setLayout(new BoxLayout(mainPanel, BoxLayout.Y_AXIS));
        mainPanel.add(createChartPanel(true));
        mainPanel.add(createChartPanel(false));
        setContentPane(mainPanel);
    }

    private JPanel createChartPanel(boolean isSparse) {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        int[] sizes = {100, 1000, 10000};
        for (int size : sizes) {
            int[][] graph = generateGraph(size, isSparse);

            long startTime = System.nanoTime();
            Result result = dijkstra(graph, 0, size - 1);
            long time = (System.nanoTime() - startTime) / 1_000_000;

            dataset.addValue(time, "Time (ms)", String.valueOf(size));
            dataset.addValue(result.nodesExtended, "Nodes Extended", String.valueOf(size));
            dataset.addValue(result.pathLength, "Path Length", String.valueOf(size));

            System.out.println("Graph Size: " + size + (isSparse ? " (Sparse)" : " (Dense)"));
            System.out.println("Time: " + time + "ms, Nodes Extended: " + result.nodesExtended + ", Path Length: " + result.pathLength);
           
            System.out.println();
        }

        String chartTitle = isSparse ? "Dijkstra's Algorithm: Sparse Graph" : "Dijkstra's Algorithm: Dense Graph";
        JFreeChart chart = ChartFactory.createBarChart(
            chartTitle,
            "Graph Size",
            "Value",
            dataset,
            PlotOrientation.VERTICAL,
            true, true, false);

        CategoryPlot plot = (CategoryPlot) chart.getPlot();
        BarRenderer renderer = (BarRenderer) plot.getRenderer();

        // Set colors for different series
        renderer.setSeriesPaint(0, Color.RED);
        renderer.setSeriesPaint(1, Color.BLUE);
        renderer.setSeriesPaint(2, Color.GREEN);

        ChartPanel chartPanel = new ChartPanel(chart);
        chartPanel.setPreferredSize(new java.awt.Dimension(800, 400));
        return chartPanel;
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            DijkstraComparison example = new DijkstraComparison("Dijkstra's Algorithm Comparison");
            example.pack();
            RefineryUtilities.centerFrameOnScreen(example);
            example.setVisible(true);
        });
    }
}

