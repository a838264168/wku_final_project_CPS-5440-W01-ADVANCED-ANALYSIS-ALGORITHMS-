//Program to implement the Bi-Directional Dijkstra Algorithm.

import java.util.Scanner;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Comparator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import javafx.scene.chart.NumberAxis;

import javax.swing.JFrame;
public class FriendSuggestion {

    //class Vertex of Graph.
    static class Vertex{
	int vertexNum;                  //int vertexNum
	ArrayList<Integer> adjList;     //list of adjacent vertices.
	ArrayList<Integer> costList;    //list of cost or distance of adjacent vertices.

	int queuePos;                   //pos of vertex in the priorityqueue.
	long dist;                      //distance from start vetex.    
	boolean processed;              //is processed while traversing the graph.

	public Vertex(){
	}


	//Vertex Constructor.
	public Vertex(int vertexNum){
		this.vertexNum=vertexNum;
		this.adjList = new ArrayList<Integer>();
		this.costList = new ArrayList<Integer>();
	}


	//function to create the graph and reverse graph. forwPriorityQ for graph and revPriorityQ for reverse graph.
	public void createGraph(Vertex [] graph, Vertex [] reverseGraph, int [] forwPriorityQ, int [] revPriorityQ){
		for(int i=0;i<graph.length;i++){
			graph[i].queuePos = i;
			graph[i].processed = false;
			graph[i].dist = Integer.MAX_VALUE;
		
			reverseGraph[i].queuePos = i;
			reverseGraph[i].processed = false;
			reverseGraph[i].dist = Integer.MAX_VALUE;

			forwPriorityQ[i]=i;
			revPriorityQ[i]=i;
		}
	}
    }


    //Implementing PrioirtyQueue data structure by myself using min_heap property. 
    static class PriorityQueue{

	//function to swap elements in the PriorityQueue
	public void swap(Vertex [] graph, int [] priorityQ, int index1, int index2){
		int temp = priorityQ[index1];

		priorityQ[index1]=priorityQ[index2];
		graph[priorityQ[index2]].queuePos=index1;

		priorityQ[index2]=temp;
		graph[temp].queuePos=index2;
	}

	//function to swap start vertex and first vertex in the priorityQ.		
	public void makeQueue(Vertex [] graph,int [] forwpriorityQ, int source, int target){
		swap(graph, forwpriorityQ,0,source);
	}


	//function to extract the min value from the PriorityQueue	
	public int extractMin(Vertex [] graph, int [] priorityQ, int extractNum){
		int vertex = priorityQ[0];
		int size = priorityQ.length-1-extractNum;
		swap(graph,priorityQ,0,size);
		siftDown(0,graph,priorityQ,size);
		return vertex;
	}

	//function to siftdown the element at the given index in the PriorityQueue.
	public void siftDown(int index, Vertex [] graph, int [] priorityQ, int size){
		int min = index;
		if(2*index+1<size && graph[priorityQ[index]].dist > graph[priorityQ[2*index+1]].dist){
			min = 2*index+1;
		}
		if(2*index+2<size && graph[priorityQ[min]].dist > graph[priorityQ[2*index+2]].dist){
			min = 2*index+2;
		}
		if(min!=index){
			swap(graph,priorityQ,min,index);
			siftDown(min,graph,priorityQ,size);
		}
	}
	
	//function to change priority of an element.(priority can only be decreased.)
	public void changePriority(Vertex [] graph, int [] priorityQ, int index){
		if((index-1)/2 > -1 && graph[priorityQ[index]].dist < graph[priorityQ[(index-1)/2]].dist){
			swap(graph,priorityQ,index,(index-1)/2);
			changePriority(graph,priorityQ,(index-1)/2);
		}
	}
    }
   
    //function to relax edges i.e traverse only the adjacent vertices of the given vertex.
    private static void relaxEdges(Vertex [] graph, int vertex, int [] priorityQ, PriorityQueue queue,int queryId){
	ArrayList<Integer> vertexList = graph[vertex].adjList;   //get the adjacent vertices list.
	ArrayList<Integer> costList = graph[vertex].costList;    //get the cost list of adjacent vertices.
	graph[vertex].processed = true;  			 //mark processed true.
	
	for(int i=0;i<vertexList.size();i++){
		int temp = vertexList.get(i);
		int cost = costList.get(i);
	
		if(graph[temp].dist>graph[vertex].dist + cost){
			graph[temp].dist = graph[vertex].dist + cost;	
			queue.changePriority(graph,priorityQ,graph[temp].queuePos);
		}
	}
    }

    
    //function to compute distance between start vertex s and target vertex t.   
    public static long computeDistance(Vertex [] graph, Vertex [] reverseGraph, int s, int t,int queryId){
	
	//create two PriorityQueues forwQ for forward graph and revQ for reverse graph. 
	PriorityQueue queue = new PriorityQueue();
	int [] forwPriorityQ = new int[graph.length];  //for forward propagation.
	int [] revPriorityQ = new int[graph.length];   //for reverse propagation.

	//create graph.
	Vertex vertex = new Vertex();
	vertex.createGraph(graph,reverseGraph,forwPriorityQ,revPriorityQ);

	//dist of s from s is 0.
	//in rev graph dist of t from t is 0.
	graph[s].dist=0;
	reverseGraph[t].dist=0;
	queue.makeQueue(graph,forwPriorityQ,s,t);
	queue.makeQueue(reverseGraph,revPriorityQ,t,s);

	//store the processed vertices while traversing.
	ArrayList<Integer> forgraphprocessedVertices = new ArrayList<Integer>();  //for forward propagation.
	ArrayList<Integer> revgraphprocessedVertices = new ArrayList<Integer>();  //for reverse propagation.
	
	
	for(int i=0;i<graph.length;i++){

		//extract the vertex with min dist from forwQ.
		int vertex1 = queue.extractMin(graph,forwPriorityQ,i);
		if(graph[vertex1].dist==Integer.MAX_VALUE){
			continue;
		}

		//relax the edges of the extracted vertex.
		relaxEdges(graph,vertex1,forwPriorityQ,queue,queryId);

		//store into the processed vertices list.
		forgraphprocessedVertices.add(vertex1);

		//check if extratced vertex also processed in the reverse graph. If yes find the shortest distance.
		if(reverseGraph[vertex1].processed){
			return shortestPath(graph,reverseGraph,forgraphprocessedVertices,revgraphprocessedVertices,queryId);
		}
		

		//extract the vertex with min dist from revQ.
		int revVertex = queue.extractMin(reverseGraph,revPriorityQ,i);
		if(reverseGraph[revVertex].dist==Integer.MAX_VALUE){
			continue;
		}

		//relax the edges of the extracted vertex.
		relaxEdges(reverseGraph,revVertex,revPriorityQ,queue,queryId);

		//store in the processed vertices list of reverse graph.
		revgraphprocessedVertices.add(revVertex);
		
		//check if extracted vertex is also processed in the forward graph. If yes find the shortest distance.
		if(graph[revVertex].processed){
			return shortestPath(graph,reverseGraph,forgraphprocessedVertices,revgraphprocessedVertices,queryId);
		}
		
	}
	
	//if no path between s and t.
	return -1;
    }

	
    //function to find the shortest distance from the stored processed vertices of both forward and reverse graph.
    private static long shortestPath(Vertex [] graph, Vertex [] reverseGraph, ArrayList<Integer> forgraphprocessedVertices, ArrayList<Integer> revgraphprocessedVertices,int queryId){
	long distance = Integer.MAX_VALUE;
	
	//process the forward list.
	for(int i=0;i<forgraphprocessedVertices.size();i++){
		int vertex = forgraphprocessedVertices.get(i);
		if(reverseGraph[vertex].dist + graph[vertex].dist>=Integer.MAX_VALUE){
			continue;
		}
		long tempdist = graph[vertex].dist + reverseGraph[vertex].dist;	
		if(distance>tempdist){
			distance=tempdist;
		}
	}
	
	//process the reverse list.
	for(int i=0;i<revgraphprocessedVertices.size();i++){
		int vertex = revgraphprocessedVertices.get(i);
		if(reverseGraph[vertex].dist + graph[vertex].dist>=Integer.MAX_VALUE){
			continue;
		}
		long tempdist = reverseGraph[vertex].dist + graph[vertex].dist;
		if(distance>tempdist){
			distance=tempdist;
		}
	
	}
	return distance;
    }
    ////
    ////
    private static void printGraph(Vertex[] graph) {
        for (int i = 0; i < graph.length; i++) {
            System.out.println("Vertex " + i + ":");
            for (int j = 0; j < graph[i].adjList.size(); j++) {
                System.out.println("  -> " + graph[i].adjList.get(j) + " (Cost: " + graph[i].costList.get(j) + ")");
            }
        }
    }

    public static void main(String[] args) {
        int[] numVerticesArray = {100, 1000, 10000};
        String[] graphTypes = {"Sparse", "Dense"};
        long seed = System.nanoTime(); // 使用当前时间纳秒值作为种子

        for (String graphType : graphTypes) {
            double[] averageTimes = new double[numVerticesArray.length];
            double[] expandedNodes = new double[numVerticesArray.length];
            double[] distances = new double[numVerticesArray.length];
            double[] theoreticalComplexities = new double[numVerticesArray.length];

            System.out.println("Testing " + graphType.substring(0, 1).toUpperCase() + graphType.substring(1) + " Graphs");
            for (int i = 0; i < numVerticesArray.length; i++) {
                int numVertices = numVerticesArray[i];
                int numEdges = graphType.equals("Sparse") ? (numVertices / 10) : (numVertices * (numVertices - 1) / 2);
                double totalAverageTime = 0;
                double totalExpandedNodes = 0;
                double totalDistance = 0;

                for (int j = 0; j < 1; j++) { // 对每种大小的图测试1次
                    long localSeed = seed + i + j;
                    double[] results = simulateGraphType(graphType, numVertices, numEdges, true, localSeed);
                    totalAverageTime += results[0];
                    totalExpandedNodes += results[1];
                    totalDistance += results[2];
                }

                averageTimes[i] = totalAverageTime ; // 计算平均时间
                expandedNodes[i] = totalExpandedNodes ; // 计算平均扩展节点数
                distances[i] = totalDistance ; // 计算平均距离
                theoreticalComplexities[i] = calculateTheoreticalComplexity(numVertices, numEdges);

                // 打印结构化输出
                System.out.printf("Size: %d, Time: %.2f ms, Expanded Nodes: %.0f, Distance: %.0f, Theoretical Complexity: %.2f%n",
                                  numVertices, averageTimes[i], expandedNodes[i], distances[i], theoreticalComplexities[i]);
            }

            // 为当前图类型创建并显示图表
            createAndShowLineChart(graphType, numVerticesArray, averageTimes, theoreticalComplexities);
            createAndShowBarChart(graphType, numVerticesArray, averageTimes, expandedNodes, distances, theoreticalComplexities);
        }
    }

    private static double[] simulateGraphType(String graphType, int numVertices, int numEdges, boolean isWeighted, long seed) {
        // 生成随机图
        Vertex[] graph = new Vertex[numVertices];
        Vertex[] reverseGraph = new Vertex[numVertices];
        for (int i = 0; i < numVertices; i++) {
            graph[i] = new Vertex(i);
            reverseGraph[i] = new Vertex(i);
        }
        generateRandomGraph(graph, reverseGraph, numVertices, numEdges, isWeighted, seed);

        // 选择随机的起始顶点和目标顶点
        Random rand = new Random(seed);
        int source = rand.nextInt(numVertices);
        int target = rand.nextInt(numVertices);
        while (source == target) {
            target = rand.nextInt(numVertices);
        }

        // 执行双向 Dijkstra 算法
        long startTime = System.nanoTime();
        long distance = computeDistance(graph, reverseGraph, source, target, 0);
//        System.out.println(distance);
        long endTime = System.nanoTime();
        long duration = endTime - startTime;

        // 收集结果
        double[] results = new double[3];
        results[0] = duration / 1_000_000.0; // 转换为毫秒
        results[1] = countExpandedNodes(graph); // 扩展节点数
        results[2] = distance; // 最短路径距离

        return results;
    }

    private static int countExpandedNodes(Vertex[] graph) {
        int count = 0;
        for (Vertex vertex : graph) {
            if (vertex.processed) {
                count++;
            }
        }
        return count;
    }

    private static void generateRandomGraph(Vertex[] graph, Vertex[] reverseGraph, int numVertices, int numEdges, boolean isWeighted, long seed) {
        Random rand = new Random(seed);
        for (int i = 0; i < numEdges; i++) {
            int u = rand.nextInt(numVertices);
            int v = rand.nextInt(numVertices);
            // 避免自环
            while (u == v) {
                v = rand.nextInt(numVertices);
            }
            // 生成非负边权重
            int w = isWeighted ? rand.nextInt(10) + 1 : 1; // 确保权重至少为1
            graph[u].adjList.add(v);
            graph[u].costList.add(w);
            // 由于是无向图，需要添加反向边
            reverseGraph[v].adjList.add(u);
            reverseGraph[v].costList.add(w);
        	//get the edges.
//            printGraph(graph);

        }
    }
    private static double calculateTheoreticalComplexity(int numVertices, int numEdges) {
        // 双向Dijkstra算法的时间复杂度为O(ElogV  )
        return numEdges * Math.log(numVertices)  ;
    }

//    private static void createAndShowLineChart(String graphType, int[] numVerticesArray, double[] averageTimes, double[] theoreticalComplexities) {
//        XYSeries seriesTime = new XYSeries("Actual Time (ms)");
//        XYSeries seriesComplexity = new XYSeries("Theoretical Complexity");
//
//        for (int i = 0; i < numVerticesArray.length; i++) {
//            seriesTime.add(numVerticesArray[i], averageTimes[i]);
//            seriesComplexity.add(numVerticesArray[i], theoreticalComplexities[i]);
//        }
//
//        XYSeriesCollection dataset = new XYSeriesCollection();
//        dataset.addSeries(seriesTime);
//        dataset.addSeries(seriesComplexity);
//
//        JFreeChart chart = ChartFactory.createXYLineChart(
//                graphType.substring(0, 1).toUpperCase() + graphType.substring(1) + " Graphs - Time vs Complexity",
//                "Graph Size", "Time / Complexity (ms)",
//                dataset,
//                PlotOrientation.VERTICAL,
//                true, true, false);
//
//        ChartPanel panel = new ChartPanel(chart);
//        JFrame frame = new JFrame(graphType.substring(0, 1).toUpperCase() + graphType.substring(1) + " Graphs - Time vs Complexity");
//        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//        frame.add(panel);
//        frame.setSize(800, 600);
//        frame.setVisible(true);
//    }
//
//
//    private static void createAndShowBarChart(String graphType, int[] numVerticesArray, double[] averageTimes, double[] expandedNodes, double[] distances, double[] theoreticalComplexities) {
//        DefaultCategoryDataset dataset = new DefaultCategoryDataset();
//
//        for (int i = 0; i < numVerticesArray.length; i++) {
//            dataset.addValue(averageTimes[i], "Actual Time (ms)", "Vertices " + numVerticesArray[i] + " - " + graphType);
//            dataset.addValue(expandedNodes[i], "Expanded Nodes", "Vertices " + numVerticesArray[i] + " - " + graphType);
//            dataset.addValue(distances[i], "Distance", "Vertices " + numVerticesArray[i] + " - " + graphType);
//            dataset.addValue(theoreticalComplexities[i], "Theoretical Complexity", "Vertices " + numVerticesArray[i] + " - " + graphType);
//        }
//
//        JFreeChart chart = ChartFactory.createBarChart(
//                "Query Execution Analysis - " + graphType, // 图表标题
//                "Graph Size and Type", // 域轴标签
//                "Value", // 范围轴标签
//                dataset, // 数据集
//                PlotOrientation.VERTICAL, // 绘图方向
//                true, // 是否显示图例
//                true, // 是否生成工具
//                false // 是否生成URL链接
//        );
//
//        ChartPanel panel = new ChartPanel(chart);
//        JFrame frame = new JFrame("Query Execution Analysis - " + graphType);
//        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//        frame.add(panel);
//        frame.setSize(800, 600);
//        frame.setVisible(true);
//    }
    // 对数变换：避免负无穷，log(0+1) = 0
    private static double logTransform(double value) {
        return Math.log(value + 1); // +1 防止对数为负无穷
    }

    // 创建并显示线性图：时间 vs 理论复杂度
    private static void createAndShowLineChart(String graphType, int[] numVerticesArray, double[] averageTimes, double[] theoreticalComplexities) {
        XYSeries seriesTime = new XYSeries("Actual Time (ms)");
        XYSeries seriesComplexity = new XYSeries("Theoretical Complexity");

        for (int i = 0; i < numVerticesArray.length; i++) {
            // 对数据进行对数变换
            seriesTime.add(numVerticesArray[i], logTransform(averageTimes[i]));
            seriesComplexity.add(numVerticesArray[i], logTransform(theoreticalComplexities[i]));
        }

        XYSeriesCollection dataset = new XYSeriesCollection();
        dataset.addSeries(seriesTime);
        dataset.addSeries(seriesComplexity);

        // 创建图表
        JFreeChart chart = ChartFactory.createXYLineChart(
                graphType.substring(0, 1).toUpperCase() + graphType.substring(1) + " Graphs - Time vs Complexity",
                "Graph Size", "Log(Time / Complexity)",
                dataset,
                PlotOrientation.VERTICAL,
                true, true, false);

        // 显示图表
        ChartPanel panel = new ChartPanel(chart);
        JFrame frame = new JFrame(graphType.substring(0, 1).toUpperCase() + graphType.substring(1) + " Graphs - Time vs Complexity");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(panel);
        frame.setSize(800, 600);
        frame.setVisible(true);
    }

    // 创建并显示条形图：时间、扩展节点数、距离、理论复杂度
    private static void createAndShowBarChart(String graphType, int[] numVerticesArray, double[] averageTimes, double[] expandedNodes, double[] distances, double[] theoreticalComplexities) {
        DefaultCategoryDataset dataset = new DefaultCategoryDataset();

        for (int i = 0; i < numVerticesArray.length; i++) {
            // 对数据进行对数变换
            dataset.addValue(logTransform(averageTimes[i]), "Actual Time (ms)", "Vertices " + numVerticesArray[i] + " - " + graphType);
            dataset.addValue(logTransform(expandedNodes[i]), "Expanded Nodes", "Vertices " + numVerticesArray[i] + " - " + graphType);
            dataset.addValue(logTransform(distances[i]), "Distance", "Vertices " + numVerticesArray[i] + " - " + graphType);
            dataset.addValue(logTransform(theoreticalComplexities[i]), "Theoretical Complexity", "Vertices " + numVerticesArray[i] + " - " + graphType);
        }

        // 创建条形图
        JFreeChart chart = ChartFactory.createBarChart(
                "Query Execution Analysis - " + graphType, // 图表标题
                "Graph Size and Type", // 域轴标签
                "Log(Value)", // 范围轴标签
                dataset, // 数据集
                PlotOrientation.VERTICAL, // 绘图方向
                true, // 是否显示图例
                true, // 是否生成工具
                false // 是否生成URL链接
        );

        // 显示图表
        ChartPanel panel = new ChartPanel(chart);
        JFrame frame = new JFrame("Query Execution Analysis - " + graphType);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.add(panel);
        frame.setSize(800, 600);
        frame.setVisible(true);
    }


}
