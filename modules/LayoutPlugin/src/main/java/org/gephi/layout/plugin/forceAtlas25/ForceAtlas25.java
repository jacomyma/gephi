/*
Copyright 2008-2011 Gephi
Authors : Mathieu Jacomy <mathieu.jacomy@gmail.com>
Website : http://www.gephi.org

This file is part of Gephi.

DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS HEADER.

Copyright 2011 Gephi Consortium. All rights reserved.

The contents of this file are subject to the terms of either the GNU
General Public License Version 3 only ("GPL") or the Common
Development and Distribution License("CDDL") (collectively, the
"License"). You may not use this file except in compliance with the
License. You can obtain a copy of the License at
http://gephi.org/about/legal/license-notice/
or /cddl-1.0.txt and /gpl-3.0.txt. See the License for the
specific language governing permissions and limitations under the
License.  When distributing the software, include this License Header
Notice in each file and include the License files at
/cddl-1.0.txt and /gpl-3.0.txt. If applicable, add the following below the
License Header, with the fields enclosed by brackets [] replaced by
your own identifying information:
"Portions Copyrighted [year] [name of copyright owner]"

If you wish your version of this file to be governed by only the CDDL
or only the GPL Version 3, indicate your decision by adding
"[Contributor] elects to include this software in this distribution
under the [CDDL or GPL Version 3] license." If you do not indicate a
single choice of license, a recipient has the option to distribute
your version of this file under either the CDDL, the GPL Version 3 or
to extend the choice of license to its licensees as provided above.
However, if you add GPL Version 3 code and therefore, elected the GPL
Version 3 license, then the option applies only if the new code is
made subject to such option by the copyright holder.

Contributor(s):

Portions Copyrighted 2011 Gephi Consortium.
 */
package org.gephi.layout.plugin.forceAtlas25;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import org.gephi.data.attributes.type.TimeInterval;
import org.gephi.dynamic.DynamicUtilities;
import org.gephi.dynamic.api.DynamicController;
import org.gephi.dynamic.api.DynamicModel;
import org.gephi.graph.api.Edge;
import org.gephi.graph.api.GraphModel;
import org.gephi.graph.api.HierarchicalGraph;
import org.gephi.graph.api.Node;
import org.gephi.graph.api.NodeData;
import org.gephi.layout.plugin.forceAtlas25.ForceFactory.AttractionForce;
import org.gephi.layout.plugin.forceAtlas25.ForceFactory.RepulsionForce;
import org.gephi.layout.spi.Layout;
import org.gephi.layout.spi.LayoutBuilder;
import org.gephi.layout.spi.LayoutProperty;
import org.gephi.project.api.Workspace;
import org.openide.util.Exceptions;
import org.openide.util.Lookup;
import org.openide.util.NbBundle;

/**
 * ForceAtlas 2 Layout, manages each layout_step of the computations.
 * @author Mathieu Jacomy
 */
public class ForceAtlas25 implements Layout {

    private GraphModel graphModel;
    private HierarchicalGraph graph;
    private final ForceAtlas25Builder layoutBuilder;
    private DynamicModel dynamicModel;
    private double edgeWeightInfluence;
    private double scalingRatio;
    private double gravity;
    private boolean outboundAttractionDistribution;
    private boolean adjustSizes;
    private boolean barnesHutOptimize;
    private double barnesHutTheta;
    private boolean linLogMode;
    private boolean strongGravityMode;
    private int threadCount;
    private int currentThreadCount;
    private Region rootRegion;
    double outboundAttCompensation = 1;
    //Dynamic Weight
    private TimeInterval timeInterval;
    private ExecutorService pool;
    private boolean converged;
    
    private int layout_step;

    public ForceAtlas25(ForceAtlas25Builder layoutBuilder) {
        this.layoutBuilder = layoutBuilder;
        this.threadCount = Math.min(4, Math.max(1, Runtime.getRuntime().availableProcessors() - 1));
    }

    @Override
    public void initAlgo() {
        graph = graphModel.getHierarchicalGraphVisible();
        this.timeInterval = DynamicUtilities.getVisibleInterval(dynamicModel);

        graph.readLock();
        Node[] nodes = graph.getNodes().toArray();

        // Initialise layout data
        for (Node n : nodes) {
            if (n.getNodeData().getLayoutData() == null || !(n.getNodeData().getLayoutData() instanceof ForceAtlas25LayoutData)) {
                ForceAtlas25LayoutData nLayout = new ForceAtlas25LayoutData();
                n.getNodeData().setLayoutData(nLayout);
            }
            NodeData nData = n.getNodeData();
            ForceAtlas25LayoutData nLayout = nData.getLayoutData();
            nLayout.mass = 1 + graph.getDegree(n);
            nLayout.old_dx = 0;
            nLayout.old_dy = 0;
            nLayout.dx = 0;
            nLayout.dy = 0;
            nLayout.convergenceEstimation = 1;
        }

        pool = Executors.newFixedThreadPool(threadCount);
        currentThreadCount = threadCount;
        
        layout_step = 0;
        
        setConverged(false);

        // Benchmark
        benchmark();

    }

    @Override
    public void goAlgo() {
        // Initialize graph data
        if (graphModel == null) {
            return;
        }
        graph = graphModel.getHierarchicalGraphVisible();
        this.timeInterval = DynamicUtilities.getVisibleInterval(dynamicModel);

        graph.readLock();
        Node[] nodes = graph.getNodes().toArray();
        Edge[] edges = graph.getEdgesAndMetaEdges().toArray();

        // Initialise layout data
        for (Node n : nodes) {
            if (n.getNodeData().getLayoutData() == null || !(n.getNodeData().getLayoutData() instanceof ForceAtlas25LayoutData)) {
                ForceAtlas25LayoutData nLayout = new ForceAtlas25LayoutData();
                n.getNodeData().setLayoutData(nLayout);
            }
            NodeData nData = n.getNodeData();
            ForceAtlas25LayoutData nLayout = nData.getLayoutData();
            nLayout.mass = 1 + graph.getDegree(n);
            nLayout.old_dx = nLayout.dx;
            nLayout.old_dy = nLayout.dy;
            nLayout.dx = 0;
            nLayout.dy = 0;
        }

        // If Barnes Hut active, initialize root region
        if (isBarnesHutOptimize()) {
            rootRegion = new Region(nodes);
            rootRegion.buildSubRegions();
        }

        // If outboundAttractionDistribution active, compensate.
        if (isOutboundAttractionDistribution()) {
            outboundAttCompensation = 0;
            for (Node n : nodes) {
                NodeData nData = n.getNodeData();
                ForceAtlas25LayoutData nLayout = nData.getLayoutData();
                outboundAttCompensation += nLayout.mass;
            }
            outboundAttCompensation /= nodes.length;
        }

        // Repulsion (and gravity)
        // NB: Muti-threaded
        RepulsionForce Repulsion = ForceFactory.builder.buildRepulsion(isAdjustSizes(), getScalingRatio());

        int taskCount = 8 * currentThreadCount;  // The threadPool Executor Service will manage the fetching of tasks and threads.
        // We make more tasks than threads because some tasks may need more time to compute.
        ArrayList<Future> threads = new ArrayList();
        for (int t = taskCount; t > 0; t--) {
            int from = (int) Math.floor(nodes.length * (t - 1) / taskCount);
            int to = (int) Math.floor(nodes.length * t / taskCount);
            Future future = pool.submit(new NodesThread(nodes, from, to, isBarnesHutOptimize(), getBarnesHutTheta(), getGravity(), (isStrongGravityMode()) ? (ForceFactory.builder.getStrongGravity(getScalingRatio())) : (Repulsion), getScalingRatio(), rootRegion, Repulsion));
            threads.add(future);
        }
        for (Future future : threads) {
            try {
                future.get();
            } catch (InterruptedException ex) {
                Exceptions.printStackTrace(ex);
            } catch (ExecutionException ex) {
                Exceptions.printStackTrace(ex);
            }
        }

        // Attraction
        AttractionForce Attraction = ForceFactory.builder.buildAttraction(isLinLogMode(), isOutboundAttractionDistribution(), isAdjustSizes(), 1 * ((isOutboundAttractionDistribution()) ? (outboundAttCompensation) : (1)));
        if (getEdgeWeightInfluence() == 0) {
            for (Edge e : edges) {
                Attraction.apply(e.getSource(), e.getTarget(), 1);
            }
        } else if (getEdgeWeightInfluence() == 1) {
            for (Edge e : edges) {
                Attraction.apply(e.getSource(), e.getTarget(), getWeight(e));
            }
        } else {
            for (Edge e : edges) {
                Attraction.apply(e.getSource(), e.getTarget(), Math.pow(getWeight(e), getEdgeWeightInfluence()));
            }
        }
        
        // Apply forces
        if (isAdjustSizes()) {
            // If nodes overlap prevention is active, it's not possible to trust the swinging mesure.
            for (Node n : nodes) {
                NodeData nData = n.getNodeData();
                ForceAtlas25LayoutData nLayout = nData.getLayoutData();
                if (!nData.isFixed()) {
                    // Limit force
                    double force = Math.sqrt(nLayout.dx*nLayout.dx + nLayout.dy*nLayout.dy);
                    double maxForce = 10;
                    if(force > maxForce){
                        nLayout.dx = nLayout.dx * maxForce / force;
                        nLayout.dy = nLayout.dy * maxForce / force;
                    }
                    
                    // Adaptive auto-speed: the speed of each node is lowered
                    // when the node swings.
                    double swinging = nLayout.mass * Math.sqrt((nLayout.old_dx - nLayout.dx) * (nLayout.old_dx - nLayout.dx) + (nLayout.old_dy - nLayout.dy) * (nLayout.old_dy - nLayout.dy));
                    double traction = Math.sqrt((nLayout.old_dx + nLayout.dx) * (nLayout.old_dx + nLayout.dx) + (nLayout.old_dy + nLayout.dy) * (nLayout.old_dy + nLayout.dy)) / 2;
                    
                    double nodespeed = 0.1 * Math.log(1 + traction) / (1 + Math.sqrt(swinging));
                    
                    double x = nData.x() + nLayout.dx * nodespeed;
                    double y = nData.y() + nLayout.dy * nodespeed;
                    
                    nData.setX((float) x);
                    nData.setY((float) y);
                }
            }
        } else {
            for (Node n : nodes) {
                NodeData nData = n.getNodeData();
                ForceAtlas25LayoutData nLayout = nData.getLayoutData();
                if (!nData.isFixed()) {
                    
                    // Adaptive auto-speed: the speed of each node is lowered
                    // when the node swings.
                    double swinging = nLayout.mass * Math.sqrt((nLayout.old_dx - nLayout.dx) * (nLayout.old_dx - nLayout.dx) + (nLayout.old_dy - nLayout.dy) * (nLayout.old_dy - nLayout.dy));
                    double traction = Math.sqrt((nLayout.old_dx + nLayout.dx) * (nLayout.old_dx + nLayout.dx) + (nLayout.old_dy + nLayout.dy) * (nLayout.old_dy + nLayout.dy)) / 2;
                    
                    double nodespeed = nLayout.convergenceEstimation * Math.log(1 + traction) / (1 + Math.sqrt(swinging));
                    
                    nLayout.convergenceEstimation = Math.min(1, Math.sqrt(nodespeed * (nLayout.dx * nLayout.dx + nLayout.dy * nLayout.dy)) / (1 + Math.sqrt(swinging)));
                    
                    double x = nData.x() + nLayout.dx * nodespeed;
                    double y = nData.y() + nLayout.dy * nodespeed;
                    
                    nData.setX((float) x);
                    nData.setY((float) y);
                }
            }
        }
        
        layout_step++;
        
        // Benchmark
        benchmark();
        
        if(layout_step >= 100){
            setConverged(true);
        }
        
        graph.readUnlockAll();
    }
    
    public void benchmark(){
        Node[] nodes = graph.getNodes().toArray();
        Edge[] edges = graph.getEdgesAndMetaEdges().toArray();

        
        // We compute Noack's normalized^endv atedge length
        double card_e = graph.getEdgeCount();
        double card_n2 = graph.getNodeCount() * graph.getNodeCount();
        double sum_edges_distances = 0;
        for(Edge e : edges){
            NodeData sourceData = e.getSource().getNodeData();
            NodeData targetData = e.getTarget().getNodeData();
            double distance = Math.sqrt((sourceData.x() - targetData.x())*(sourceData.x() - targetData.x()) + (sourceData.y() - targetData.y())*(sourceData.y() - targetData.y()));
            sum_edges_distances += distance;
        }
        double sum_npairs_distances = 0;
        for (Node n1 : nodes) {
            NodeData nData1 = n1.getNodeData();
            for (Node n2 : nodes) {
                NodeData nData2 = n2.getNodeData();
                if(n1.getId() < n2.getId()){
                    double distance = Math.sqrt((nData1.x() - nData2.x())*(nData1.x() - nData2.x()) + (nData1.y() - nData2.y())*(nData1.y() - nData2.y()));
                    sum_npairs_distances += distance;
                }
            }
        }
        double neal = (sum_edges_distances / card_e) / (sum_npairs_distances / card_n2);

        // We compute the number of edge crossings
        // http://www.dcs.gla.ac.uk/publications/PAPERS/6621/final.pdf
        double c_all = (card_e * (card_e - 1)) / 2;
        double c_impossible = 0;
        for (Node n : nodes) {
            double degree = graph.getDegree(n);
            c_impossible += degree * (degree - 1);
        }
        c_impossible = c_impossible/2;
        double c_max = c_all - c_impossible;
        double aleph_c;
        if(c_max > 0){
            double c = 0;
            for(Edge e1 : edges){
                NodeData sourceData1 = e1.getSource().getNodeData();
                NodeData targetData1 = e1.getTarget().getNodeData();
                for(Edge e2 : edges){
                    if(e1.getId() < e2.getId()){
                        NodeData sourceData2 = e2.getSource().getNodeData();
                        NodeData targetData2 = e2.getTarget().getNodeData();
                        if(doLineSegmentsIntersect(sourceData1.x(), sourceData1.y(), targetData1.x(), targetData1.y(), sourceData2.x(), sourceData2.y(), targetData2.x(), targetData2.y())){
                            c += 1;
                        }
                    }
                }
            }
            aleph_c = 1 - c / c_max;
        } else {
            aleph_c = 0;
        }


        String layout_signature;
        if(linLogMode){
            layout_signature = "FA25_LL";
        } else {
            layout_signature = "FA25";
        }
        System.out.println("#benchmark,"+layout_signature+","+layout_step+","+neal + "," + (1-aleph_c));
    }
    
    public boolean doLineSegmentsIntersect(double px, double py, double p2x, double p2y, double qx, double qy, double q2x, double q2y){
        double rx = p2x - px;
        double ry = p2y - py;
        double sx = q2x - qx;
        double sy = q2y - qy;
        
        if(((px == qx) && (py == qy)) || ((px == q2x) && (py == q2y)) || ((p2x == qx) && (p2y == qy)) || ((p2x == q2x) && (p2y == q2y))){
            return false;
        }
        
        double sub_qp_x = qx - px;
        double sub_qp_y = qy - py;
        
        double uNumerator = sub_qp_x * ry - sub_qp_y * rx;
        double denominator = rx * sy - ry * sx;
        
        if(uNumerator == 0 && denominator == 0){
            // Colinear, so do they overlap ?
            return ((qx - px < 0) != (qx - p2x < 0) != (q2x - px < 0) != (q2x - p2x < 0)) || ((qy - py < 0) != (qy - p2y < 0) != (q2y - py < 0) != (q2y - p2y < 0));
        }
        
        if(denominator == 0){
            // Parallel
            return false;
        }
        
        double u = uNumerator / denominator;
        double t = (sub_qp_x * sy - sub_qp_y * sx) / denominator;
        
        return (t >= 0) && (t <= 1) && (u >= 0) && (u <= 1);
    }

    @Override
    public boolean canAlgo() {
        return !isConverged() && graphModel != null;
    }

    public void setConverged(boolean converged) {
        this.converged = converged;
    }

    public boolean isConverged() {
        return converged;
    }

    @Override
    public void endAlgo() {
        for (Node n : graph.getNodes()) {
            n.getNodeData().setLayoutData(null);
        }
        pool.shutdown();
        graph.readUnlockAll();
    }

    @Override
    public LayoutProperty[] getProperties() {
        List<LayoutProperty> properties = new ArrayList<LayoutProperty>();
        final String FORCEATLAS25_TUNING = NbBundle.getMessage(getClass(), "ForceAtlas25.tuning");
        final String FORCEATLAS25_BEHAVIOR = NbBundle.getMessage(getClass(), "ForceAtlas25.behavior");
        final String FORCEATLAS25_PERFORMANCE = NbBundle.getMessage(getClass(), "ForceAtlas25.performance");
        final String FORCEATLAS25_THREADS = NbBundle.getMessage(getClass(), "ForceAtlas25.threads");

        try {
            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas25.scalingRatio.name"),
                    FORCEATLAS25_TUNING,
                    "ForceAtlas25.scalingRatio.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas25.scalingRatio.desc"),
                    "getScalingRatio", "setScalingRatio"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas25.strongGravityMode.name"),
                    FORCEATLAS25_TUNING,
                    "ForceAtlas25.strongGravityMode.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas25.strongGravityMode.desc"),
                    "isStrongGravityMode", "setStrongGravityMode"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas25.gravity.name"),
                    FORCEATLAS25_TUNING,
                    "ForceAtlas25.gravity.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas25.gravity.desc"),
                    "getGravity", "setGravity"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas25.distributedAttraction.name"),
                    FORCEATLAS25_BEHAVIOR,
                    "ForceAtlas25.distributedAttraction.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas25.distributedAttraction.desc"),
                    "isOutboundAttractionDistribution", "setOutboundAttractionDistribution"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas25.linLogMode.name"),
                    FORCEATLAS25_BEHAVIOR,
                    "ForceAtlas25.linLogMode.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas25.linLogMode.desc"),
                    "isLinLogMode", "setLinLogMode"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas25.adjustSizes.name"),
                    FORCEATLAS25_BEHAVIOR,
                    "ForceAtlas25.adjustSizes.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas25.adjustSizes.desc"),
                    "isAdjustSizes", "setAdjustSizes"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas25.edgeWeightInfluence.name"),
                    FORCEATLAS25_BEHAVIOR,
                    "ForceAtlas25.edgeWeightInfluence.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas25.edgeWeightInfluence.desc"),
                    "getEdgeWeightInfluence", "setEdgeWeightInfluence"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas25.barnesHutOptimization.name"),
                    FORCEATLAS25_PERFORMANCE,
                    "ForceAtlas25.barnesHutOptimization.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas25.barnesHutOptimization.desc"),
                    "isBarnesHutOptimize", "setBarnesHutOptimize"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas25.barnesHutTheta.name"),
                    FORCEATLAS25_PERFORMANCE,
                    "ForceAtlas25.barnesHutTheta.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas25.barnesHutTheta.desc"),
                    "getBarnesHutTheta", "setBarnesHutTheta"));

            properties.add(LayoutProperty.createProperty(
                    this, Integer.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas25.threads.name"),
                    FORCEATLAS25_THREADS,
                    "ForceAtlas25.threads.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas25.threads.desc"),
                    "getThreadsCount", "setThreadsCount"));

        } catch (Exception e) {
            e.printStackTrace();
        }

        return properties.toArray(new LayoutProperty[0]);
    }

    @Override
    public void resetPropertiesValues() {
        int nodesCount = 0;

        if (graphModel != null) {
            nodesCount = graphModel.getHierarchicalGraphVisible().getNodeCount();
        }

        // Tuning
        setScalingRatio(2.0);
        setStrongGravityMode(false);
        setGravity(1.);

        // Behavior
        setOutboundAttractionDistribution(false);
        setLinLogMode(false);
        setAdjustSizes(false);
        setEdgeWeightInfluence(0.);

        // Performance
        if (nodesCount >= 1000) {
            setBarnesHutOptimize(true);
        } else {
            setBarnesHutOptimize(false);
        }
        setBarnesHutTheta(1.2);
        setThreadsCount(4);
    }

    @Override
    public LayoutBuilder getBuilder() {
        return layoutBuilder;
    }

    @Override
    public void setGraphModel(GraphModel graphModel) {
        this.graphModel = graphModel;
        Workspace workspace = graphModel.getWorkspace();
        DynamicController dynamicController = Lookup.getDefault().lookup(DynamicController.class);
        if (dynamicController != null && workspace != null) {
            dynamicModel = dynamicController.getModel(workspace);
        }
        // Trick: reset here to take the profile of the graph in account for default values
        resetPropertiesValues();
    }

    public Double getBarnesHutTheta() {
        return barnesHutTheta;
    }

    public void setBarnesHutTheta(Double barnesHutTheta) {
        this.barnesHutTheta = barnesHutTheta;
    }

    public Double getEdgeWeightInfluence() {
        return edgeWeightInfluence;
    }

    public void setEdgeWeightInfluence(Double edgeWeightInfluence) {
        this.edgeWeightInfluence = edgeWeightInfluence;
    }

    public Boolean isLinLogMode() {
        return linLogMode;
    }

    public void setLinLogMode(Boolean linLogMode) {
        this.linLogMode = linLogMode;
    }

    public Double getScalingRatio() {
        return scalingRatio;
    }

    public void setScalingRatio(Double scalingRatio) {
        this.scalingRatio = scalingRatio;
    }

    public Boolean isStrongGravityMode() {
        return strongGravityMode;
    }

    public void setStrongGravityMode(Boolean strongGravityMode) {
        this.strongGravityMode = strongGravityMode;
    }

    public Double getGravity() {
        return gravity;
    }

    public void setGravity(Double gravity) {
        this.gravity = gravity;
    }

    public Integer getThreadsCount() {
        return threadCount;
    }

    public void setThreadsCount(Integer threadCount) {
        if (threadCount < 1) {
            setThreadsCount(1);
        } else {
            this.threadCount = threadCount;
        }
    }

    public Boolean isOutboundAttractionDistribution() {
        return outboundAttractionDistribution;
    }

    public void setOutboundAttractionDistribution(Boolean outboundAttractionDistribution) {
        this.outboundAttractionDistribution = outboundAttractionDistribution;
    }

    public Boolean isAdjustSizes() {
        return adjustSizes;
    }

    public void setAdjustSizes(Boolean adjustSizes) {
        this.adjustSizes = adjustSizes;
    }

    public Boolean isBarnesHutOptimize() {
        return barnesHutOptimize;
    }

    public void setBarnesHutOptimize(Boolean barnesHutOptimize) {
        this.barnesHutOptimize = barnesHutOptimize;
    }

    private float getWeight(Edge edge) {
        if (timeInterval != null) {
            return edge.getWeight(timeInterval.getLow(), timeInterval.getHigh());
        } else {
            return edge.getWeight();
        }
    }
}
