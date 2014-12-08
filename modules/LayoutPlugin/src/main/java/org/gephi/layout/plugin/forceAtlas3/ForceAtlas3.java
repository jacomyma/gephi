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
package org.gephi.layout.plugin.forceAtlas3;

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
import org.gephi.layout.plugin.forceAtlas3.ForceFactory.AttractionForce;
import org.gephi.layout.plugin.forceAtlas3.ForceFactory.RepulsionForce;
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
public class ForceAtlas3 implements Layout {

    private GraphModel graphModel;
    private HierarchicalGraph graph;
    private final ForceAtlas3Builder layoutBuilder;
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

    public ForceAtlas3(ForceAtlas3Builder layoutBuilder) {
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
            if (n.getNodeData().getLayoutData() == null || !(n.getNodeData().getLayoutData() instanceof ForceAtlas3LayoutData)) {
                ForceAtlas3LayoutData nLayout = new ForceAtlas3LayoutData();
                n.getNodeData().setLayoutData(nLayout);
            }
            NodeData nData = n.getNodeData();
            ForceAtlas3LayoutData nLayout = nData.getLayoutData();
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
        // benchmark();

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
            if (n.getNodeData().getLayoutData() == null || !(n.getNodeData().getLayoutData() instanceof ForceAtlas3LayoutData)) {
                ForceAtlas3LayoutData nLayout = new ForceAtlas3LayoutData();
                n.getNodeData().setLayoutData(nLayout);
            }
            NodeData nData = n.getNodeData();
            ForceAtlas3LayoutData nLayout = nData.getLayoutData();
            nLayout.mass = 1 + graph.getDegree(n);
            nLayout.old_dx = nLayout.dx;
            nLayout.old_dy = nLayout.dy;
            nLayout.dx = 0;
            nLayout.dy = 0;
            nLayout.average_dx = 0;
            nLayout.average_dy = 0;
            nLayout.average_swinging = 0;
            nLayout.average_traction = 0;
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
                ForceAtlas3LayoutData nLayout = nData.getLayoutData();
                outboundAttCompensation += nLayout.mass;
            }
            outboundAttCompensation /= nodes.length;
        }

        // Repulsion (and gravity)
        // NB: Muti-threaded
        double linLogModifier = 1.;
        if(isLinLogMode()){
            linLogModifier = 0.05;
        }
        RepulsionForce Repulsion = ForceFactory.builder.buildRepulsion(isAdjustSizes(), getScalingRatio() * linLogModifier);

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
        
        // Average dx and dy
        for (Edge e : edges) {
            NodeData nSourceData = e.getSource().getNodeData();
            ForceAtlas3LayoutData nSourceLayout = nSourceData.getLayoutData();
            NodeData nTargetData = e.getTarget().getNodeData();
            ForceAtlas3LayoutData nTargetLayout = nTargetData.getLayoutData();
            
            nSourceLayout.average_dx += nTargetLayout.dx;
            nTargetLayout.average_dx += nSourceLayout.dx;
            
            nSourceLayout.average_dy += nTargetLayout.dy;
            nTargetLayout.average_dy += nSourceLayout.dy;
        }
        for(Node n : nodes){
            NodeData nData = n.getNodeData();
            ForceAtlas3LayoutData nLayout = nData.getLayoutData();
            nLayout.average_dx = (nLayout.average_dx) / nLayout.mass;
            nLayout.average_dy = (nLayout.average_dy) / nLayout.mass;
        }
        
        // Share displacement
        for(Node n : nodes){
            NodeData nData = n.getNodeData();
            ForceAtlas3LayoutData nLayout = nData.getLayoutData();

            double flocknessInertia = 0.5;
            nLayout.flockness = flocknessInertia * nLayout.flockness + (1 - flocknessInertia) * (nLayout.average_dx * nLayout.dy + nLayout.dx * nLayout.average_dy);

            double averageDisplacement = Math.sqrt(nLayout.average_dx * nLayout.average_dx + nLayout.average_dy * nLayout.average_dy);
            double ratio = 1;
            nLayout.dx += ratio * nLayout.average_dx;
            nLayout.dy += ratio * nLayout.average_dy;
            nData.setSize((float) (10f + Math.sqrt(averageDisplacement)));
        }

        // Compute and store swinging info
        for(Node n : nodes){
            NodeData nData = n.getNodeData();
            ForceAtlas3LayoutData nLayout = nData.getLayoutData();
            nLayout.swinging = Math.sqrt((nLayout.old_dx - nLayout.dx) * (nLayout.old_dx - nLayout.dx) + (nLayout.old_dy - nLayout.dy) * (nLayout.old_dy - nLayout.dy));
            nLayout.traction = Math.sqrt((nLayout.old_dx + nLayout.dx) * (nLayout.old_dx + nLayout.dx) + (nLayout.old_dy + nLayout.dy) * (nLayout.old_dy + nLayout.dy)) / 2;
            
        }
        
        // Average these
        for (Edge e : edges) {
            NodeData nSourceData = e.getSource().getNodeData();
            ForceAtlas3LayoutData nSourceLayout = nSourceData.getLayoutData();
            NodeData nTargetData = e.getTarget().getNodeData();
            ForceAtlas3LayoutData nTargetLayout = nTargetData.getLayoutData();
            
            nSourceLayout.average_swinging += nTargetLayout.swinging;
            nTargetLayout.average_swinging += nSourceLayout.swinging;
            
            nSourceLayout.average_traction += nTargetLayout.traction;
            nTargetLayout.average_traction += nSourceLayout.traction;
        }
        for(Node n : nodes){
            NodeData nData = n.getNodeData();
            ForceAtlas3LayoutData nLayout = nData.getLayoutData();
            nLayout.average_swinging = (nLayout.average_swinging + nLayout.swinging) / nLayout.mass;
            nLayout.average_traction = (nLayout.average_traction + nLayout.traction) / nLayout.mass;
        }
        
        // Apply forces
        if (isAdjustSizes()) {
            // If nodes overlap prevention is active, it's not possible to trust the swinging mesure.
            for (Node n : nodes) {
                NodeData nData = n.getNodeData();
                ForceAtlas3LayoutData nLayout = nData.getLayoutData();
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
                ForceAtlas3LayoutData nLayout = nData.getLayoutData();
                if (!nData.isFixed()) {
                                        
                    //double nodespeed = nLayout.convergenceEstimation * Math.log(1 + nLayout.traction) / (1 + Math.sqrt(nLayout.swinging));
                    //nLayout.convergenceEstimation = Math.min(1, Math.sqrt(nodespeed * (nLayout.dx * nLayout.dx + nLayout.dy * nLayout.dy)) / (1 + Math.sqrt(nLayout.swinging)));
                    
                    double nodespeed = nLayout.convergenceEstimation * Math.log(1 + nLayout.traction) / (1 + Math.sqrt(nLayout.mass * nLayout.swinging));
                    nLayout.convergenceEstimation = Math.min(1, Math.sqrt(nodespeed * (nLayout.dx * nLayout.dx + nLayout.dy * nLayout.dy)) / (1 + Math.sqrt(nLayout.mass * nLayout.swinging)));
                    
                    // nData.setSize(0.1f * Math.min((float)Math.exp(0.1 * nLayout.flockness), 200));
                    /*if(nLayout.flockness > 1){
                        nData.setColor(1, 0, 0);
                    } else {
                        nData.setColor(0.6f, 0.6f, 0.6f);
                    }*/
                    
                    double x = nData.x() + nLayout.dx * nodespeed;
                    double y = nData.y() + nLayout.dy * nodespeed;
                    
                    nData.setX((float) x);
                    nData.setY((float) y);
                }
            }
        }
        
        layout_step++;
        
        graph.readUnlockAll();
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
        final String FORCEAtlas3_TUNING = NbBundle.getMessage(getClass(), "ForceAtlas3.tuning");
        final String FORCEAtlas3_BEHAVIOR = NbBundle.getMessage(getClass(), "ForceAtlas3.behavior");
        final String FORCEAtlas3_PERFORMANCE = NbBundle.getMessage(getClass(), "ForceAtlas3.performance");
        final String FORCEAtlas3_THREADS = NbBundle.getMessage(getClass(), "ForceAtlas3.threads");

        try {
            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas3.scalingRatio.name"),
                    FORCEAtlas3_TUNING,
                    "ForceAtlas3.scalingRatio.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas3.scalingRatio.desc"),
                    "getScalingRatio", "setScalingRatio"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas3.strongGravityMode.name"),
                    FORCEAtlas3_TUNING,
                    "ForceAtlas3.strongGravityMode.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas3.strongGravityMode.desc"),
                    "isStrongGravityMode", "setStrongGravityMode"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas3.gravity.name"),
                    FORCEAtlas3_TUNING,
                    "ForceAtlas3.gravity.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas3.gravity.desc"),
                    "getGravity", "setGravity"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas3.distributedAttraction.name"),
                    FORCEAtlas3_BEHAVIOR,
                    "ForceAtlas3.distributedAttraction.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas3.distributedAttraction.desc"),
                    "isOutboundAttractionDistribution", "setOutboundAttractionDistribution"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas3.linLogMode.name"),
                    FORCEAtlas3_BEHAVIOR,
                    "ForceAtlas3.linLogMode.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas3.linLogMode.desc"),
                    "isLinLogMode", "setLinLogMode"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas3.adjustSizes.name"),
                    FORCEAtlas3_BEHAVIOR,
                    "ForceAtlas3.adjustSizes.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas3.adjustSizes.desc"),
                    "isAdjustSizes", "setAdjustSizes"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas3.edgeWeightInfluence.name"),
                    FORCEAtlas3_BEHAVIOR,
                    "ForceAtlas3.edgeWeightInfluence.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas3.edgeWeightInfluence.desc"),
                    "getEdgeWeightInfluence", "setEdgeWeightInfluence"));

            properties.add(LayoutProperty.createProperty(
                    this, Boolean.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas3.barnesHutOptimization.name"),
                    FORCEAtlas3_PERFORMANCE,
                    "ForceAtlas3.barnesHutOptimization.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas3.barnesHutOptimization.desc"),
                    "isBarnesHutOptimize", "setBarnesHutOptimize"));

            properties.add(LayoutProperty.createProperty(
                    this, Double.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas3.barnesHutTheta.name"),
                    FORCEAtlas3_PERFORMANCE,
                    "ForceAtlas3.barnesHutTheta.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas3.barnesHutTheta.desc"),
                    "getBarnesHutTheta", "setBarnesHutTheta"));

            properties.add(LayoutProperty.createProperty(
                    this, Integer.class,
                    NbBundle.getMessage(getClass(), "ForceAtlas3.threads.name"),
                    FORCEAtlas3_THREADS,
                    "ForceAtlas3.threads.name",
                    NbBundle.getMessage(getClass(), "ForceAtlas3.threads.desc"),
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
