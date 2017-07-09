//   RTree.java
//   Java Spatial Index Library
//   Copyright (C) 2002-2005 Infomatiq Limited
//   Copyright (C) 2008-2010 aled@users.sourceforge.net
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

package net.sf.jsi.rtree;

import gnu.trove.list.array.TIntArrayList;
import gnu.trove.map.hash.TIntObjectHashMap;
import gnu.trove.procedure.TIntProcedure;
import gnu.trove.stack.TIntStack;
import gnu.trove.stack.array.TIntArrayStack;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Properties;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import net.sf.jsi.BuildProperties;
import net.sf.jsi.Point;
import net.sf.jsi.Rectangle;
import net.sf.jsi.PriorityQueue;
import net.sf.jsi.SpatialIndex;

/**
 * <p>
 * This is a lightweight RTree implementation, specifically designed for the
 * following features (in order of importance):
 * <ul>
 * <li>Fast intersection query performance. To achieve this, the RTree uses only
 * main memory to store entries. Obviously this will only improve performance if
 * there is enough physical memory to avoid paging.</li>
 * <li>Low memory requirements.</li>
 * <li>Fast add performance.</li>
 * </ul>
 * </p>
 *
 * <p>
 * The main reason for the high speed of this RTree implementation is the
 * avoidance of the creation of unnecessary objects, mainly achieved by using
 * primitive collections from the trove4j library.
 * </p>
 */
public class RTree implements SpatialIndex, Serializable {
	private static final long serialVersionUID = 5946232781609920309L;
	private static final Logger log = LoggerFactory.getLogger(RTree.class);
	private static final Logger deleteLog = LoggerFactory.getLogger(RTree.class.getName() + "-delete");

	// parameters of the tree
	private final static int DEFAULT_MAX_NODE_ENTRIES = 50;
	private final static int DEFAULT_MIN_NODE_ENTRIES = 20;
	int maxNodeEntries;
	int minNodeEntries;

	// map of nodeId -> node object
	// TODO eliminate this map - it should not be needed. Nodes
	// can be found by traversing the tree.
	private TIntObjectHashMap<Node> nodeMap = new TIntObjectHashMap<Node>();

	// internal consistency checking - set to true if debugging tree corruption
	private final static boolean INTERNAL_CONSISTENCY_CHECKING = false;

	// used to mark the status of entries during a node split
	private final static int ENTRY_STATUS_ASSIGNED = 0;
	private final static int ENTRY_STATUS_UNASSIGNED = 1;
	private byte[] entryStatus = null;
	private byte[] initialEntryStatus = null;

	// stacks used to store nodeId and entry index of each node
	// from the root down to the leaf. Enables fast lookup
	// of nodes when a split is propagated up the tree.
	private TIntStack parents = new TIntArrayStack();
	private TIntStack parentsEntry = new TIntArrayStack();

	// initialisation
	private int treeHeight = 1; // leaves are always level 1
	private int rootNodeId = 0;
	private int size = 0;

	// Enables creation of new nodes
	private int highestUsedNodeId = rootNodeId;

	// Deleted node objects are retained in the nodeMap,
	// so that they can be reused. Store the IDs of nodes
	// which can be reused.
	private TIntStack deletedNodeIds = new TIntArrayStack();
	private int dim;

	/**
	 * Constructor. Use init() method to initialize parameters of the RTree.
	 */
	public RTree() {
		return; // NOP
	}

	// -------------------------------------------------------------------------
	// public implementation of SpatialIndex interface:
	// init(Properties)
	// add(Rectangle, int)
	// delete(Rectangle, int)
	// nearest(Point, TIntProcedure, float)
	// intersects(Rectangle, TIntProcedure)
	// contains(Rectangle, TIntProcedure)
	// size()
	// -------------------------------------------------------------------------
	/**
	 * <p>
	 * Initialize implementation dependent properties of the RTree. Currently
	 * implemented properties are:
	 * <ul>
	 * <li>MaxNodeEntries</li> This specifies the maximum number of entries in a
	 * node. The default value is 10, which is used if the property is not
	 * specified, or is less than 2.
	 * <li>MinNodeEntries</li> This specifies the minimum number of entries in a
	 * node. The default value is half of the MaxNodeEntries value (rounded
	 * down), which is used if the property is not specified or is less than 1.
	 * </ul>
	 * </p>
	 *
	 * @see net.sf.jsi.SpatialIndex#init(Properties)
	 */
	public void init(Properties props) {
		if (props == null) {
			// use sensible defaults if null is passed in.
			maxNodeEntries = DEFAULT_MAX_NODE_ENTRIES;
			minNodeEntries = DEFAULT_MIN_NODE_ENTRIES;
		} else {
			maxNodeEntries = Integer.parseInt(props.getProperty("MaxNodeEntries", "0"));
			minNodeEntries = Integer.parseInt(props.getProperty("MinNodeEntries", "0"));
			dim = Integer.parseInt(props.getProperty("dim"));

			// Obviously a node with less than 2 entries cannot be split.
			// The node splitting algorithm will work with only 2 entries
			// per node, but will be inefficient.
			if (maxNodeEntries < 2) {
				log.warn("Invalid MaxNodeEntries = " + maxNodeEntries + " Resetting to default value of "
						+ DEFAULT_MAX_NODE_ENTRIES);
				maxNodeEntries = DEFAULT_MAX_NODE_ENTRIES;
			}

			// The MinNodeEntries must be less than or equal to (int)
			// (MaxNodeEntries / 2)
			if (minNodeEntries < 1 || minNodeEntries > maxNodeEntries / 2) {
				log.warn("MinNodeEntries must be between 1 and MaxNodeEntries / 2");
				minNodeEntries = maxNodeEntries / 2;
			}
		}

		entryStatus = new byte[maxNodeEntries];
		initialEntryStatus = new byte[maxNodeEntries];

		for (int i = 0; i < maxNodeEntries; i++) {
			initialEntryStatus[i] = ENTRY_STATUS_UNASSIGNED;
		}

		Node root = new Node(rootNodeId, 1, maxNodeEntries, dim);
		nodeMap.put(rootNodeId, root);

		log.debug("init() " + " MaxNodeEntries = " + maxNodeEntries + ", MinNodeEntries = " + minNodeEntries);
	}

	/**
	 * @see net.sf.jsi.SpatialIndex#add(Rectangle, int)
	 */
	public void add(Rectangle r, int id) {
		if (log.isDebugEnabled()) {
			log.debug("Adding rectangle " + r + ", id " + id);
		}

		add(r.minCoords, r.maxCoords, id, 1);

		size++;

		if (INTERNAL_CONSISTENCY_CHECKING) {
			checkConsistency();
		}
	}

	/**
	 * Adds a new entry at a specified level in the tree
	 */
	private void add(float[] minCoords, float[] maxCoords, int id, int level) {
		// I1 [Find position for new record] Invoke ChooseLeaf to select a
		// leaf node L in which to place r
		Node n = chooseNode(minCoords, maxCoords, level);
		Node newLeaf = null;

		// I2 [Add record to leaf node] If L has room for another entry,
		// install E. Otherwise invoke SplitNode to obtain L and LL containing
		// E and all the old entries of L
		if (n.entryCount < maxNodeEntries) {
			n.addEntry(minCoords.clone(), maxCoords.clone(), id);
		} else {
			newLeaf = splitNode(n, minCoords.clone(), maxCoords.clone(), id);
		}

		// I3 [Propagate changes upwards] Invoke AdjustTree on L, also passing
		// LL
		// if a split was performed
		Node newNode = adjustTree(n, newLeaf);

		// I4 [Grow tree taller] If node split propagation caused the root to
		// split, create a new root whose children are the two resulting nodes.
		if (newNode != null) {
			int oldRootNodeId = rootNodeId;
			Node oldRoot = getNode(oldRootNodeId);

			rootNodeId = getNextNodeId();
			treeHeight++;
			Node root = new Node(rootNodeId, treeHeight, maxNodeEntries, dim);
			root.addEntry(newNode.mbb.minCoords.clone(), newNode.mbb.maxCoords.clone(), newNode.nodeId);
			root.addEntry(oldRoot.mbb.minCoords.clone(), oldRoot.mbb.maxCoords.clone(), oldRoot.nodeId);
			nodeMap.put(rootNodeId, root);
		}
	}

	/**
	 * @see net.sf.jsi.SpatialIndex#delete(Rectangle, int)
	 */
	public boolean delete(Rectangle r, int id) {
		// FindLeaf algorithm inlined here. Note the "official" algorithm
		// searches all overlapping entries. This seems inefficient to me,
		// as an entry is only worth searching if it contains (NOT overlaps)
		// the rectangle we are searching for.
		//
		// Also the algorithm has been changed so that it is not recursive.

		// FL1 [Search subtrees] If root is not a leaf, check each entry
		// to determine if it contains r. For each entry found, invoke
		// findLeaf on the node pointed to by the entry, until r is found or
		// all entries have been checked.
		parents.clear();
		parents.push(rootNodeId);

		parentsEntry.clear();
		parentsEntry.push(-1);
		Node n = null;
		int foundIndex = -1; // index of entry to be deleted in leaf

		while (foundIndex == -1 && parents.size() > 0) {
			n = getNode(parents.peek());
			int startIndex = parentsEntry.peek() + 1;

			if (!n.isLeaf()) {
				deleteLog.debug("searching node " + n.nodeId + ", from index " + startIndex);
				boolean contains = false;
				for (int i = startIndex; i < n.entryCount; i++) {
					Rectangle rn = new Rectangle(n.entriesMin[i], n.entriesMax[i]);
					if (rn.contains(r)) {
						parents.push(n.ids[i]);
						parentsEntry.pop();
						parentsEntry.push(i); // this becomes the start index
												// when the child has been
												// searched
						parentsEntry.push(-1);
						contains = true;
						break; // ie go to next iteration of while()
					}
				}
				if (contains) {
					continue;
				}
			} else {
				foundIndex = n.findEntry(r.minCoords, r.maxCoords, id);
			}

			parents.pop();
			parentsEntry.pop();
		} // while not found

		if (foundIndex != -1 && n != null) {
			n.deleteEntry(foundIndex);
			condenseTree(n);
			size--;
		}

		// shrink the tree if possible (i.e. if root node has exactly one
		// entry,and that
		// entry is not a leaf node, delete the root (it's entry becomes the new
		// root)
		Node root = getNode(rootNodeId);
		while (root.entryCount == 1 && treeHeight > 1) {
			deletedNodeIds.push(rootNodeId);
			root.entryCount = 0;
			rootNodeId = root.ids[0];
			treeHeight--;
			root = getNode(rootNodeId);
		}

		// if the tree is now empty, then set the MBR of the root node back to
		// it's original state
		// (this is only needed when the tree is empty, as this is the only
		// state where an empty node
		// is not eliminated)
		if (size == 0) {
			Arrays.fill(root.mbb.minCoords, Float.POSITIVE_INFINITY);
			Arrays.fill(root.mbb.maxCoords, Float.NEGATIVE_INFINITY);
		}

		if (INTERNAL_CONSISTENCY_CHECKING) {
			checkConsistency();
		}

		return (foundIndex != -1);
	}

	/**
	 * @see net.sf.jsi.SpatialIndex#nearest(Point, TIntProcedure, float)
	 */
	public void nearest(Point p, TIntProcedure v, float furthestDistance) {
		Node rootNode = getNode(rootNodeId);

		float furthestDistanceSq = furthestDistance * furthestDistance;
		TIntArrayList nearestIds = new TIntArrayList();
		nearest(p, rootNode, furthestDistanceSq, nearestIds);

		nearestIds.forEach(v);
		nearestIds.reset();
	}

	private void createNearestNDistanceQueue(Point p, int count, PriorityQueue distanceQueue, float furthestDistance) {
		// return immediately if given an invalid "count" parameter
		if (count <= 0) {
			return;
		}

		TIntStack parents = new TIntArrayStack();
		parents.push(rootNodeId);

		TIntStack parentsEntry = new TIntArrayStack();
		parentsEntry.push(-1);

		TIntArrayList savedValues = new TIntArrayList();
		float savedPriority = 0;

		// TODO: possible shortcut here - could test for intersection with the
		// MBR of the root node. If no intersection, return immediately.

		float furthestDistanceSq = furthestDistance * furthestDistance;

		while (parents.size() > 0) {
			Node n = getNode(parents.peek());
			int startIndex = parentsEntry.peek() + 1;

			if (!n.isLeaf()) {
				// go through every entry in the index node to check
				// if it could contain an entry closer than the farthest entry
				// currently stored.
				boolean near = false;
				for (int i = startIndex; i < n.entryCount; i++) {
					if (Rectangle.distanceSq(n.entriesMin[i], n.entriesMax[i], p) <= furthestDistanceSq) {
						parents.push(n.ids[i]);
						parentsEntry.pop();
						parentsEntry.push(i); // this becomes the start index
												// when the child has been
												// searched
						parentsEntry.push(-1);
						near = true;
						break; // ie go to next iteration of while()
					}
				}
				if (near) {
					continue;
				}
			} else {
				// go through every entry in the leaf to check if
				// it is currently one of the nearest N entries.
				for (int i = 0; i < n.entryCount; i++) {
					float entryDistanceSq = Rectangle.distanceSq(n.entriesMin[i], n.entriesMax[i], p);
					int entryId = n.ids[i];

					if (entryDistanceSq <= furthestDistanceSq) {
						distanceQueue.insert(entryId, entryDistanceSq);

						while (distanceQueue.size() > count) {
							// normal case - we can simply remove the lowest
							// priority (highest distance) entry
							int value = distanceQueue.getValue();
							float distanceSq = distanceQueue.getPriority();
							distanceQueue.pop();

							// rare case - multiple items of the same priority
							// (distance)
							if (distanceSq == distanceQueue.getPriority()) {
								savedValues.add(value);
								savedPriority = distanceSq;
							} else {
								savedValues.reset();
							}
						}

						// if the saved values have the same distance as the
						// next one in the tree, add them back in.
						if (savedValues.size() > 0 && savedPriority == distanceQueue.getPriority()) {
							for (int svi = 0; svi < savedValues.size(); svi++) {
								distanceQueue.insert(savedValues.get(svi), savedPriority);
							}
							savedValues.reset();
						}

						// narrow the search, if we have already found N items
						if (distanceQueue.getPriority() < furthestDistanceSq && distanceQueue.size() >= count) {
							furthestDistanceSq = distanceQueue.getPriority();
						}
					}
				}
			}
			parents.pop();
			parentsEntry.pop();
		}
	}

	/**
	 * @see net.sf.jsi.SpatialIndex#nearestNUnsorted(Point, TIntProcedure, int,
	 *      float)
	 */
	public void nearestNUnsorted(Point p, TIntProcedure v, int count, float furthestDistance) {
		// This implementation is designed to give good performance
		// where
		// o N is high (100+)
		// o The results do not need to be sorted by distance.
		//
		// Uses a priority queue as the underlying data structure.
		//
		// Note that more than N items will be returned if items N and N+x have
		// the
		// same priority.
		PriorityQueue distanceQueue = new PriorityQueue(PriorityQueue.SORT_ORDER_DESCENDING);
		createNearestNDistanceQueue(p, count, distanceQueue, furthestDistance);

		while (distanceQueue.size() > 0) {
			v.execute(distanceQueue.getValue());
			distanceQueue.pop();
		}
	}

	/**
	 * @see net.sf.jsi.SpatialIndex#nearestN(Point, TIntProcedure, int, float)
	 */
	public void nearestN(Point p, TIntProcedure v, int count, float furthestDistance) {
		PriorityQueue distanceQueue = new PriorityQueue(PriorityQueue.SORT_ORDER_DESCENDING);
		createNearestNDistanceQueue(p, count, distanceQueue, furthestDistance);
		distanceQueue.setSortOrder(PriorityQueue.SORT_ORDER_ASCENDING);

		while (distanceQueue.size() > 0) {
			v.execute(distanceQueue.getValue());
			distanceQueue.pop();
		}
	}

	/**
	 * @see net.sf.jsi.SpatialIndex#intersects(Rectangle, TIntProcedure)
	 */
	public void intersects(Rectangle r, TIntProcedure v) {
		Node rootNode = getNode(rootNodeId);
		intersects(r, v, rootNode);
	}

	/**
	 * @see net.sf.jsi.SpatialIndex#contains(Rectangle, TIntProcedure)
	 */
	public void contains(Rectangle r, TIntProcedure v) {
		// find all rectangles in the tree that are contained by the passed
		// rectangle
		// written to be non-recursive (should model other searches on this?)
		TIntStack parents = new TIntArrayStack();
		parents.push(rootNodeId);

		TIntStack parentsEntry = new TIntArrayStack();
		parentsEntry.push(-1);

		// TODO: possible shortcut here - could test for intersection with the
		// MBR of the root node. If no intersection, return immediately.

		while (parents.size() > 0) {
			Node n = getNode(parents.peek());
			int startIndex = parentsEntry.peek() + 1;

			if (!n.isLeaf()) {
				// go through every entry in the index node to check
				// if it intersects the passed rectangle. If so, it
				// could contain entries that are contained.
				boolean intersects = false;
				for (int i = startIndex; i < n.entryCount; i++) {
					Rectangle rn = new Rectangle(n.entriesMin[i].clone(), n.entriesMax[i].clone());
					if (r.intersects(rn)) {
						parents.push(n.ids[i]);
						parentsEntry.pop();
						parentsEntry.push(i); // this becomes the start index
												// when the child has been
												// searched
						parentsEntry.push(-1);
						intersects = true;
						break; // ie go to next iteration of while()
					}
				}
				if (intersects) {
					continue;
				}
			} else {
				// go through every entry in the leaf to check if
				// it is contained by the passed rectangle
				for (int i = 0; i < n.entryCount; i++) {
					Rectangle rn = new Rectangle(n.entriesMin[i].clone(), n.entriesMax[i].clone());
					if (r.contains(rn)) {
						if (!v.execute(n.ids[i])) {
							return;
						}
					}
				}
			}
			parents.pop();
			parentsEntry.pop();
		}
	}

	/**
	 * @see net.sf.jsi.SpatialIndex#size()
	 */
	public int size() {
		return size;
	}

	/**
	 * @see net.sf.jsi.SpatialIndex#getBounds()
	 */
	public Rectangle getBounds() {
		Rectangle bounds = null;

		Node n = getNode(getRootNodeId());
		if (n != null && n.entryCount > 0) {
			bounds = new Rectangle(dim);
			bounds.minCoords = n.mbb.minCoords.clone();
			bounds.maxCoords = n.mbb.maxCoords.clone();
		}
		return bounds;
	}

	/**
	 * @see net.sf.jsi.SpatialIndex#getVersion()
	 */
	public String getVersion() {
		return "RTree-" + BuildProperties.getVersion();
	}
	// -------------------------------------------------------------------------
	// end of SpatialIndex methods
	// -------------------------------------------------------------------------

	/**
	 * Get the next available node ID. Reuse deleted node IDs if possible
	 */
	private int getNextNodeId() {
		int nextNodeId = 0;
		if (deletedNodeIds.size() > 0) {
			nextNodeId = deletedNodeIds.pop();
		} else {
			nextNodeId = 1 + highestUsedNodeId++;
		}
		return nextNodeId;
	}

	/**
	 * Get a node object, given the ID of the node.
	 */
	public Node getNode(int id) {
		return nodeMap.get(id);
	}

	/**
	 * Get the highest used node ID
	 */
	public int getHighestUsedNodeId() {
		return highestUsedNodeId;
	}

	/**
	 * Get the root node ID
	 */
	public int getRootNodeId() {
		return rootNodeId;
	}

	/**
	 * Split a node. Algorithm is taken pretty much verbatim from Guttman's
	 * original paper.
	 *
	 * @return new node object.
	 */
	private Node splitNode(Node n, float[] newRectMin, float[] newRectMax, int newId) {
		// [Pick first entry for each group] Apply algorithm pickSeeds to
		// choose two entries to be the first elements of the groups. Assign
		// each to a group.

		// debug code
		float initialArea = 0;
		// if (log.isDebugEnabled()) {
		// float unionMinX = Math.min(n.mbrMinX, newRectMinX);
		// float unionMinY = Math.min(n.mbrMinY, newRectMinY);
		// float unionMaxX = Math.max(n.mbrMaxX, newRectMaxX);
		// float unionMaxY = Math.max(n.mbrMaxY, newRectMaxY);
		//
		// initialArea = (unionMaxX - unionMinX) * (unionMaxY - unionMinY);
		// }

		System.arraycopy(initialEntryStatus, 0, entryStatus, 0, maxNodeEntries);

		Node newNode = null;
		newNode = new Node(getNextNodeId(), n.level, maxNodeEntries, dim);
		nodeMap.put(newNode.nodeId, newNode);

		pickSeeds(n, newRectMin, newRectMax, newId, newNode); // this also sets
																// the
																// entryCount to
																// 1

		// [Check if done] If all entries have been assigned, stop. If one
		// group has so few entries that all the rest must be assigned to it in
		// order for it to have the minimum number m, assign them and stop.
		while (n.entryCount + newNode.entryCount < maxNodeEntries + 1) {
			if (maxNodeEntries + 1 - newNode.entryCount == minNodeEntries) {
				// assign all remaining entries to original node
				for (int i = 0; i < maxNodeEntries; i++) {
					if (entryStatus[i] == ENTRY_STATUS_UNASSIGNED) {
						entryStatus[i] = ENTRY_STATUS_ASSIGNED;

						for (int e = 0; e < dim; e++) {
							if (n.entriesMin[i][e] < n.mbb.minCoords[e])
								n.mbb.minCoords[e] = n.entriesMin[i][e];
							if (n.entriesMax[i][e] > n.mbb.maxCoords[e])
								n.mbb.maxCoords[e] = n.entriesMax[i][e];
						}

						n.entryCount++;
					}
				}
				break;
			}
			if (maxNodeEntries + 1 - n.entryCount == minNodeEntries) {
				// assign all remaining entries to new node
				for (int i = 0; i < maxNodeEntries; i++) {
					if (entryStatus[i] == ENTRY_STATUS_UNASSIGNED) {
						entryStatus[i] = ENTRY_STATUS_ASSIGNED;
						newNode.addEntry(n.entriesMin[i], n.entriesMax[i], n.ids[i]);
						n.ids[i] = -1; // an id of -1 indicates the entry is not
										// in use
					}
				}
				break;
			}

			// [Select entry to assign] Invoke algorithm pickNext to choose the
			// next entry to assign. Add it to the group whose covering
			// rectangle
			// will have to be enlarged least to accommodate it. Resolve ties
			// by adding the entry to the group with smaller area, then to the
			// the one with fewer entries, then to either. Repeat from S2
			pickNext(n, newNode);
		}

		n.reorganize(this);

		// check that the MBR stored for each node is correct.
		if (INTERNAL_CONSISTENCY_CHECKING) {
			Rectangle nMBR = new Rectangle(n.mbb.minCoords.clone(), n.mbb.maxCoords.clone());
			if (!nMBR.equals(calculateMBR(n))) {
				log.error("Error: splitNode old node MBR wrong");
			}
			Rectangle newNodeMBR = new Rectangle(newNode.mbb.minCoords.clone(), newNode.mbb.maxCoords.clone());
			if (!newNodeMBR.equals(calculateMBR(newNode))) {
				log.error("Error: splitNode new node MBR wrong");
			}
		}

		// debug code
		if (log.isDebugEnabled()) {
			float newArea = Rectangle.area(n.mbb) + Rectangle.area(newNode.mbb);
			float percentageIncrease = (100 * (newArea - initialArea)) / initialArea;
			log.debug("Node " + n.nodeId + " split. New area increased by " + percentageIncrease + "%");
		}

		return newNode;
	}

	/**
	 * Pick the seeds used to split a node. Select two entries to be the first
	 * elements of the groups
	 */
	private void pickSeeds(Node n, float[] newRectMin, float[] newRectMax, int newId, Node newNode) {
		// Find extreme rectangles along all dimension. Along each dimension,
		// find the entry whose rectangle has the highest low side, and the one
		// with the lowest high side. Record the separation.
		float maxNormalizedSeparation = -1; // initialize to -1 so that even
											// overlapping rectangles will be
											// considered for the seeds
		int highestLowIndex = -1;
		int lowestHighIndex = -1;

		// for the purposes of picking seeds, take the MBR of the node to
		// include
		// the new rectangle aswell.
		for (int e = 0; e < dim; e++) {
			if (newRectMin[e] < n.mbb.minCoords[e])
				n.mbb.minCoords[e] = newRectMin[e];
			if (newRectMax[e] > n.mbb.maxCoords[e])
				n.mbb.maxCoords[e] = newRectMax[e];
		}

		float[] mbrLen = new float[dim];
		for (int e = 0; e < dim; e++) {
			mbrLen[e] = n.mbb.maxCoords[e] - n.mbb.minCoords[e];
		}

		if (log.isDebugEnabled()) {
			log.debug("pickSeeds(): NodeId = " + n.nodeId);
		}

		for (int e = 0; e < dim; e++) {

			float tempHighestLow = newRectMin[e];
			int tempHighestLowIndex = -1; // -1 indicates the new rectangle is
											// the seed

			float tempLowestHigh = newRectMax[e];
			int tempLowestHighIndex = -1; // -1 indicates the new rectangle is
											// the seed
			
			for (int i = 0; i < n.entryCount; i++) {
				float tempLow = n.entriesMin[i][e];
				if (tempLow >= tempHighestLow) {
					tempHighestLow = tempLow;
					tempHighestLowIndex = i;
				} else { // ensure that the same index cannot be both lowestHigh
							// and highestLow
					float tempHigh = n.entriesMax[i][e];
					if (tempHigh <= tempLowestHigh) {
						tempLowestHigh = tempHigh;
						tempLowestHighIndex = i;
					}
				}

				// PS2 [Adjust for shape of the rectangle cluster] Normalize the
				// separations
				// by dividing by the widths of the entire set along the
				// corresponding
				// dimension
				float normalizedSeparation = mbrLen[e] == 0 ? 1 : (tempHighestLow - tempLowestHigh) / mbrLen[e];
				if (normalizedSeparation > 1.0f || normalizedSeparation < -1.0f) {
					log.error("Invalid normalized separation " + e);
				}

				if (log.isDebugEnabled()) {
					log.debug("Entry " + i + ", dimension " + e + ": HighestLow = " + tempHighestLow + " (index "
							+ tempHighestLowIndex + ")" + ", LowestHigh = " + tempLowestHigh + " (index "
							+ tempLowestHighIndex + ", NormalizedSeparation = " + normalizedSeparation);
				}

				// PS3 [Select the most extreme pair] Choose the pair with the
				// greatest
				// normalized separation along any dimension.
				// Note that if negative it means the rectangles overlapped.
				// However still include
				// overlapping rectangles if that is the only choice available.
				if (normalizedSeparation >= maxNormalizedSeparation) {
					highestLowIndex = tempHighestLowIndex;
					lowestHighIndex = tempLowestHighIndex;
					maxNormalizedSeparation = normalizedSeparation;
				}
			}
		}

		// At this point it is possible that the new rectangle is both
		// highestLow and lowestHigh.
		// This can happen if all rectangles in the node overlap the new
		// rectangle.
		// Resolve this by declaring that the highestLowIndex is the lowest Y
		// and,
		// the lowestHighIndex is the largest X (but always a different
		// rectangle)
		if (highestLowIndex == lowestHighIndex) {
			highestLowIndex = -1;
			float[] tempMinY = newRectMin.clone();
			lowestHighIndex = 0;
			float[] tempMaxX = n.entriesMax[0].clone();

			for (int i = 1; i < n.entryCount; i++) {
				for (int e = 1; e < dim; e++)
				{
				if (n.entriesMin[i][e] < tempMinY[e]) {
					tempMinY[e] = n.entriesMin[i][e];
					highestLowIndex = i;
				} else if (n.entriesMax[i][e-1] > tempMaxX[e-1]) {
					tempMaxX[e-1] = n.entriesMax[i][e-1];
					lowestHighIndex = i;
				}
				}
			}
		}

		// highestLowIndex is the seed for the new node.
		if (highestLowIndex == -1) {
			newNode.addEntry(newRectMin.clone(), newRectMax.clone(), newId);
		} else {
			newNode.addEntry(n.entriesMin[highestLowIndex].clone(), n.entriesMax[highestLowIndex].clone(), n.ids[highestLowIndex]);
			n.ids[highestLowIndex] = -1;

			// move the new rectangle into the space vacated by the seed for the
			// new node
			n.entriesMin[highestLowIndex] = newRectMin.clone();
			n.entriesMax[highestLowIndex] = newRectMax.clone();

			n.ids[highestLowIndex] = newId;
		}

		// lowestHighIndex is the seed for the original node.
		if (lowestHighIndex == -1) {
			lowestHighIndex = highestLowIndex;
		}

		entryStatus[lowestHighIndex] = ENTRY_STATUS_ASSIGNED;
		n.entryCount = 1;
		for (int i = 0; i < dim; i++)
		{
			n.mbb.minCoords[i] = n.entriesMin[lowestHighIndex][i];
			n.mbb.maxCoords[i] = n.entriesMax[lowestHighIndex][i];
		}
	}

	/**
	 * Pick the next entry to be assigned to a group during a node split.
	 *
	 * [Determine cost of putting each entry in each group] For each entry not
	 * yet in a group, calculate the area increase required in the covering
	 * rectangles of each group
	 */
	private int pickNext(Node n, Node newNode) {
		float maxDifference = Float.NEGATIVE_INFINITY;
		int next = 0;
		int nextGroup = 0;

		maxDifference = Float.NEGATIVE_INFINITY;

		if (log.isDebugEnabled()) {
			log.debug("pickNext()");
		}

		for (int i = 0; i < maxNodeEntries; i++) {
			if (entryStatus[i] == ENTRY_STATUS_UNASSIGNED) {

				if (n.ids[i] == -1) {
					log.error("Error: Node " + n.nodeId + ", entry " + i + " is null");
				}

				Rectangle entrieRecty = new Rectangle(n.entriesMin[i].clone(), n.entriesMax[i].clone());
				float nIncrease = n.mbb.enlargement(entrieRecty);
				float newNodeIncrease = newNode.mbb.enlargement(entrieRecty);

				float difference = Math.abs(nIncrease - newNodeIncrease);

				if (difference > maxDifference) {
					next = i;

					if (nIncrease < newNodeIncrease) {
						nextGroup = 0;
					} else if (newNodeIncrease < nIncrease) {
						nextGroup = 1;
					} else if (n.mbb.area() < newNode.mbb.area()) {
						nextGroup = 0;
					} else if (newNode.mbb.area() < n.mbb.area()) {
						nextGroup = 1;
					} else if (newNode.entryCount < maxNodeEntries / 2) {
						nextGroup = 0;
					} else {
						nextGroup = 1;
					}
					maxDifference = difference;
				}
				if (log.isDebugEnabled()) {
					log.debug("Entry " + i + " group0 increase = " + nIncrease + ", group1 increase = "
							+ newNodeIncrease + ", diff = " + difference + ", MaxDiff = " + maxDifference + " (entry "
							+ next + ")");
				}
			}
		}

		entryStatus[next] = ENTRY_STATUS_ASSIGNED;

		if (nextGroup == 0) {
			for (int i = 0; i < dim; i++) {
				if (n.entriesMin[next][i] < n.mbb.minCoords[i])
					n.mbb.minCoords[i] = n.entriesMin[next][i];
				if (n.entriesMax[next][i] > n.mbb.maxCoords[i])
					n.mbb.maxCoords[i] = n.entriesMax[next][i];
			}
			n.entryCount++;
		} else {
			// move to new node.
			newNode.addEntry(n.entriesMin[next].clone(), n.entriesMax[next].clone(), n.ids[next]);
			n.ids[next] = -1;
		}

		return next;
	}

	/**
	 * Recursively searches the tree for the nearest entry. Other queries call
	 * execute() on an IntProcedure when a matching entry is found; however
	 * nearest() must store the entry Ids as it searches the tree, in case a
	 * nearer entry is found. Uses the member variable nearestIds to store the
	 * nearest entry IDs (it is an array, rather than a single value, in case
	 * multiple entries are equally near)
	 */
	private float nearest(Point p, Node n, float furthestDistanceSq, TIntArrayList nearestIds) {
		for (int i = 0; i < n.entryCount; i++) {
			float tempDistanceSq = Rectangle.distanceSq(n.entriesMin[i], n.entriesMax[i], p);
			if (n.isLeaf()) { // for leaves, the distance is an actual nearest
								// distance
				if (tempDistanceSq < furthestDistanceSq) {
					furthestDistanceSq = tempDistanceSq;
					nearestIds.reset();
				}
				if (tempDistanceSq <= furthestDistanceSq) {
					nearestIds.add(n.ids[i]);
				}
			} else { // for index nodes, only go into them if they potentially
						// could have
						// a rectangle nearer than actualNearest
				if (tempDistanceSq <= furthestDistanceSq) {
					// search the child node
					furthestDistanceSq = nearest(p, getNode(n.ids[i]), furthestDistanceSq, nearestIds);
				}
			}
		}
		return furthestDistanceSq;
	}

	/**
	 * Recursively searches the tree for all intersecting entries. Immediately
	 * calls execute() on the passed IntProcedure when a matching entry is
	 * found.
	 *
	 * TODO rewrite this to be non-recursive? Make sure it doesn't slow it down.
	 */
	private boolean intersects(Rectangle r, TIntProcedure v, Node n) {
		for (int i = 0; i < n.entryCount; i++) {
			Rectangle rn = new Rectangle(n.entriesMin[i], n.entriesMax[i]);
			if (r.intersects(rn)) {
				if (n.isLeaf()) {
					if (!v.execute(n.ids[i])) {
						return false;
					}
				} else {
					Node childNode = getNode(n.ids[i]);
					if (!intersects(r, v, childNode)) {
						return false;
					}
				}
			}
		}
		return true;
	}

	/**
	 * Used by delete(). Ensures that all nodes from the passed node up to the
	 * root have the minimum number of entries.
	 *
	 * Note that the parent and parentEntry stacks are expected to contain the
	 * nodeIds of all parents up to the root.
	 */
	private void condenseTree(Node l) {
		// CT1 [Initialize] Set n=l. Set the list of eliminated
		// nodes to be empty.
		Node n = l;
		Node parent = null;
		int parentEntry = 0;

		TIntStack eliminatedNodeIds = new TIntArrayStack();

		// CT2 [Find parent entry] If N is the root, go to CT6. Otherwise
		// let P be the parent of N, and let En be N's entry in P
		while (n.level != treeHeight) {
			parent = getNode(parents.pop());
			parentEntry = parentsEntry.pop();

			// CT3 [Eliminiate under-full node] If N has too few entries,
			// delete En from P and add N to the list of eliminated nodes
			if (n.entryCount < minNodeEntries) {
				parent.deleteEntry(parentEntry);
				eliminatedNodeIds.push(n.nodeId);
			} else {
				// CT4 [Adjust covering rectangle] If N has not been eliminated,
				// adjust EnI to tightly contain all entries in N
				float[] deletedMin = new float[dim];
				float[] deletedMax = new float[dim];
				boolean needsRecalc = false;
				int i = 0;
				while (!needsRecalc && i < dim) {
					if (n.mbb.minCoords[i] != parent.entriesMin[parentEntry][i]
							|| n.mbb.maxCoords[i] != parent.entriesMax[parentEntry][i]) {
						needsRecalc = true;
					}
					i++;
				}
				if (needsRecalc) {
					for (i = 0; i < dim; i++) {
						deletedMin[i] = parent.entriesMin[parentEntry][i];
						deletedMax[i] = parent.entriesMax[parentEntry][i];
						parent.entriesMin[parentEntry][i] = n.mbb.minCoords[i];
						parent.entriesMax[parentEntry][i] = n.mbb.maxCoords[i];
					}
					parent.recalculateMBRIfInfluencedBy(deletedMin, deletedMax);
				}
			}
			// CT5 [Move up one level in tree] Set N=P and repeat from CT2
			n = parent;
		}

		// CT6 [Reinsert orphaned entries] Reinsert all entries of nodes in set
		// Q.
		// Entries from eliminated leaf nodes are reinserted in tree leaves as
		// in
		// Insert(), but entries from higher level nodes must be placed higher
		// in
		// the tree, so that leaves of their dependent subtrees will be on the
		// same
		// level as leaves of the main tree
		while (eliminatedNodeIds.size() > 0) {
			Node e = getNode(eliminatedNodeIds.pop());
			for (int j = 0; j < e.entryCount; j++) {
				add(e.entriesMin[j].clone(), e.entriesMax[j].clone(), e.ids[j], e.level);
				e.ids[j] = -1;
			}
			e.entryCount = 0;
			deletedNodeIds.push(e.nodeId);
		}
	}

	/**
	 * Used by add(). Chooses a leaf to add the rectangle to.
	 */
	private Node chooseNode(float[] minCoords, float[] maxCoords, int level) {
		// CL1 [Initialize] Set N to be the root node
		Node n = getNode(rootNodeId);
		parents.clear();
		parentsEntry.clear();

		Rectangle coordsrect = new Rectangle(minCoords, maxCoords);

		// CL2 [Leaf check] If N is a leaf, return N
		while (true) {
			if (n == null) {
				log.error("Could not get root node (" + rootNodeId + ")");
			}

			if (n.level == level) {
				return n;
			}

			// CL3 [Choose subtree] If N is not at the desired level, let F be
			// the entry in N
			// whose rectangle FI needs least enlargement to include EI. Resolve
			// ties by choosing the entry with the rectangle of smaller area.
			Rectangle leastRect = new Rectangle(n.entriesMin[0], n.entriesMax[0]);
			float leastEnlargement = leastRect.enlargement(coordsrect);
			int index = 0; // index of rectangle in subtree
			for (int i = 1; i < n.entryCount; i++) {
				float[] tempMinX = n.entriesMin[i].clone();
				float[] tempMaxX = n.entriesMax[i].clone();
				Rectangle tmpRect = new Rectangle(tempMinX, tempMaxX);
				float tempEnlargement = tmpRect.enlargement(coordsrect);
				Rectangle indexRect = new Rectangle(n.entriesMin[index], n.entriesMax[index]);
				if ((tempEnlargement < leastEnlargement) || ((tempEnlargement == leastEnlargement)
						&& (Rectangle.area(tmpRect) < Rectangle.area(indexRect)))) {
					index = i;
					leastEnlargement = tempEnlargement;
				}
			}

			parents.push(n.nodeId);
			parentsEntry.push(index);

			// CL4 [Descend until a leaf is reached] Set N to be the child node
			// pointed to by Fp and repeat from CL2
			n = getNode(n.ids[index]);
		}
	}

	/**
	 * Ascend from a leaf node L to the root, adjusting covering rectangles and
	 * propagating node splits as necessary.
	 */
	private Node adjustTree(Node n, Node nn) {
		// AT1 [Initialize] Set N=L. If L was split previously, set NN to be
		// the resulting second node.

		// AT2 [Check if done] If N is the root, stop
		while (n.level != treeHeight) {

			// AT3 [Adjust covering rectangle in parent entry] Let P be the
			// parent
			// node of N, and let En be N's entry in P. Adjust EnI so that it
			// tightly
			// encloses all entry rectangles in N.
			Node parent = getNode(parents.pop());
			int entry = parentsEntry.pop();

			if (parent.ids[entry] != n.nodeId) {
				log.error("Error: entry " + entry + " in node " + parent.nodeId + " should point to node " + n.nodeId
						+ "; actually points to node " + parent.ids[entry]);
			}

			boolean needsrecalc = false;
			int i = 0;
			while (!needsrecalc && i < dim) {
				if (parent.entriesMin[entry][i] != n.mbb.minCoords[i]
						|| parent.entriesMax[entry][i] != n.mbb.maxCoords[i]) {
					needsrecalc = true;
				}
				i++;
			}
			if (needsrecalc) {
				for (i = 0; i < dim; i++) {
					parent.entriesMin[entry][i] = n.mbb.minCoords[i];
					parent.entriesMax[entry][i] = n.mbb.maxCoords[i];
				}

				parent.recalculateMBR();
			}

			// AT4 [Propagate node split upward] If N has a partner NN resulting
			// from
			// an earlier split, create a new entry Enn with Ennp pointing to NN
			// and
			// Enni enclosing all rectangles in NN. Add Enn to P if there is
			// room.
			// Otherwise, invoke splitNode to produce P and PP containing Enn
			// and
			// all P's old entries.
			Node newNode = null;
			if (nn != null) {
				if (parent.entryCount < maxNodeEntries) {
					parent.addEntry(nn.mbb.minCoords.clone(), nn.mbb.maxCoords.clone(), nn.nodeId);
				} else {
					newNode = splitNode(parent, nn.mbb.minCoords.clone(), nn.mbb.maxCoords.clone(), nn.nodeId);
				}
			}

			// AT5 [Move up to next level] Set N = P and set NN = PP if a split
			// occurred. Repeat from AT2
			n = parent;
			nn = newNode;

			parent = null;
			newNode = null;
		}

		return nn;
	}

	/**
	 * Check the consistency of the tree.
	 *
	 * @return false if an inconsistency is detected, true otherwise.
	 */
	public boolean checkConsistency() {
		return checkConsistency(rootNodeId, treeHeight, null);
	}

	private boolean checkConsistency(int nodeId, int expectedLevel, Rectangle expectedMBR) {
		// go through the tree, and check that the internal data structures of
		// the tree are not corrupted.
		Node n = getNode(nodeId);

		if (n == null) {
			log.error("Error: Could not read node " + nodeId);
			return false;
		}

		// if tree is empty, then there should be exactly one node, at level 1
		// TODO: also check the MBR is as for a new node
		if (nodeId == rootNodeId && size() == 0) {
			if (n.level != 1) {
				log.error("Error: tree is empty but root node is not at level 1");
				return false;
			}
		}

		if (n.level != expectedLevel) {
			log.error("Error: Node " + nodeId + ", expected level " + expectedLevel + ", actual level " + n.level);
			return false;
		}

		Rectangle calculatedMBR = calculateMBR(n);
		Rectangle actualMBR = new Rectangle(dim);
		actualMBR.minCoords = n.mbb.minCoords.clone();
		actualMBR.maxCoords = n.mbb.maxCoords.clone();
		if (!actualMBR.equals(calculatedMBR) && log.isErrorEnabled()) {
			log.error("Error: Node " + nodeId + ", calculated MBR does not equal stored MBR");
			for (int i = 0; i < dim; i++) {
				//if (actualMBR.minCoords[i] != n.mbb.minCoords[i])
					log.error("  actualMin[" + i + "]=" + actualMBR.minCoords[i] + ", calc="
							+ calculatedMBR.minCoords[i]);
					log.error("  actualMax[" + i + "]=" + actualMBR.maxCoords[i] + ", calc="
							+ calculatedMBR.maxCoords[i]);
			}
			return false;
		}

		if (expectedMBR != null && !actualMBR.equals(expectedMBR)) {
			log.error("Error: Node " + nodeId + ", expected MBR (from parent) does not equal stored MBR");
			return false;
		}

		// Check for corruption where a parent entry is the same object as the
		// child MBR
		if (expectedMBR != null && actualMBR.sameObject(expectedMBR)) {
			log.error("Error: Node " + nodeId + " MBR using same rectangle object as parent's entry");
			return false;
		}

		for (int i = 0; i < n.entryCount; i++) {
			if (n.ids[i] == -1) {
				log.error("Error: Node " + nodeId + ", Entry " + i + " is null");
				return false;
			}

			if (n.level > 1) { // if not a leaf
				if (!checkConsistency(n.ids[i], n.level - 1, new Rectangle(n.entriesMin[i], n.entriesMax[i]))) {
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * Given a node object, calculate the node MBR from it's entries. Used in
	 * consistency checking
	 */
	private Rectangle calculateMBR(Node n) {
		Rectangle mbr = new Rectangle(dim);

		for (int i = 0; i < n.entryCount; i++) {
			for (int e = 0; e < dim; e++) {
				if (n.entriesMin[i][e] < mbr.minCoords[e])
					mbr.minCoords[e] = n.entriesMin[i][e];
				if (n.entriesMax[i][e] > mbr.maxCoords[e])
					mbr.maxCoords[e] = n.entriesMax[i][e];
			}
		}
		return mbr;
	}
}
