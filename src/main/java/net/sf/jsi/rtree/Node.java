//   Node.java
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

import java.io.Serializable;
import java.util.Arrays;

import net.sf.jsi.Rectangle;

/**
 * <p>
 * Used by RTree. There are no public methods in this class.
 * </p>
 */
public class Node implements Serializable {
	private static final long serialVersionUID = -2823316966528817396L;
	int nodeId = 0;
	final Rectangle mbb;

	float[][] entriesMin = null;
	float[][] entriesMax = null;

	int[] ids = null;
	int level;
	int entryCount;
	int dim;

	Node(int nodeId, int level, int maxNodeEntries, int dim) {
		this.nodeId = nodeId;
		this.level = level;
		this.dim = dim;
		mbb = new Rectangle(dim);
		entriesMin = new float[maxNodeEntries][];
		entriesMax = new float[maxNodeEntries][];
		for (int i = 0; i < maxNodeEntries; i++) {
			entriesMin[i] = new float[dim];
			Arrays.fill(entriesMin[i], Float.POSITIVE_INFINITY);
			entriesMax[i] = new float[dim];
			Arrays.fill(entriesMax[i], Float.NEGATIVE_INFINITY);
		}
		ids = new int[maxNodeEntries];
	}

	void addEntry(float[] minCoords, float[] maxCoords, int id) {
		ids[entryCount] = id;
		for (int i = 0; i < dim; i++) {
			entriesMin[entryCount][i] = minCoords[i];
			entriesMax[entryCount][i] = maxCoords[i];

			if (minCoords[i] < mbb.minCoords[i]) {
				mbb.minCoords[i] = minCoords[i];
			}
			if (maxCoords[i] > mbb.maxCoords[i]) {
				mbb.maxCoords[i] = maxCoords[i];
			}
		}
		entryCount++;
	}

	// Return the index of the found entry, or -1 if not found
	int findEntry(float[] minCoords, float[] maxCoords, int id) {

		for (int i = 0; i < entryCount; i++) {
			if (id == ids[i]) {
				int e = 0;
				boolean match = true;
				while (match && e < dim) {
					if (entriesMin[i][e] != minCoords[e] || entriesMax[i][e] != maxCoords[e]) {
						match = false;
					}
					e++;
				}
				if (match) {
					return i;
				}
			}
		}
		return -1;
	}

	// delete entry. This is done by setting it to null and copying the last
	// entry into its space.
	void deleteEntry(int i) {
		int lastIndex = entryCount - 1;
		float[] deletedMin = new float[dim];
		float[] deletedMax = new float[dim];
		for (int e = 0; e < dim; e++) {
			deletedMax[e] = entriesMax[i][e];
			deletedMin[e] = entriesMin[i][e];
		}

		if (i != lastIndex) {
			for (int e = 0; e < dim; e++) {
				entriesMax[i][e] = entriesMax[lastIndex][e];
				entriesMin[i][e] = entriesMin[lastIndex][e];
			}
			ids[i] = ids[lastIndex];
		}
		entryCount--;

		// adjust the MBR
		recalculateMBRIfInfluencedBy(deletedMin, deletedMax);
	}

	// deletedMin/MaxX/Y is a rectangle that has just been deleted or made
	// smaller.
	// Thus, the MBR is only recalculated if the deleted rectangle influenced
	// the old MBR
	void recalculateMBRIfInfluencedBy(float[] deletedMin, float[] deletedMax) {
		boolean needsRecalc = false;
		int e = 0;
		while (!needsRecalc && e < dim) {
			if (mbb.minCoords[e] == deletedMin[e] || mbb.maxCoords[e] == deletedMax[e]) {
				needsRecalc = true;
			}
			e++;
		}
		if (needsRecalc) {
			recalculateMBR();
		}
	}

	void recalculateMBR() {
		Arrays.fill(mbb.minCoords, Float.POSITIVE_INFINITY);
		Arrays.fill(mbb.maxCoords, Float.NEGATIVE_INFINITY);

		for (int e = 0; e < entryCount; e++) {
			for (int i = 0; i < dim; i++) {
				if (mbb.minCoords[i] > entriesMin[e][i])
					mbb.minCoords[i] = entriesMin[e][i];
				if (mbb.maxCoords[i] < entriesMax[e][i])
					mbb.maxCoords[i] = entriesMax[e][i];

			}
		}
	}

	/**
	 * eliminate null entries, move all entries to the start of the source node
	 */
	void reorganize(RTree rtree) {
		int countdownIndex = rtree.maxNodeEntries - 1;
		for (int index = 0; index < entryCount; index++) {
			if (ids[index] == -1) {
				while (ids[countdownIndex] == -1 && countdownIndex > index) {
					countdownIndex--;
				}
				for (int i = 0; i < dim; i++) {
					entriesMin[index][i] = entriesMin[countdownIndex][i];
					entriesMax[index][i] = entriesMax[countdownIndex][i];
				}
				ids[index] = ids[countdownIndex];
				ids[countdownIndex] = -1;
			}
		}
	}

	public int getEntryCount() {
		return entryCount;
	}

	public int getId(int index) {
		if (index < entryCount) {
			return ids[index];
		}
		return -1;
	}

	boolean isLeaf() {
		return (level == 1);
	}

	public int getLevel() {
		return level;
	}
}
