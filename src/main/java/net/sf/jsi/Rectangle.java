//   Rectangle.java
//   Java Spatial Index Library
//   Copyright (C) 2002-2005 Infomatiq Limited
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
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

package net.sf.jsi;

import java.util.Arrays;

/**
 * Currently hardcoded to 2 dimensions, but could be extended.
 */
public class Rectangle {

	/**
	 * use primitives instead of arrays for the coordinates of the rectangle, to
	 * reduce memory requirements.
	 */
	public float[] minCoords, maxCoords;

	public Rectangle(int dim) {
		minCoords = new float[dim];
		Arrays.fill(minCoords, Float.POSITIVE_INFINITY);
		maxCoords = new float[dim];
		Arrays.fill(maxCoords, Float.NEGATIVE_INFINITY);
	}

	/**
	 * Constructor.
	 *
	 * @param x1
	 *            coordinate of any corner of the rectangle
	 * @param y1
	 *            (see x1)
	 * @param x2
	 *            coordinate of the opposite corner
	 * @param y2
	 *            (see x2)
	 */
	public Rectangle(float[] minCoords, float[] maxCoords) {
		this.minCoords = minCoords.clone();
		this.maxCoords = maxCoords.clone();
	}

	/**
	 * Make a copy of this rectangle
	 *
	 * @return copy of this rectangle
	 */
	public Rectangle copy() {
		return new Rectangle(minCoords.clone(), maxCoords.clone());
	}

	/**
	 * Determine whether an edge of this rectangle overlies the equivalent edge
	 * of the passed rectangle
	 */
	public boolean edgeOverlaps(Rectangle r) {
		boolean overlaps = false;
		for (int i = 0; i < minCoords.length; i++)
		{
			if (minCoords[i] == r.minCoords[i] || maxCoords[i] == r.maxCoords[i])
			{
				overlaps |= true;
			}
		}
		return overlaps;
	}

	/**
	 * Determine whether this rectangle intersects the passed rectangle
	 *
	 * @param r
	 *            The rectangle that might intersect this rectangle
	 *
	 * @return true if the rectangles intersect, false if they do not intersect
	 */
	public boolean intersects(Rectangle r) {
		boolean intersects = true;
		int i = 0;
		while (i < this.minCoords.length && intersects)
		{
			// überprüfen, ob in der jeweiligen Dimension eine Überschneidiung
			// vorliegt
			// Bedingung ist erfüllt, wenn keine vorliegt
			if (r.minCoords[i] > this.maxCoords[i] || r.maxCoords[i] < this.minCoords[i])
			{
				intersects &= false;
			}
			i++;
		}
		return intersects;
	}

	/**
	 * Determine whether or not two rectangles intersect
	 *
	 * @param r1MinX
	 *            minimum X coordinate of rectangle 1
	 * @param r1MinY
	 *            minimum Y coordinate of rectangle 1
	 * @param r1MaxX
	 *            maximum X coordinate of rectangle 1
	 * @param r1MaxY
	 *            maximum Y coordinate of rectangle 1
	 * @param r2MinX
	 *            minimum X coordinate of rectangle 2
	 * @param r2MinY
	 *            minimum Y coordinate of rectangle 2
	 * @param r2MaxX
	 *            maximum X coordinate of rectangle 2
	 * @param r2MaxY
	 *            maximum Y coordinate of rectangle 2
	 *
	 * @return true if r1 intersects r2, false otherwise.
	 */
	static public boolean intersects(Rectangle r1, Rectangle r2) {
		return r1.intersects(r2);
	}

	/**
	 * Determine whether this rectangle contains the passed rectangle
	 *
	 * @param r
	 *            The rectangle that might be contained by this rectangle
	 *
	 * @return true if this rectangle contains the passed rectangle, false if it
	 *         does not
	 */
	public boolean contains(Rectangle r) {
		return this.contains(r.minCoords, r.maxCoords);
	}
	
	/**
	 * Determine whether this rectangle contains the passed rectangle
	 *
	 * @param r
	 *            The rectangle that might be contained by this rectangle
	 *
	 * @return true if this rectangle contains the passed rectangle, false if it
	 *         does not
	 */
	public boolean contains(float[] minCoords, float[] maxCoords) {
		boolean isContained = true;
		int i = 0;
		while (i < this.minCoords.length && isContained)
		{
			// Achtung: inverse Bedingung gegenüber obigem Kommentar!
			// es wird auf nicht Enthaltensein geprüft!
			if (minCoords[i] < this.minCoords[i] || maxCoords[i] > this.maxCoords[i])
			{
				isContained &= false;
			}
			i++;
		}
		return isContained;
	}

	/**
	 * Determine whether this rectangle is contained by the passed rectangle
	 *
	 * @param r
	 *            The rectangle that might contain this rectangle
	 *
	 * @return true if the passed rectangle contains this rectangle, false if it
	 *         does not
	 */
	public boolean containedBy(Rectangle r) {
		boolean contained = true;
		int i = 0;
		while (i < minCoords.length)
		{
			if(r.maxCoords[i] < maxCoords[i] || r.minCoords[i] > minCoords[i])
			{
				return false;
			}
		}
		return contained;
	}

	/**
	 * Return the distance between this rectangle and the passed point. If the
	 * rectangle contains the point, the distance is zero.
	 *
	 * @param p
	 *            Point to find the distance to
	 *
	 * @return distance beween this rectangle and the passed point.
	 */
	public float distance(Point p) {

		return (float) Math.sqrt(Rectangle.distanceSq(minCoords, maxCoords, p));
	}

	/**
	 * Return the distance between a rectangle and a point. If the rectangle
	 * contains the point, the distance is zero.
	 *
	 * @param minX
	 *            minimum X coordinate of rectangle
	 * @param minY
	 *            minimum Y coordinate of rectangle
	 * @param maxX
	 *            maximum X coordinate of rectangle
	 * @param maxY
	 *            maximum Y coordinate of rectangle
	 * @param pX
	 *            X coordinate of point
	 * @param pY
	 *            Y coordinate of point
	 *
	 * @return distance beween this rectangle and the passed point.
	 */
	static public float distance(float[] minCoords, float[] maxCoords, Point p) {
		return (float) Math.sqrt(distanceSq(minCoords, maxCoords, p));
	}

	static public float distanceSq(float[] minCoords, float[] maxCoords, Point p) {

		float distanceSq = 0;
		for (int i = 0; i < minCoords.length; i++) {
			float distanceSqX = 0.0f;
			if (minCoords[i] > p.coords[i]) {
				distanceSqX = minCoords[i] - p.coords[i];
				distanceSqX *= distanceSqX;
			} else if (p.coords[i] > maxCoords[i]) {
				distanceSqX = p.coords[i] - maxCoords[i];
				distanceSqX *= distanceSqX;
			}
			distanceSq += distanceSqX;
		}

		return distanceSq;
	}

	/**
	 * Return the distance between this rectangle and the passed rectangle. If
	 * the rectangles overlap, the distance is zero.
	 *
	 * @param r
	 *            Rectangle to find the distance to
	 *
	 * @return distance between this rectangle and the passed rectangle
	 */

	public float distance(Rectangle r) {
		float distanceSq = 0;
		
		for (int i = 0; i < minCoords.length; i++) {
			float greatestMin = Math.max(minCoords[i], r.minCoords[i]);
			float leastMax = Math.min(maxCoords[i], r.maxCoords[i]);
			if (greatestMin > leastMax) {
				distanceSq += ((greatestMin - leastMax) * (greatestMin - leastMax));
			}
		}
		return (float) Math.sqrt(distanceSq);
	}

	/**
	 * Calculate the area by which this rectangle would be enlarged if added to
	 * the passed rectangle. Neither rectangle is altered.
	 *
	 * @param r
	 *            Rectangle to union with this rectangle, in order to compute
	 *            the difference in area of the union and the original rectangle
	 *
	 * @return enlargement
	 */
	public float enlargement(Rectangle r) {
		float enlargedArea = 1;
		for (int i = 0; i < minCoords.length; i++) {
			enlargedArea *= (Math.max(maxCoords[i], r.maxCoords[i]) - Math.min(minCoords[i], r.minCoords[i]));
		}
		return enlargedArea - area();
	}

	/**
	 * Compute the area of this rectangle.
	 *
	 * @return The area of this rectangle
	 */
	public float area() {
		return area(this);
	}

	/**
	 * Compute the area of a rectangle.
	 *
	 * @param minX
	 *            the minimum X coordinate of the rectangle
	 * @param minY
	 *            the minimum Y coordinate of the rectangle
	 * @param maxX
	 *            the maximum X coordinate of the rectangle
	 * @param maxY
	 *            the maximum Y coordinate of the rectangle
	 *
	 * @return The area of the rectangle
	 */
	static public float area(Rectangle r) {
		float vol = r.maxCoords[0] - r.minCoords[0];
		for (int i = 1; i < r.minCoords.length; i++)
		{
			vol *= (r.maxCoords[i] - r.minCoords[i]);
		}
		return vol;
	}

	/**
	 * Computes the union of this rectangle and the passed rectangle, storing
	 * the result in this rectangle.
	 *
	 * @param r
	 *            Rectangle to add to this rectangle
	 */
	public void add(Rectangle r) {
		for (int i = 0; i < minCoords.length; i++)
		{
			if (r.minCoords[i] < minCoords[i])
			{
				minCoords[i] = r.minCoords[i];
			}
			if (r.maxCoords[i] > maxCoords[i])
			{
				maxCoords[i] = r.maxCoords[i];
			}
		}
	}

	/**
	 * Computes the union of this rectangle and the passed point, storing the
	 * result in this rectangle.
	 *
	 * @param p
	 *            Point to add to this rectangle
	 */
	public void add(Point p) {

		for (int i = 0; i < minCoords.length; i++)
		{
			if (p.coords[i] < minCoords[i])
			{
				minCoords[i] = p.coords[i];
			}
			if (p.coords[i] > maxCoords[i])
			{
				maxCoords[i] = p.coords[i];
			}
		}
	}

	/**
	 * Find the the union of this rectangle and the passed rectangle. Neither
	 * rectangle is altered
	 *
	 * @param r
	 *            The rectangle to union with this rectangle
	 */
	public Rectangle union(Rectangle r) {
		Rectangle union = this.copy();
		union.add(r);
		return union;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		for (int i = 0; i < minCoords.length; i++)
		{
			result = prime * result + Float.floatToIntBits(this.maxCoords[i]);
			result = prime * result + Float.floatToIntBits(this.minCoords[i]);
		}
		return result;
	}

	/**
	 * Determine whether this rectangle is equal to a given object. Equality is
	 * determined by the bounds of the rectangle.
	 *
	 * @param o
	 *            The object to compare with this rectangle
	 */
	@Override
	public boolean equals(Object o) {
		boolean equals = true;
		if (o instanceof Rectangle) {
			Rectangle r = (Rectangle) o;
			for (int i = 0; i < minCoords.length; i++)
			{
				if (minCoords[i] != r.minCoords[i] || maxCoords[i] != r.maxCoords[i])
				{
					return false;
				}
			}
		}
		return equals;
	}

	/**
	 * Determine whether this rectangle is the same as another object
	 *
	 * Note that two rectangles can be equal but not the same object, if they
	 * both have the same bounds.
	 *
	 * @param o
	 *            The object to compare with this rectangle.
	 */
	public boolean sameObject(Object o) {
		return super.equals(o);
	}

	public Point centre() {
		float[] centre = new float[minCoords.length];
		for (int i = 0; i < minCoords.length; i++)
		{
			centre[i] = (minCoords[i] + maxCoords[i]) / 2.0f;
		}
		return new Point(centre);
	}

}
