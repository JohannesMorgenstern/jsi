//   ReferenceCompareTest.java
//   Java Spatial Index Library
//   Copyright (C) 2002-2005 Infomatiq Limited.
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

import java.util.Properties;

import junit.framework.TestCase;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * ReferenceCompareTest
 *  
 * Generates results used for comparing the performance of the Java Spatial 
 * Index library against alternative implementations.
 */
public class ReferenceCompareTest extends TestCase {

  private static final Logger log = LoggerFactory.getLogger(ReferenceCompareTest.class);
  
  private Script script = new Script();
  private Properties linear_3_6 = props("Linear", 3, 6);
  private Properties linear_5_10 = props("Linear", 5, 10);
  private Properties rstar_1_13 = props("RStar", 1, 13);
  private Properties rstar_6_13 = props("RStar", 6, 13);
  
  protected int entriesToTest = 100;
  
  public ReferenceCompareTest(String s) {
    super(s);
  }

  public int GetNumEntriesToTest() {
	  return 100;
  }

  public Properties props(String treeVariant, int minNodeEntries, int maxNodeEntries) {
    Properties p = new Properties();
    p.setProperty("MinNodeEntries", Integer.toString(minNodeEntries));
    p.setProperty("MaxNodeEntries", Integer.toString(maxNodeEntries)); 
    p.setProperty("TreeVariant", treeVariant);
    p.setProperty("dim", "2");
    return p;
  }
  
  private void runComparisonTest(String scriptName, String referenceType, Properties refProps, String testType, Properties testProps) {
    log.info(scriptName + " - creating reference test results");
    script.run(referenceType, refProps, scriptName, Script.REFERENCE_GENERATE);
    
    log.info(scriptName + " - running comparison test");
    script.run(testType, testProps, scriptName, Script.REFERENCE_COMPARISON);  
  }
  
  
  // 100, 1000, 10,000  
  //    Reference results are generated by the SimpleIndex implementation.
  //    Therefore, compare results for both SIL and JSI implementations.
  //
  // 100,000
  //    Reference result generated by SIL, therefore compare results
  //    for JSI only
  public void testReferenceCompareAllFunctions() {
    log.debug("testReferenceCompareAllFunctions()");
    
    if (entriesToTest >= 100) {
      runComparisonTest("allfunctions-100", "SimpleIndex", null, "SILWrapper", linear_3_6);
      runComparisonTest("allfunctions-100", "SimpleIndex", null, "rtree.RTree", linear_3_6);
    }
    
    if (entriesToTest >= 1000) {      
      runComparisonTest("allfunctions-1000", "SimpleIndex", null, "SILWrapper", linear_3_6);
      runComparisonTest("allfunctions-1000", "SimpleIndex", null, "rtree.RTree", linear_3_6);
    }

    if (entriesToTest >= 10000) {
      runComparisonTest("allfunctions-10000", "SimpleIndex", null, "SILWrapper", linear_3_6);
      runComparisonTest("allfunctions-10000", "SimpleIndex", null, "rtree.RTree", linear_3_6);
    }

    if (entriesToTest >= 100000) {  
      runComparisonTest("allfunctions-100000", "SILWrapper", rstar_1_13, "rtree.RTree", linear_3_6);
    }
  }
  
  public void testReferenceCompareDelete() {
    log.debug("testReferenceCompareDelete()");
    
    if (entriesToTest >= 100) {
      runComparisonTest("delete-100", "SimpleIndex", null, "SILWrapper", linear_3_6);
      runComparisonTest("delete-100", "SimpleIndex", null, "rtree.RTree", linear_3_6);
    }
    
    if (entriesToTest >= 1000) {
      runComparisonTest("delete-1000", "SimpleIndex", null, "SILWrapper", linear_3_6);
      runComparisonTest("delete-1000", "SimpleIndex", null, "rtree.RTree", linear_3_6);
    }
    
    if (entriesToTest >= 10000) {
      runComparisonTest("delete-10000", "SimpleIndex", null, "SILWrapper", linear_3_6);
      runComparisonTest("delete-10000", "SimpleIndex", null, "rtree.RTree", linear_3_6);
    }

    if (entriesToTest >= 100000) {
      runComparisonTest("delete-100000", "SILWrapper", rstar_1_13, "rtree.RTree", linear_3_6);
    }
  }

  public void testReferenceCompareIntersect() {
    log.debug("testReferenceCompareIntersect()");
    
    if (entriesToTest >= 100) {
      runComparisonTest("intersect-100", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("intersect-100", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }
    
    if (entriesToTest >= 1000) {
      runComparisonTest("intersect-1000", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("intersect-1000", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }
    
    if (entriesToTest >= 10000) {
      runComparisonTest("intersect-10000", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("intersect-10000", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }

    if (entriesToTest >= 100000) {
      runComparisonTest("intersect-100000", "SILWrapper", rstar_1_13, "rtree.RTree", linear_5_10);
    }
  }
  
  public void testReferenceCompareNearest() {
    log.debug("testReferenceCompareNearest()");
    
    if (entriesToTest >= 100) {
      runComparisonTest("nearest-100", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("nearest-100", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }
    
    if (entriesToTest >= 1000) {
      runComparisonTest("nearest-1000", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("nearest-1000", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }
    
    if (entriesToTest >= 10000) {
      runComparisonTest("nearest-10000", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("nearest-10000", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }

    if (entriesToTest >= 100000) {
      runComparisonTest("nearest-100000", "SILWrapper", rstar_1_13, "rtree.RTree", linear_5_10);
    }    
  }
  
  public void testReferenceCompareNearestN() {
    log.debug("testReferenceCompareNearestN()");
    
    if (entriesToTest >= 100) {
      runComparisonTest("nearestN-100", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("nearestN-100", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }
    
    if (entriesToTest >= 1000) {
      runComparisonTest("nearestN-1000", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("nearestN-1000", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }
    
    if (entriesToTest >= 10000) {
      runComparisonTest("nearestN-10000", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("nearestN-10000", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }

    if (entriesToTest >= 100000) {
      runComparisonTest("nearestN-100000", "SILWrapper", rstar_1_13, "rtree.RTree", linear_5_10);
    } 
  }

  public void testReferenceCompareContains() {
    log.debug("testReferenceCompareContains()");
    
    if (entriesToTest >= 100) {
      runComparisonTest("contains-100", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("contains-100", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }
    
    if (entriesToTest >= 1000) {
      runComparisonTest("contains-1000", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("contains-1000", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }
    
    if (entriesToTest >= 10000) {
      runComparisonTest("contains-10000", "SimpleIndex", null, "SILWrapper", linear_5_10);
      runComparisonTest("contains-10000", "SimpleIndex", null, "rtree.RTree", linear_5_10);
    }

    if (entriesToTest >= 100000) {
      runComparisonTest("contains-100000", "SILWrapper", rstar_6_13, "rtree.RTree", linear_5_10);
    } 
  }
}
