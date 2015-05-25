# Running Reuleaux #

From the command line, navigate to the folder containing the binary Reuleaux file. Execute this file followed by the full filename with extension of the PSC file to be meshed. This file must end with the extension .psc. If a configuration file is being used, it must be in the same directory as the PSC file, and must end with the extension .cfg. The results will be stored in files located in the same directory as the PSC file.

The following is a key for the file extensions of the output created by running Reuleaux on a PSC. SVG files can be viewed with a standard web browser. Node and .ele files can be viewed using Triangle. The .log file can be viewed with a standard text editor.

  * .log The logfile containing information about the mesh from several steps of execution.
  * .inp.svg The SVG file of the input provided. Often will not display all features because of later initialization of the bounding box.
  * .pre.svg The SVG file of the input with presplit points added but before initialization of the bounding box.
  * .fsa.svg The SVG file of the input after detecting small angles and initializing the bounding box.
  * .ini.svg The SVG file of the input after initializing the Delaunay triangulation.
  * .out.svg The SVG file of the final mesh triangulation.
  * .out.node The Triangle file of the final mesh triangulation.
  * .out.ele The Triangle file of the final mesh triangulation.
  * .int.svg The SVG file of the final mesh triangulation after removing the mesh outside of exterior features.
  * .rem.svg The SVG file of the final mesh triangulation after removing the mesh outside of exterior features and inside of small angle protection circles.
  * .zoom.svg The SVG file of the final triangulation zoomed in on a selected portion.

The script reuleauxAll.sh is run by execution from a command line. It will run all data files in the ‘data’ folder that is located in the same directory as the directory containing reuleauxAll.sh and reuleauxDir.sh.

# PSC File Format #

The .psc file format is an ASCII text file that specifies a piecewise smooth complex. Its goal is to fully describe a PSC including the intersection points between features. Note that the parser for a PSC file ignores any whitespace except to use it for separating characters. The text formatting, in particular newlines, is only for the user. Integer values can be entered instead of double values without the need for decimal point and following zero. It is critical to ensure that in the input PSC that all intersection points are explicitly specified. Failure to do so will cause the code to crash.

A .psc file consists of 4 sections:
  1. PSC
  1. Dimension
  1. Points
  1. Curves

## Section 1: PSC section ##

This section identifies the file type.

First line: PSC

## Section 2: Dimension section ##

This section specifies the dimensions of the PSC; the Reuleaux code only accepts two-dimensional PSCs.

First line: DIMENSION 2

## Section 3: Points section ##

This section gives the coordinates of all the points in the PSC. It must include the points of all intersections of the curves defined in the PSC. The program requires at least one point to be entered to generate a mesh.

First line: POINTS _#points_

Remaining lines: _pointid xcoord ycoord_


_#points_ is an integer specify the number of points in the PSC.

_pointid_ should be a unique integer identifier.

_xcoord ycoord_ are the double values for the x and y coordinates of the point.

**Example:** Specifies 8 points, 00 is on the origin, 01 is one unit to the right along the x axis, and 02 is one unit up and two to the right of 00.

POINTS 8

00 0.0 0.0

01 1.0 0.0

02 2.0 1.0

03 2.0 0.0

04 0.0 1.5

05 3.0 1.5

06 0.0 4.0

07 3.0 4.0

## Section 4: Curves section ##

This section defines all of the curves and curve segments in the PSC.

First line: CURVES _#curves_

Remaining lines A: _curveid curvetype curvedata_

Remaining lines B: _curvedatacontinued_

Remaining lines C: _curvedatacontinued_

_#curves_ is an integer specify the number of curves in the PSC.

_curveid_ should be a integer identifier unique for each curve.

_curvetype_ is one of the following types that are enumerated below: linesegment, circle, arc, bezier, bezier3, or nbezier and are described in their corresponding section.

_curvedata_ is dependent on the curve type and is specified in the following allowed curve types. Line C is only used with n-degree Bézier curves.

**Example:** Specifies 6 curves.

CURVES 6

### Line Segments ###

A line segment is characterized by:

Line A: _curveid_ linesegment _#segments_

Line B: _pointlist_

_#segments_ is a positive integer specifying the number of segments in a straight line separated by points.

_pointlist_ is a sequence of integer identifiers from the specified list of points in the points section A. There must be a minimum of two points listed, and the number of points specified is equal to one more than the number of segments.

**Example:** Specifies one straight line 100 with two segments, from point 00 to point 01 to point 03.

100 linesegment 2 00 01 03

### Circles ###

A circle is characterized by:

Line A: _curveid_ circle _xcoord ycoord radius_

Line B: _#points pointlist_

_xcoord ycoord_ are double values for the coordinates of the center point of the circle.

_radius_ is a double value for the radius of the circle.

_#points_ is an integer specifying the number of intersection points with other features on the circle.

_pointlist_ is a sequence of integer identifiers for intersection points listed in the points section.

**Example:** Specifies a circle with id 101, centered on the origin with a radius of one unit. There is one point on the circle, 01.

101 circle 0.0 0.0 1.0 1 01

### Circular Arcs ###

A circular arc is characterized by:

Line A: _curveid_ arc _leftpoint rightpoint xcoord ycoord_

Line B: _#points pointlist_

_leftpoint rightpoint_ are the integer identifiers for left and right endpoints of the arc respectively. The right endpoint is clock-wise along the circle from the left endpoint.

_xcoord ycoord_ are double values for the coordinates of the center point of the circle describing the arc.

_#points_ is an integer specifying the number of intersection points with other features on the circular arc.

_pointlist_ is a sequence of integer identifiers for intersection points listed in the points section.

The radius of the arc does not need to be specified because it can be calculated by the distance between either of the endpoints and the center point of the arc.

**Example:** Specifies an arc from point 01 to 02, centered at (2.0, 0.0) thus giving it a radius of one unit. No points other than the endpoints are added to the curve.

102 arc 01 02 2.0 0.0 0

### Quadratic Bézier Curves ###

A quadratic Bézier curve is characterized by:

Line A: _curveid_ bezier _leftpoint rightpoint xctrl yctrl_

Line B: _#points pointlist_

_leftpoint rightpoint_ are the integer identifiers for left and right endpoints of the quadratic Bézier curve respectively.

_xctrl yctrl_ are the double values for the coordinates of the control point of the quadratic Bézier curve.

_#points_ is an integer specifying the number of intersection points with other features on the quadratic Bézier curve.

_pointlist_ is a sequence of integer identifiers for intersection points listed in the points section.

**Example:** Specifies a quadratic Bézier curve from point 02 to 03 that has a control point at (5.0, 0.5). No points other than the endpoints were added to the curve.

103 bezier 02 03 5.0 0.5 0

### Cubic Bézier Curves ###

A cubic Bézier curve is characterized by:

Line A: _curveid_ bezier3 _leftpoint rightpoint xyctrl1 xyctrl2_

Line B: _#points pointlist_

_leftpoint rightpoint_ are the integer identifiers for left and right endpoints of the cubic Bézier curve respectively.

_xyctrl1_ contains the two double values representing the coordinates of the location of the first control point of the cubic Bézier curve.

_xyctrl2_ has the coordinates of the second control point.

_#points_ is an integer specifying the number of intersection points with other features on the cubic Bézier curve.

_pointlist_ is a sequence of integer identifiers for intersection points listed in the points section.

**Example:** Specifies a cubic Bézier curve from 04 to 05 with first control point (1.0,4.0), and second control point (2.00.5). No points other than the endpoints were added to the curve.

104 bezier3 04 05 1.0 4.0 2.0 0.5

0

### n-Degree Bézier Curves ###

A n-degree Bézier curve is characterized by:

Line A: _curveid_ nbezier _leftpoint rightpoint degree #points_

Line B: _x, ycoordlist_

Line C: _pointlist_

_leftpoint rightpoint_ are the integer identifiers for left and right endpoints of the n-degree Bézier curve respectively.

_degree_ is the degree of the curve and specifies the number of control points for the Bézier curve. The number of control points is equal to one less than the degree of the Bézier curve.

_#points_ is an integer specifying the number of intersection points with other features on the n-degree Bézier curve.

_x, ycoordlist_ is a list of pairs of doubles that specify the locations of the control points for the n-degree Bézier curve.

_pointlist_ is a sequence of integer identifiers for intersection points listed in the points section.

**Example:** Specifies a Bézier curve of degree 5, from point 06 to 07, and with control points (0.5, 8.0), (1.0, 2.0), (1.5, 2.0), and (2.7, 8.0). No points other than the endpoints were added to the curve.

105 nbezier 06 07 5 0

0.5 8.0 1.0 2.0 1.5 2.0 2.7 8.0

# Configuration File Format #

The .cfg file format is used to set different configuration options for the meshing of a PSC by Reuleaux (the name for our code). It must be located in the same directory and have the same name of its associated .psc file except for the file extension.

## Configuration Options ##

Each specifies a type and the default option.

```
Epsilon = (double, 1.0e-12)
```

A threshold for numerical zero.

```
Infinity = (double, 1.0e12)
```

A threshold for numerical infinity.

```
MaxRadiusEdgeRatio = (double, 1.5)
```

The maximum allowable circumradius to shortest edge ratio for any triangle in the triangulation.

```
MaxOutputRadius = (double, -1.0)
```

A user specified size parameter for the maximum allowable output circumradius. The final mesh will contain no triangle with circumradius larger than the specified value. The default value of -1.0 means no restriction.

```
BoundingBoxRatio = (double, 6.0)
```

The size of the bounding box used when generating the mesh. Bounding box will be that many times larger than the farthest distance between input features in the x or y directions respectively. To generate a mesh within the bounding box, use at least 6.

```
PresplitVariation = (double, -1.0)
```

The maximum allowed total variation in orientation in radians when presplitting curves. Default value of -1.0 splits circles and arcs with total variation greater than π/2 radians and Bézier curves with total variation greater than π/4 radians.

```
Maxiterations = (integer, -1)
```

Maximum number of iterations allowed during Ruppert’s Algorithm. The default value of -1 means there is no limitation.

```
ExactInsphere = (boolean, true)
```

Use exact arithmetic instead of standard floating point arithmetic. The code is more stable with the default of exact arithmetic.

```
RemoveBoundingBox = (boolean, false)
```

If this option is set to true, the algorithm will remove all triangles associated with the bounding box after processing. If the input is not closed curves this will remove all of the mesh.

```
printNBezierWSeg = (boolean, false)
```

If this option is true, n-degree Bézier curves are displayed in the .svg files using 100, small straight line segments to approximate the curve. If false, the n-degree Bézier curves are output using quadratic Béziers curves between points on the n-degree Bézier curve.

# Full Example #

## From file example.psc ##

PSC

DIMENSION 2

POINTS 8

00 0.0 0.0

01 1.0 0.0

02 2.0 1.0

03 2.0 0.0

04 0.0 1.5

05 3.0 1.5

06 0.0 4.0

07 3.0 4.0

CURVES 6

100 linesegment 1

00 01

101 circle 0.0 0.0 1.0

1 01

102 arc 01 02 2.0 0.0

0

103 bezier 02 03 5.0 0.5

0

104 bezier3 04 05 1.0 4.0 2.0 0.5

0

105 nbezier 06 07 5 0

0.5 8.0 1.0 2.0 1.5 2.0 2.7 8.0

## From file example.cfg ##

```
[Reuleaux]

Epsilon = .000000001

Infinity = 1000000000.0

BoundingBoxRatio = 6.0 ; size of the bounding box relative to the data

MaxRadiusEdgeRatio = 1.2 ; output quality

MaxOutputRadius = -1

Maxiterations = -1

ExactInsphere = false

RemoveBoundingBox = false

PresplitVariation = -1.0

printNBezierWSeg = false
```