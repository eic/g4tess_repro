# G4TessellatedSolid::Inside() returns kInside for a point outside the solid

Standalone reproducer, no input files needed. Confirmed on Geant4 11.4.1 and 11.4.2.

## Build and run

```sh
g++ -std=c++17 repro.cc $(geant4-config --cflags --libs) -o repro && ./repro
# or
cmake -B build && cmake --build build && ./build/repro
```

Exit code 1 means the bug was reproduced.

## Expected vs. actual

The solid is a closed, outward-oriented tetrahedron plus a small disjoint cube
(the cube is only there so the test point is inside the bounding box; without
it `OutsideOfExtent()` returns early). `GetCubicVolume()` is positive.

The test point lies on the outer side of two of the tetrahedron's face planes
(+0.245 mm and +2.72 mm), so it is at least 2.72 mm outside the convex
tetrahedron, and it is trivially outside the cube.

```
Ground truth           : OUTSIDE
G4TessellatedSolid     : Inside(p) = kInside      <-- wrong
  facet 0 Intersect(outgoing=1): HIT at 2.749801 mm
  exact ray/plane(ABC) crossing at t=2.749801 mm, barycentrics (-4.528e-07, 2.062e-01, 7.938e-01)
  -> crossing lies 49.9 nm outside edge BC
```

## Mechanism

* `Inside()` (both `InsideVoxels` and `InsideNoVoxels`) casts a ray along the
  first fixed direction `fRandir[0]` and classifies the point by the nearest
  crossing: if the nearest crossing is outgoing, the point is inside.
* Facet 0 is a sliver triangle taken verbatim from a CAD-exported mesh: two
  edges are ~129 mm long and one is 0.09 mm long (aspect ratio ~10^4).
* The ray hits the facet's plane ~50 nm outside the short edge BC, with
  barycentric weight -4.5e-7. That is far beyond `kCarTolerance` (1e-9 m),
  but `G4TriangularFacet::Intersect(p, v, outgoing=true, ...)` still reports a
  hit.
* No other facet is crossed, so the only crossing is a spurious "outgoing"
  one, and `Inside()` returns `kInside`.

In a real CAD mesh this produces sparse false `kInside` results, clustered
along lines (the set of points whose `fRandir[0]` ray grazes a given sliver
edge). These show up as spurious `GeomVol1002` overlaps in
`/geometry/test/run`. The reported "overlap depth" is just the point's
clearance to the solid, because `DistanceToOut()` is evaluated at a point
that is actually outside.

## Origin

Found in the ePIC silicon outer barrel geometry (`hybrid_OB` branch of
eic/epic). In the Geant4 overlap check, sensor boxes were flagged as
overlapping a tessellated carbon-composite support by 113–341 µm. In fact
every flagged point was outside the support, and its distance to the support
equalled the reported overlap depth.
