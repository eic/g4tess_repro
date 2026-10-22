// Standalone reproducer: G4TessellatedSolid::Inside() returns kInside for a
// point that is at least 2.7 mm outside a closed, correctly oriented solid
// (a tetrahedron plus a small disjoint cube).
//
// One face (A,B,C) is a long sliver triangle (two ~129 mm edges, one 0.09 mm
// edge), taken verbatim from a CAD-exported mesh. The first fixed ray
// direction used by Inside() (fRandir[0]) misses this triangle by ~50 nm,
// past its short edge, but G4TriangularFacet::Intersect() reports an
// "outgoing" hit anyway. With no matching incoming hit, Inside() concludes
// the point is inside.
//
// Build:  g++ -std=c++17 repro.cc $(geant4-config --cflags --libs) -o repro
//    or:  cmake -B build && cmake --build build
// Exit code 1 means the bug was reproduced.
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include <cstdio>

static const char* name(EInside e)
{ return e == kInside ? "kInside" : e == kSurface ? "kSurface" : "kOutside"; }

int main()
{
  // Sliver triangle, verbatim from the CAD mesh (mm)
  const G4ThreeVector A(-18.78619995117188,  -2.1279968261718523, -256.50000152587893);
  const G4ThreeVector B(-17.95416870117188,  -2.3019165039062273, -127.31544647216796);
  const G4ThreeVector C(-17.879125976562506, -2.3171752929687273, -127.36312255859374);
  // Apex placed 1 mm behind the sliver (opposite to its outward normal)
  const G4ThreeVector n = (B - A).cross(C - A).unit();
  const G4ThreeVector D = (A + B + C) / 3. - 1.0 * mm * n;

  auto* solid = new G4TessellatedSolid("tet");
  // Outward-facing (anticlockwise seen from outside) faces
  solid->AddFacet(new G4TriangularFacet(A, B, C, ABSOLUTE));
  solid->AddFacet(new G4TriangularFacet(A, D, B, ABSOLUTE));
  solid->AddFacet(new G4TriangularFacet(B, D, C, ABSOLUTE));
  solid->AddFacet(new G4TriangularFacet(C, D, A, ABSOLUTE));

  const G4ThreeVector p(-15.261, -3.06546, -127.6);

  // A second, disjoint closed shell (small cube) whose only purpose is to
  // enlarge the bounding box so that p is inside the solid's extent;
  // otherwise Inside() returns kOutside early via OutsideOfExtent().
  // The cube lies at +x from p, i.e. behind the ray direction used below.
  {
    const G4double x0 = -14.0, x1 = -13.0, y0 = -2.5, y1 = -2.4, z0 = -200.0, z1 = -199.0;
    G4ThreeVector c[8] = {{x0,y0,z0},{x1,y0,z0},{x1,y1,z0},{x0,y1,z0},
                          {x0,y0,z1},{x1,y0,z1},{x1,y1,z1},{x0,y1,z1}};
    auto quad = [&](int a, int b, int cc, int d) {   // anticlockwise from outside
      solid->AddFacet(new G4TriangularFacet(c[a], c[b], c[cc], ABSOLUTE));
      solid->AddFacet(new G4TriangularFacet(c[a], c[cc], c[d], ABSOLUTE));
    };
    quad(0,3,2,1); quad(4,5,6,7); quad(0,1,5,4);
    quad(2,3,7,6); quad(1,2,6,5); quad(0,4,7,3);
  }
  solid->SetSolidClosed(true);

  // Independent ground truth for a convex solid: the point is outside if it
  // lies on the outer side of any face plane.
  // (p is trivially outside the cube: p.x < -14)
  G4bool outside = false;
  const G4ThreeVector faces[4][3] = {{A,B,C},{A,D,B},{B,D,C},{C,D,A}};
  for (auto& f : faces) {
    G4ThreeVector fn = (f[1]-f[0]).cross(f[2]-f[0]).unit();
    G4double s = fn.dot(p - f[0]);
    printf("  signed distance to face plane: %+.6f mm\n", s);
    if (s > 0) outside = true;
  }
  printf("Ground truth           : %s\n", outside ? "OUTSIDE" : "INSIDE");
  printf("Cubic volume (sanity)  : %g mm3 (positive => outward normals)\n",
         solid->GetCubicVolume());
  printf("G4TessellatedSolid     : Inside(p) = %s\n", name(solid->Inside(p)));

  // Show the offending facet intersection along fRandir[0]
  const G4ThreeVector v(-0.9577428892113370, 0.2732676269591740, 0.0897405271949221);
  G4double dist, distFromSurface; G4ThreeVector normal;
  for (G4int i = 0; i < solid->GetNumberOfFacets(); ++i) {
    G4VFacet* f = solid->GetFacet(i);
    for (G4bool outgoing : {true, false})
      if (f->Intersect(p, v, outgoing, dist, distFromSurface, normal))
        printf("  facet %d Intersect(outgoing=%d): HIT at %.6f mm\n", i, outgoing, dist);
  }
  // Exact barycentric coordinates of the ray/plane crossing for facet ABC
  G4double t = n.dot(A - p) / n.dot(v);
  G4ThreeVector x = p + t * v;
  G4double area2 = (B - A).cross(C - A).mag();
  G4double bu = ((C - B).cross(x - B)).dot(n) / area2;   // barycentric weight of A
  G4double bv = ((A - C).cross(x - C)).dot(n) / area2;   // of B
  G4double bw = ((B - A).cross(x - A)).dot(n) / area2;   // of C
  printf("  exact ray/plane(ABC) crossing at t=%.6f mm, barycentrics (%.3e, %.3e, %.3e)\n",
         t, bu, bv, bw);
  printf("  -> crossing lies %.1f nm outside edge BC\n",
         -bu * area2 / (C - B).mag() * 1e6);
  return (outside && solid->Inside(p) == kInside) ? 1 : 0;
}
