// CurveGenerator.h - turns the mathematics into points that can be drawn.
//
//  * family curves: explicit  y = f(x, c)  -> sample x and compute y
//                   implicit  F(x, y, c) = 0 -> MARCHING SQUARES on a grid
//  * orthogonal trajectories: RK4Tracer started from well spread seed points
//  * orthogonality check: finds where the two families cross and measures
//    the angle between the tangents (cos θ should be ≈ 0).
#pragma once
#include <array>
#include <vector>

#include "FamilyAnalyzer.h"
#include "Solver.h"
#include "Types.h"

namespace ot {

struct OrthoCheck {
    int crossings = 0;
    double maxAbsCos = 0;   // worst |cos θ| (0 means exactly 90°)
    double meanAbsCos = 0;
    int skippedSingular = 0;  // crossings at points where the slope is undefined
    // Every tested crossing, for right-angle marks on the graph:
    // x, y and the unit tangent (ux, uy) of the family curve, four numbers each.
    std::vector<double> marks;
};

class CurveGenerator {
public:
    // view: the visible area.   area: the (bigger) area we compute in.
    CurveGenerator(const FamilyAnalyzer& fa, bool isExplicit, Ptr explicitRhs, const Box& view, const Box& area);

    std::vector<Polyline> family(const std::vector<double>& cValues) const;
    std::vector<Polyline> trajectories(const std::vector<Polyline>& familyCurves, int count) const;
    OrthoCheck checkOrthogonality(const std::vector<Polyline>& fam, const std::vector<Polyline>& traj) const;

    double stepSize() const { return h_; }

private:
    const FamilyAnalyzer& fa_;
    bool isExplicit_;
    Ptr rhs_;
    Box view_, area_;
    double h_;  // RK4 step length

    void explicitCurve(double c, std::vector<Polyline>& out) const;
    void marchingSquares(double c, std::vector<Polyline>& out) const;
    std::vector<std::array<double, 4>> symmetries() const;
};

}  // namespace ot
