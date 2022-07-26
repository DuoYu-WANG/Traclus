//
// Created by 张开顺 on 2018/12/17.
//

#ifndef TRACLUS_SEGMENT_H
#define TRACLUS_SEGMENT_H
#include "Point.h"
#include <cstdio>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <math.h>
class Segment {
private:
    Point s, e;
    const int traj_id;
    int seg_id;
    int cluster_id;

public:
    Segment() = delete;
    Segment(const Point& p1, const Point& p2, const int& traj_id, const int& seg_id = -2);
    std::pair<Point,Point> getSegment();
    void const setClusterId(int cluster_id);
    int const getClusterId();
    int const getTrajId();
    int const getSegmentId();
    double const perpendicular_dist(Segment seg); // 两个线段的垂直距离
    double const parallel_dist(Segment seg);
    double const angle_dist(Segment seg);
    void const setSegment(Point p,Point q);
    double getSegmentLength();
    double const getAllDistance(Segment seg);
};
#endif //TRACLUS_SEGMENT_H
