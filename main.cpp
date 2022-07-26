#include <iostream>
#include <sstream>
#include <stdio.h>
#include "Point.h"
#include "Point.cpp"
#include "Segment.h"
#include "Segment.cpp"
#include <algorithm>
#include <string>
#include <vector>
#include <queue>
#include <map>
#include <set>

using namespace std;
int dim,traj_n;
vector<vector<Point>> originLines;
vector<Segment*> allSegments;                     // 线段的集合
map<int,set<int>> cluster2Segments;           
vector<vector<Point*>> representive_point;
vector<vector<Segment*>> allClusteredSegment;
const double eplison = 40;
const int minLins = 3;
const int minLines = 1;
const int minClusterSegCnt = 4;
const double min_dist = 3;
const double constant = 1;


/*
 * 测试一下距离公式
 */

void test_dist() {
//    Point pi_s(1,1),pi_e(3,1);
//    Point pj_s(1,2),pj_e(3,4);
//    Segment si(pi_s,pi_e),sj(pj_s,pj_e);
//    cout << min(1,2);/**/
//    cout << si.perpendicular_dist(sj) << "  " << si.parallel_dist(sj) <<" "<<si.angle_dist(sj)<< endl;
}

// cost 计算函数
double MDL(int traj_id, int start_index, int currIndex, string t){
    Segment se = Segment(originLines[traj_id][start_index], originLines[traj_id][currIndex], traj_id);
    double l_h = 0.0, l_d_h_l = 0.0, l_d_h_r = 0.0;
    if(t == "par") {
        if(se.getSegmentLength() < eps) {
            l_h = 0;
        }
        else {
            l_h = log2(se.getSegmentLength());
        }
    }
    for (int i = start_index; i < currIndex; i++) {
        Segment sub = Segment(originLines[traj_id][i] , originLines[traj_id][i + 1], traj_id);
        if (t == "par") {
            l_d_h_l += se.perpendicular_dist(sub);
            // l_d_h_l += 0.5 * se.perpendicular_dist(sub);
            l_d_h_r += se.angle_dist(sub);
        } else if (t == "nopar"){
            l_h += sub.getSegmentLength();
        }
    }
    if (t == "par") {
        if (l_d_h_l > eps) {
            l_h += log2(l_d_h_l);
        }
        if (l_d_h_r > eps) {
            l_h += log2(l_d_h_r);
        }
        return l_h;
    } else {
        if (l_h < eps) {
            return 0;
        } else {
            return log2(l_h);
        }
    }
}

// 切分算法
void approximate_trajectory_partitioning() {
    int seg_id = 0;
    freopen("/home/nio/traclus_test/test_0722/car_segment_result.tra","w", stdout);
    for(int i = 0; i < originLines.size(); i++) {
        vector<Point> seq = originLines[i];
        size_t size = seq.size();
        int start_index = 0, length = 1;
        cout << start_index;
        while (start_index + length < size) {
            int currIndex = start_index + length;
            double cost_par = MDL(i, start_index, currIndex, "par");
            double cost_nopar = MDL(i, start_index, currIndex, "nopar");

            if (cost_par > cost_nopar + constant || length == size - 1) {
                cout << " " << currIndex - 1;
                Segment *seg = new Segment(originLines[i][start_index],originLines[i][currIndex - 1], i, seg_id++);
                allSegments.push_back(seg);
                start_index = currIndex - 1;
            } else {
                ++length;
            }
        }
        // conner case
        Segment *seg = new Segment(originLines[i][start_index], originLines[i][seq.size() - 1], i, seg_id++);
        cout << " " << seq.size() - 1 << endl;
        allSegments.push_back(seg);
    }
}

// 选取所有与输入线段距离小于eplison线段
set<Segment*> n_eplison_l(Segment *seg) {
    set<Segment*> tmp;
    Segment *seg_i, *seg_j;
    for (Segment *curSeg: allSegments) {
        if (seg == curSeg) {           // 同一线段则跳过
            continue;
        }
        if (curSeg->getSegmentLength() > seg->getSegmentLength()) {     // 取较长的
            seg_i = curSeg;
            seg_j = seg;
        } else {
            seg_i = seg;
            seg_j = curSeg;
        }
        if (seg_i->getAllDistance(*seg_j) <= eplison) {                 // 距离满足则聚类
//            cout <<"all distance: "<< seg_i->getAllDistance(*seg_j) << " "<<endl;
            tmp.insert(curSeg);
        }
    }
    return tmp;
}

// 寻找与输入线段距离较近但未被聚类的线段，并聚类
void expand_cluster(queue<Segment*> que, int cluster_id){
//    cout<< "start  -------- "<<endl;
    while (!que.empty()) {
        Segment *front = que.front();
        que.pop();
        set<Segment*> tmp_set = n_eplison_l(front);
        if (tmp_set.size() >= minLins) {
            for (Segment *seg: tmp_set) {
                if (seg->getClusterId() == -1) {
                    que.push(seg);
                    seg->setClusterId(cluster_id);
                }
            }
        }
    }
}

struct cmp {
    bool operator()(Point p1, Point p2) {
        if (p1.getX() != p2.getX() ){
            return p1.getX() < p2.getX();
        } else {
            return p1.getY() < p2.getY();
        }
    }
};

void representative_trajectory_generation(){
    representive_point.resize(allClusteredSegment.size());
    for (int i = 0; i < allClusteredSegment.size() ; i++) {
        int point_size = (int)allClusteredSegment[i].size() * 2;
        Point *sort_point = new Point[point_size];
        // vector<Point> sort_point(point_size);
        Point *rep_point = new Point(0, 0, -1), *zero_point = new Point(1, 0, -1);

        // average direction vector of current cluster 
        for (int j = 0; j < allClusteredSegment[i].size(); j++) {
            *rep_point = *rep_point + (allClusteredSegment[i][j]->getSegment().second - allClusteredSegment[i][j]->getSegment().first);
        }
        *rep_point = *rep_point / allClusteredSegment[i].size();

        // rotation angle
        double cos_theta = rep_point->dot(*zero_point) / rep_point->dist(Point(0,0,-1)); // cos(theta)
        double sin_theta = sqrt(1 - pow(cos_theta , 2)); // sin(theta)

        // rotate & sort
        for (int j = 0; j < allClusteredSegment[i].size(); j++) {
            Point s =  allClusteredSegment[i][j]->getSegment().first, e =  allClusteredSegment[i][j]->getSegment().second;
            allClusteredSegment[i][j]->setSegment(Point(s.getX() * cos_theta + s.getY() * sin_theta, s.getY() * cos_theta - s.getX() * sin_theta, -1),
                                                  Point(e.getX() * cos_theta + e.getY() * sin_theta, e.getY() * cos_theta - e.getX() * sin_theta, -1));
            sort_point[2 * j] = allClusteredSegment[i][j]->getSegment().first;
            sort_point[2 * j + 1] = allClusteredSegment[i][j]->getSegment().second;
        }
        // sort all start & end points in this cluster by x
        sort(sort_point, sort_point + point_size, cmp());
        // sort(sort_point.begin(), sort_point.end(),
        //     [](Point p1, Point p2)->bool{return p1.getX() < p2.getX();});

        for (int p = 0; p < point_size; p++) { //
            int intersect_cnt = 0;
            Point *y = new Point(0, 0 , -1);
 
            for(int q = 0 ;q < allClusteredSegment[i].size(); q++){
                // s: start point of segment in cluster, e: end point of segment in cluster
                Point s = allClusteredSegment[i][q]->getSegment().first, e = allClusteredSegment[i][q]->getSegment().second;
                if(sort_point[p].getX() <= max(e.getX(), s.getX()) && sort_point[p].getX() >= min(s.getX(),e.getX())){ // 在一个线段的中间
                    // 计算该点在线段上的点
                    if(s.getX() == e.getX()){
                        continue;
                    } else if (s.getY() == e.getY()) {
                        intersect_cnt++;
                        *y = *y + Point(sort_point[p].getX(), s.getY(), -1); //
                    } else {
                        intersect_cnt++;
                        *y = *y + Point(sort_point[p].getX(),(e.getY() - s.getY()) / (e.getX() - s.getX()) * (sort_point[p].getX() - s.getX()) + s.getY(), -1);
                    }
                }
            }

            if (intersect_cnt >= minLines) { // 取均值
                *y = *y / intersect_cnt;
                Point *tmp = new Point(y->getX() * cos_theta - sin_theta * y->getY(), sin_theta * y->getX() + cos_theta * y->getY(), -1);
                int size_point = representive_point[i].size() - 1; // 不让两点之间靠的太近
                // representive_point[i].push_back(tmp);
                if (size_point < 0 || (size_point >= 0 && tmp->dist(*representive_point[i][size_point]) > min_dist))
                    representive_point[i].push_back(tmp);
            }
            delete y;
        }
        delete []sort_point;
        delete rep_point, zero_point;
    }
    // output
    freopen("/home/nio/traclus_test/test_0722/reprezentitive.tra","w", stdout);
    for(int i = 0; i < representive_point.size(); i++) { 
        int len = ((int)representive_point[i].size()) - 1;
        for(int j = 0; j < len ;j++){
            cout << representive_point[i][j]->getX() << " "<< representive_point[i][j]->getY() << " ";
        }
        cout << endl;
    }
}

void line_segment_clustering(){
    freopen("/home/nio/traclus_test/test_0722/segment_cluster.tra","w", stdout);
    // cout << "allSegments_size: "<< allSegments.size() << endl;
    int cluster_id = 0;
    for (Segment *seg : allSegments) {
        queue<Segment*> que;
        while (!que.empty()) que.pop();
        if (seg->getClusterId() == -1) { // 该线段未被聚类
            // 寻找距离较近的线段
            set<Segment*> se = n_eplison_l(seg);
            // 判断数量是否满足
            if ((int)se.size() >= minLins) { // include itself
                // 设置id
                seg->setClusterId(cluster_id);
                for (Segment *seg1 : se) {
                    seg1->setClusterId(cluster_id);
                }
                // 扩大聚类
                expand_cluster(que, cluster_id++);
            }
        }
        if (seg->getClusterId() != -1) {    // 已被分类
            // cout << "segment_id : " << seg->getClusterId() << endl;
            cluster2Segments[seg->getClusterId()].insert(seg->getSegmentId());
        }
    }

    // for (const auto& seg:cluster2Segments) {
    //     cout << "ClusterId:" << seg.first << "\t"; 
    //     cout << "SegmentId:[";
    //     for (const auto& curSegment:seg.second)
    //         cout << curSegment << " ";
    //     cout << "]" << endl;
    // }

    set<int> remove_cluster;
    int cueClustersSize = (int)cluster2Segments.size();
    // cout << cueClustersSize <<" id  "<< cluster_id<< endl;
    for (int i = 0; i < cueClustersSize; i++) {
        if (cluster2Segments[i].size() < minClusterSegCnt) {
            remove_cluster.insert(i);
        }
    }

    // id + 1 ?
    allClusteredSegment.resize(cluster_id + 1);
    // remove cluster前cluster_id 映射至 remove后的cluster_id 
    int newClusterId = 0;
    map<int, int> preId2newId;
    preId2newId.clear();
    //remove
    try {
        for (Segment *curSeg : allSegments) {
            int curSegClusterId = curSeg->getClusterId();
            if (curSegClusterId == -1) {
                continue;
            } else if (remove_cluster.find(curSegClusterId) == remove_cluster.end()) {     // not in remove_cluster
                // get pre-cluster_id of current segment 
                // cur Id is not in pre-Cluster-Id to cur-Cluster-Id map
                if (preId2newId.find(curSegClusterId) == preId2newId.end()) {    
                    preId2newId.insert(make_pair(curSegClusterId, newClusterId));
                    allClusteredSegment[newClusterId].push_back(curSeg);
                    ++newClusterId;
                } else {
                    allClusteredSegment[preId2newId[curSegClusterId]].push_back(curSeg);
                }
            }
        }
        allClusteredSegment.resize(newClusterId);
        // output
        // cout << "newClusterSize:" << allClusteredSegment.size() << endl;
        for (size_t i = 0; i < allClusteredSegment.size(); ++i) {
            for (const auto& tmpSeg:allClusteredSegment[i]) {
                cout << i << " "
                     << tmpSeg->getSegment().first.getX() << " " << tmpSeg->getSegment().first.getY() << " "
                     << tmpSeg->getSegment().second.getX() << " " << tmpSeg->getSegment().second.getY() <<  " " << endl;
            }
        }
    }
    catch(exception e) {
        cout << " --------------------- " << endl;
        cout << e.what() <<endl;
    }
    fclose(stdout);
    representative_trajectory_generation();
}

int main() {
    // read
    // freopen("/home/nio/traclus_test/test_0722/debug_cluster_30.tra", "r", stdin);
    freopen("/home/nio/traclus_test/test_0722/origin_car_mapping.tra", "r", stdin);
    // freopen("/home/nio/traclus_test/test_0722/car_mapping_2mesh.tra", "r", stdin);
    // freopen("/home/nio/traclus_test/test_0722/debug_18.tra", "r", stdin);
    // freopen("/home/nio/traclus_test/test_0720/car_557039869_result.tra","w", stdout);

    string strLine;
    int idx = 0;
    while (getline(cin, strLine)) {
        if (strLine.size() == 0)
            continue;
        istringstream record(strLine);
        vector<Point> curLine;
        double x, y;
        int sample = 0;
        while (record >> x >> y) {
            if (sample % 3 == 0) {
                curLine.emplace_back(x, y, idx);
            } else {
                int a = x;
                int b = y;
            }
        }
        originLines.push_back(curLine);
        ++idx;
    }

    // debug representative lines
    // vector<Segment*> curSegments;
    // while (getline(cin, strLine)) {
    //     if (strLine.size() == 0)
    //         continue;
    //     istringstream record(strLine);
    //     double x1, y1, x2, y2;
    //     while(record >> x1 >> y1 >> x2 >> y2) {
    //         Segment* curSeg = new Segment(Point(x1, y1, -1), Point(x2, y2, -1), -1, -1);
    //         curSegments.emplace_back(curSeg);
    //     }
    // }
    // allClusteredSegment.push_back(curSegments);
    
    fclose(stdin);
    approximate_trajectory_partitioning();
    line_segment_clustering();

    return 0;
}