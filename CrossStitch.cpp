#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <random>
#include <algorithm>
#include <utility>
#include <limits>
#include <queue>
#include <cmath>

#define NUMBER_OF_LETTERS 26
#define LOCAL_INFINITY std::numeric_limits<double>::max() // 2^31 -1
#define ABSOLUTE_TEMPERATURE 1.0
#define LOCAL_EPSILON 1e-6


using namespace std;

typedef void (*function_pointer)(vector<int>&, vector<int>&, int, int);

void tr_bl_br_tl(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y);
void tr_bl_tl_br(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y);
void tl_br_bl_tr(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y);
void tl_br_tr_bl(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y);
void br_tl_tr_bl(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y);
void br_tl_bl_tr(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y);
void bl_tr_tl_br(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y);
void bl_tr_br_tl(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y);


inline int get_letter_index(char &letter) {
    return letter - 97;
}


inline char revert_get_letter_index(int &letter_index){
    return letter_index + 97;
}


class CPCC {
    
    // CPCF - Custom Pair Comparison Class

    public:
        bool operator() (const pair<double, int> &left, const pair<double, int> &right) {
            return left.first > right.first;
        }
};

//-------------------------------------------------------
//--------------------Point------------------------------
//-------------------------------------------------------


class Point {
    private:
        int m_x;
        int m_y;
        int m_b; // # badge
        double m_fd; // forward distance to point
        double m_bd; // backward distance to point

    public:
        Point() { m_x = -1.0; m_y = -1.0; m_b = -1; m_fd = -1.0; m_bd = -1.0; }
        Point(int x, int y): m_x(x), m_y(y) { m_b = -1; m_fd = -1.0; m_bd = -1.0; }        
        Point(int x, int y, int b): m_x(x), m_y(y), m_b(b) { m_fd = -1.0; m_bd = -1.0; }
        Point(int x, int y, int b, double fd, double bd): m_x(x), m_y(y), m_b(b), m_fd(fd), m_bd(bd) {}

        void set_x(int x) {m_x = x;}
        void set_y(int y) {m_y = y;}
        void set_b(int b) {m_b = b;}

        int get_x() const {return m_x;}
        int get_y() const {return m_y;}
        int get_b() const {return m_b;}

        bool operator==(const Point &other);

        static double euclidean_distance(Point &p1, Point &p2);
        static int manhattan_distance(Point &p1, Point &p2);

};


bool Point::operator==(const Point &other) {

    if ( (m_x == other.m_x) && (m_y == other.m_y) )
        return true;
    else
        return false;
}

double Point::euclidean_distance(Point &p1, Point &p2) {

    double x1 = (double)p1.get_x();
    double y1 = (double)p1.get_y();

    double x2 = (double)p2.get_x();
    double y2 = (double)p2.get_y();

    double d = sqrt( (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) );

    return d;
}


int Point::manhattan_distance(Point &p1, Point &p2) {

    int x1 = p1.get_x();
    int y1 = p1.get_y();

    int x2 = p2.get_x();
    int y2 = p2.get_y();

    int d = abs(x1 - x2) + abs(y1 - y2);

    return d;
}


std::ostream& operator<<(std::ostream& os, const Point& p) {

	string sx = to_string(p.get_x());
    string sy = to_string(p.get_y());
    string sb = to_string(p.get_b());

	string s = "x: " + sx + " y: " + sy + " b: " + sb;

    os << s;
    return os;
}



//-------------------------------------------------------
//---------------PointCollection-------------------------
//-------------------------------------------------------


class PointCollection {
    private:

        default_random_engine m_engine;

        vector<vector<Point>> m_point_collection;
        vector< vector<vector<double>>> m_distance_matrix;

        vector<vector<double>> get_distance_matrix(vector<Point> &path);

    public:
        PointCollection() {}
        PointCollection(vector<string> &pattern);


        void eprint_collection();
        void eprint_number_of_points_in_each_path();
        void eprint_distance_matrixes(int letter_index);
        void eprint_average_path_length_of_collection();

        double get_path_length_in_collection(vector<Point> &path, int &letter_index);

        vector<Point> nearest_neighbour_single_path_optimize(vector<Point> &path, int &letter_index);
        void nearest_neighbour_all_path_optimize();

        void minimal_spannig_tree(vector<Point> &path, int &letter_index, vector<vector<vector<double>>> &distance_matrix);
        void hungarian(vector<vector<double>> &cost_matrix, vector<int> &left_match, vector<int> &right_match);

        void two_opt_single_path_optimize(vector<Point> &path, int &letter_index);
        void two_opt_all_path_optimize();


        void point_swap(int &&p1i, int &&p2i, vector<Point> &path);
        void get_swap_points(int &n_points, pair<double, double> &p);
        vector<Point> markov_monte_carlo_like_single_path_optimize(int &n_iterations, vector<Point> path, int &letter_index);
        void markov_monte_carlo_like_all_path_optimize();

        double calculate_stitch_d_distance(vector<Point> &pb_points, vector<Point> &pf_points);
        bool check_if_stitch_order_is_valid(vector<Point> &pb_points, vector<Point> &pf_points);
        void product_stitch(vector<string> &ret, vector<Point> path, int depth);
        vector<string> make_stitchs(int depth);

        void rotate_collection();

};


vector<vector<double>> PointCollection::get_distance_matrix(vector<Point> &path){

    int n_points = path.size();
    vector<vector<double>> M = vector<vector<double>>(n_points);

    for(int i = 0; i < n_points; i++)
        M[i].resize(n_points, 0.0);

    for(int i = 0; i < n_points; i++){
        for(int j = 0; j < n_points; j++){
            if (i <= j)
                continue;            
            double d = Point::euclidean_distance(path[i], path[j]);
            //double d = Point::manhattan_distance(path[i], path[j]);

            M[i][j] = d;
            M[j][i] = d;
        }
    }

    return M;
}



PointCollection::PointCollection(vector<string> &pattern){

    random_device rd;
    m_engine.seed(rd());

    m_point_collection = vector<vector<Point>>(NUMBER_OF_LETTERS);
    m_distance_matrix = vector< vector<vector<double>>>(NUMBER_OF_LETTERS);

    int S = pattern.size();
    
    for(int i = 0; i < S; i++) {
        for(int j = 0; j < S; j++){

            if (pattern[i][j] != '.') {
                int letter_index = get_letter_index(pattern[i][j]);
                m_point_collection[letter_index].push_back(Point(i,j));
            }
        }
    }

    // Badge the point in each path.
    for(int letter_index = 0; letter_index < NUMBER_OF_LETTERS; letter_index++) {
        int n_points = m_point_collection[letter_index].size();

        if (n_points == 0)
            continue;

        for(int i = 0; i < n_points; i++)
            m_point_collection[letter_index][i].set_b(i);


        m_distance_matrix[letter_index] = get_distance_matrix(m_point_collection[letter_index]);
    }

}


void PointCollection::eprint_collection(){

    for(int letter_index = 0; letter_index < NUMBER_OF_LETTERS; letter_index++){
        int n_points = m_point_collection[letter_index].size();
        cerr << "Number of " << (char)revert_get_letter_index(letter_index) << " letters: " << n_points << endl;
        for(auto p : m_point_collection[letter_index])
            cerr << p << endl;
    }
}


void PointCollection::eprint_number_of_points_in_each_path(){

    for(int letter_index = 0; letter_index < NUMBER_OF_LETTERS; letter_index++){
        int n_points = m_point_collection[letter_index].size();
        if (n_points == 0)
            continue;
        cerr << "Number of " << (char)revert_get_letter_index(letter_index) << " letters: " << n_points << endl;
    }

}


void PointCollection::eprint_distance_matrixes(int letter_index){
    // 0 <= letter_index <= 25

    for(int i = 0; i < m_distance_matrix[letter_index].size(); i++){
        for(int j = 0; j < m_distance_matrix[letter_index].size(); j++)
            fprintf(stderr, "%2.2f ", m_distance_matrix[letter_index][i][j]);

        cerr << endl;
    }

}


void PointCollection::eprint_average_path_length_of_collection(){

    double apl = 0.0; // Average Path Length
    double nvp = 0; // number of valid paths

    for(int letter_index = 0; letter_index < NUMBER_OF_LETTERS; letter_index++){
        int n_points = m_point_collection[letter_index].size();
        if (n_points == 0)
            continue;

        apl = apl + get_path_length_in_collection(m_point_collection[letter_index], letter_index);
        nvp = nvp + 1;
    }
    cerr << "\nAverage path length: " << apl/nvp << endl << endl;
}


double PointCollection::get_path_length_in_collection(vector<Point> &path, int &letter_index){

    int n_points = path.size();

    double pl = 0.0; // path length
    for(int i = 1; i < n_points; i++){

        int b1 = path[i].get_b();
        int b2 = path[i - 1].get_b();
        pl = pl + m_distance_matrix[letter_index][b1][b2];
    }

    return pl;
}


vector<Point> PointCollection::nearest_neighbour_single_path_optimize(vector<Point> &path, int &letter_index){

    int n_points = path.size();

    vector<Point> minimized_path(n_points);
    vector<bool> visited_points(n_points, false);

    uniform_int_distribution<int> dist(0, n_points - 1);
    int cpi = dist(m_engine); // Current Point Index
    int nocp = 0; // Number Of Checked Points
 
    Point *cp = &path[cpi];
    minimized_path[nocp] = (*cp);
    visited_points[cpi] = true;

    while (nocp != n_points - 1) {

        double min_d = LOCAL_INFINITY;
        for(int i = 0; i < n_points; i++) {
            if (visited_points[i] == true)
                continue;

            int cpb = cp->get_b();
            int p2b = path[i].get_b();

            double d = m_distance_matrix[letter_index][cpb][p2b];

            if (d < min_d){
                min_d = d;
                cpi = i;
            }
        }
        cp = &path[cpi];
        visited_points[cpi] = true;
        nocp = nocp + 1;
        minimized_path[nocp] = (*cp);
    }

    return minimized_path;
}


void PointCollection::nearest_neighbour_all_path_optimize(){

    for(int letter_index = 0; letter_index < NUMBER_OF_LETTERS; letter_index++){
        if (m_point_collection[letter_index].size() == 0)
            continue;

        m_point_collection[letter_index] = nearest_neighbour_single_path_optimize(m_point_collection[letter_index], letter_index);
    }
}


void PointCollection::minimal_spannig_tree(vector<Point> &path, int &letter_index, vector<vector<vector<double>>> &distance_matrix){

    // Construct a minimal spanning tree for a single path.

    int n_points = path.size();
    
    vector<double> key(n_points, LOCAL_INFINITY);
    vector<int> mst(n_points, -1);
    vector<vector<int>> mst_adj_list(n_points);
    vector<bool> present_in_mst(n_points, false);

    priority_queue<pair<double, int>, vector<pair<double, int>>, CPCC> pq;

    uniform_int_distribution<int> random_initial_key(0, n_points - 1);
    int initial_key = random_initial_key(m_engine);

    pq.push(make_pair(0.0, initial_key));
    key[initial_key] = 0;

    while( !pq.empty() ) {

        int u = pq.top().second;
        pq.pop();

        cerr << "u: " << u << endl;

        present_in_mst[u] = true;

        for(int i = 0; i < n_points; i++){
            
            if (i == u)
                continue;

            double distance = distance_matrix[letter_index][u][i];

            if (distance < 0.001)
                continue;

            if(present_in_mst[i] == false && key[i] > distance){
                key[i] = distance;
                pq.push(make_pair(distance, i));
                mst[i] = u;
            }
        }
    }

    cerr << "Minimal spanning tree" << endl;
    for (int i = 0; i < n_points; ++i){
        if(mst[i] == -1)
            continue;
        mst_adj_list[i].push_back(mst[i]);
        mst_adj_list[mst[i]].push_back(i);
        cerr << mst[i] << " - " << i << endl;
    }
    cerr << "===V2===" << endl;

    for(int i = 0; i < n_points; i++){
        cerr << "i = " << i << " : ";
        for(int j = 0; j < mst_adj_list[i].size(); j++){
            cerr << mst_adj_list[i][j] << " - " ;
        }
        cerr << endl;
    }

    vector<int> odd_points;
    for(int i = 0; i < n_points; i++){
        if (mst_adj_list[i].size() % 2  != 0)
            odd_points.push_back(i);
    }

    int n_odd_points = odd_points.size();
    vector<vector<double>> cost_matrix(n_odd_points, vector<double>(n_odd_points, 0.0));

    for(int i = 0; i < n_odd_points; i++) {
        int di = odd_points[i];
        for(int j = 0; j < n_odd_points; j++) {
            int dj = odd_points[j];
            cost_matrix[i][j] = distance_matrix[letter_index][di][dj];
        }
    }

    cerr << "--- Cost matrix ---" << endl;
    for(int i = 0; i < n_odd_points; i++){
        for(int j = 0; j < n_odd_points; j++){
            cerr << cost_matrix[i][j] << " ";
        }
        cerr << endl;
    }
    cerr << "-------------------" << endl;


    // Minimum-Cost Perfect Matching

    vector<int> left_match;
    vector<int> right_match;

    hungarian(cost_matrix, left_match, right_match);
    //double ddd = MinCostMatching(cost_matrix, left_match, right_match);
    //cerr << "ddd: " << ddd << endl;


    // Print Matching
    for(int i = 0; i < left_match.size(); i++){
        int oi = odd_points[left_match[i]];
        int oj = odd_points[right_match[i]];

        cerr << oi << " <-> " << oj << endl;
    }

}



void PointCollection::hungarian(vector<vector<double>> &cost_matrix, vector<int> &left_match, vector<int> &right_match){

    int n_jobs = cost_matrix.size();

    vector<double> u_labels(n_jobs, 0.0);
    vector<double> v_labels(n_jobs, 0.0);

    for(int i = 0; i < n_jobs; i++){
        u_labels[i] = cost_matrix[i][0];
        for(int j = 0; j < n_jobs; j++)
            u_labels[i] = min(u_labels[i], cost_matrix[i][j]);
    }

    for(int i = 0; i < n_jobs; i++){
        v_labels[i] = cost_matrix[0][i] - u_labels[0];
        for(int j = 0; j < n_jobs; j++)
            v_labels[i] = min(v_labels[i], cost_matrix[j][i] - u_labels[j]);
    }

    left_match = vector<int>(n_jobs, -1);
    right_match = vector<int>(n_jobs, -1);

    int matched_up_to_now = 0;

    for(int i = 0; i < n_jobs; i++){
        for(int j = 0; j < n_jobs; j++){

            if (right_match[j] != -1)
                continue;

            if ( fabs(cost_matrix[i][j] - u_labels[i] - v_labels[j]) <  LOCAL_EPSILON){
                
                left_match[i] = j;
                right_match[j] = i;
                matched_up_to_now++;
                break;
            }
        }
    }

    vector<double> dijkstra(n_jobs, 0.0);
    vector<int> dad(n_jobs, -1);
    vector<bool> seen(n_jobs, false);

    while(matched_up_to_now < n_jobs) {

        int s_index = 0;
        while(left_match[s_index] != -1)
            s_index++;

        for(int k = 0; k < n_jobs; k++)
            dijkstra[k] = cost_matrix[s_index][k] - u_labels[s_index] - v_labels[k];

        int j_index = 0;
        while(true) {

            j_index = -1;
            for(int k = 0; k < n_jobs; k++){
                if (seen[k])
                    continue;

                if (j_index == -1 || dijkstra[k] < dijkstra[j_index])
                    j_index = k;
            }

            seen[j_index] = 1;

            if (right_match[j_index] == -1)
                break;

            const int i_index = right_match[j_index];
            for(int k = 0; k < n_jobs; k++){
                if (seen[k])
                    continue;
                
                const double new_djikstra = dijkstra[j_index] - cost_matrix[i_index][k] - u_labels[i_index] - v_labels[k];
                if (dijkstra[k] > new_djikstra){
                    dijkstra[k] = new_djikstra;
                    dad[k] = j_index;
                }

            }
        } // end of inner while loop
        
        for(int k = 0; k < n_jobs; k++){
            if (k == j_index || !seen[k])
                continue;

            const int i = right_match[k];
            u_labels[i] = u_labels[k] - dijkstra[k] - dijkstra[j_index];
            v_labels[k] = v_labels[k] + dijkstra[k] - dijkstra[j_index];
        }

        u_labels[s_index] = u_labels[s_index] + dijkstra[j_index];

        while(dad[j_index] >= 0) {

            const int d = dad[j_index];
            right_match[j_index] = right_match[d];
            left_match[right_match[j_index]] = j_index;

            j_index = d;
            //cerr << "While dad" << endl;
        } // end of inner while loop

        right_match[j_index] = s_index;
        left_match[s_index] = j_index;
        matched_up_to_now++;

    }
    
    cerr << "--- Cost matrix ---" << endl;
    for(int i = 0; i < n_jobs; i++){
        for(int j = 0; j < n_jobs; j++){
            cerr << cost_matrix[i][j] << " ";
        }
        cerr << endl;
    }

    cerr << "left_match: ";
    for(int i = 0; i < n_jobs; i++){
        cerr << left_match[i] << " ";
    }
    cerr << endl;

    double value = 0.0;
    for(int i = 0; i < n_jobs; i++){
        cerr << "value: " << value << endl;
        value = value + cost_matrix[right_match[i]][left_match[i]];
    }

    cerr << "Min cost perfect matching: " << value << endl;

}




void PointCollection::two_opt_single_path_optimize(vector<Point> &path, int &letter_index){

    int imp = 0;
    int n_points = path.size();

    while (imp < 1){

        for(int i = 1; i < n_points - 1; i++){
            for(int j = i + 1; j < n_points; j++){

                double d_before_prev = 0.0;
                double d_before_post = 0.0;

                double d_after_prev = 0.0;
                double d_after_post = 0.0;


                if (j != (n_points - 1)){
                    int i1b = path[i - 1].get_b();
                    int ib = path[i].get_b();
                    d_before_prev = m_distance_matrix[letter_index][i1b][ib];

                    int jb = path[j].get_b(); 
                    int j1b = path[j + 1].get_b();
                    d_before_post = m_distance_matrix[letter_index][jb][j1b];


                    i1b = path[i - 1].get_b();
                    jb = path[j].get_b();
                    d_after_prev = m_distance_matrix[letter_index][i1b][jb];

                    ib = path[i].get_b();
                    j1b = path[j + 1].get_b();
                    d_after_post = m_distance_matrix[letter_index][ib][j1b];

                } else {
                    int i1b = path[i - 1].get_b();
                    int ib = path[i].get_b();
                    d_before_prev = m_distance_matrix[letter_index][i1b][ib];

                    i1b = path[i - 1].get_b();
                    int jb = path[j].get_b();
                    d_after_prev = m_distance_matrix[letter_index][i1b][jb];
                }

                double d1 = d_before_prev + d_before_post;
                double d2 = d_after_prev + d_after_post;

                if ((d1 - d2) > 0.0){
                        imp = 0;
                        reverse(path.begin() + i, path.begin() + j + 1 );
                }
            } // second for loop end
        } // first for loop end

        imp++;
    } // while loop end


}


void PointCollection::two_opt_all_path_optimize(){

    for(int letter_index = 0; letter_index < NUMBER_OF_LETTERS; letter_index++){
        if (m_point_collection[letter_index].size() == 0)
            continue;

        two_opt_single_path_optimize(m_point_collection[letter_index], letter_index);
    }
}


void PointCollection::point_swap(int &&p1i, int &&p2i, vector<Point> &path){

    int x = path[p1i].get_x();
    int y = path[p1i].get_y();
    int b = path[p1i].get_b();

    path[p1i].set_x( path[p2i].get_x() );
    path[p1i].set_y( path[p2i].get_y() );
    path[p1i].set_b( path[p2i].get_b() );

    path[p2i].set_x( x );
    path[p2i].set_y( y );
    path[p2i].set_b( b );

}

void PointCollection::get_swap_points(int &n_points, pair<double, double> &p){

    uniform_int_distribution<int> dist(0, n_points - 1);
    int p1i = dist(m_engine);
    int p2i = 0;

    p.first = p1i;

    if (p1i == 0) {
        p2i = p1i + 1;
        p.second = p2i;
        return;
    } else if (p1i == (n_points - 1)) {
        p2i = p1i - 1;
        p.second = p2i;
        return;
    } else {
        uniform_int_distribution<int> coin(0, 1);
        int coin_toss = coin(m_engine);
        if (coin_toss == 0) {
            p2i = p1i - 1;
            p.second = p2i;
            return;
        } else {
            p2i = p1i + 1;
            p.second = p2i;
            return;
        }
    }

    return;
}


vector<Point> PointCollection::markov_monte_carlo_like_single_path_optimize(int &n_iterations, vector<Point> path, int &letter_index){

    uniform_real_distribution<double> uni(0.0, 1.0);

    int n_points = path.size();

    pair<double, double> pr(0, 0);
    double min_energy = get_path_length_in_collection(path, letter_index);
    double c_energy = min_energy;

    double T = 1.0;
    double lambda = 0.9;

    vector<Point> new_path(path);

    for(int i = 1; i < n_iterations; i++){

        get_swap_points(n_points, pr);
        point_swap(pr.first, pr.second, path);

        double n_energy = get_path_length_in_collection(path, letter_index);

        double A = exp( (c_energy - n_energy) / T );
        double p = uni(m_engine);

        if (p < A){
            c_energy = n_energy;
            if (n_energy < min_energy){
                min_energy = n_energy;           
                new_path = path; 
            }
            continue;        
        } else {
            point_swap(pr.first, pr.second, path);
        }
    }

    return new_path;
}


void PointCollection::markov_monte_carlo_like_all_path_optimize(){

    int n_iterations = 1000;
    for(int letter_index = 0; letter_index < NUMBER_OF_LETTERS; letter_index++){
        if (m_point_collection[letter_index].size() == 0)
            continue;

        m_point_collection[letter_index] = markov_monte_carlo_like_single_path_optimize(n_iterations, m_point_collection[letter_index], letter_index);
    }

}


double PointCollection::calculate_stitch_d_distance(vector<Point> &pb_points, vector<Point> &pf_points){

    // Calcualtes stitch_d = sum( distance(pb_points[n], pf_points[n-1]) ), from n = 1 to len(pb_points)

    double stitch_d = 0; 
    for(int i = 1; i < pb_points.size(); i++)
        stitch_d = stitch_d + Point::euclidean_distance(pb_points[i], pf_points[i - 1]);

    return stitch_d;
}


bool PointCollection::check_if_stitch_order_is_valid(vector<Point> &pb_points, vector<Point> &pf_points){

    // Check if the begining point of stitch i is not equal
    // to the end of stitch (i - 1).

    for(int i = 1; i < pb_points.size(); i++)
        if (pb_points[i] == pf_points[i - 1])
            return false;

    return true;
}


void PointCollection::product_stitch(vector<string> &ret, vector<Point> path, int depth){

    vector< vector<function_pointer> > stitich_functions = {
                                                            {tr_bl_br_tl,tr_bl_br_tl,tr_bl_br_tl},
                                                            {tr_bl_br_tl,tr_bl_br_tl,tr_bl_tl_br},
                                                            {tr_bl_br_tl,tr_bl_br_tl,tl_br_bl_tr},
                                                            {tr_bl_br_tl,tr_bl_br_tl,tl_br_tr_bl},
                                                            {tr_bl_br_tl,tr_bl_br_tl,br_tl_tr_bl},
                                                            {tr_bl_br_tl,tr_bl_br_tl,br_tl_bl_tr},
                                                            {tr_bl_br_tl,tr_bl_br_tl,bl_tr_tl_br},
                                                            {tr_bl_br_tl,tr_bl_br_tl,bl_tr_br_tl},
                                                            {tr_bl_br_tl,tr_bl_tl_br,tr_bl_br_tl},
                                                            {tr_bl_br_tl,tr_bl_tl_br,tr_bl_tl_br},
                                                            {tr_bl_br_tl,tr_bl_tl_br,tl_br_bl_tr},
                                                            {tr_bl_br_tl,tr_bl_tl_br,tl_br_tr_bl},
                                                            {tr_bl_br_tl,tr_bl_tl_br,br_tl_tr_bl},
                                                            {tr_bl_br_tl,tr_bl_tl_br,br_tl_bl_tr},
                                                            {tr_bl_br_tl,tr_bl_tl_br,bl_tr_tl_br},
                                                            {tr_bl_br_tl,tr_bl_tl_br,bl_tr_br_tl},
                                                            {tr_bl_br_tl,tl_br_bl_tr,tr_bl_br_tl},
                                                            {tr_bl_br_tl,tl_br_bl_tr,tr_bl_tl_br},
                                                            {tr_bl_br_tl,tl_br_bl_tr,tl_br_bl_tr},
                                                            {tr_bl_br_tl,tl_br_bl_tr,tl_br_tr_bl},
                                                            {tr_bl_br_tl,tl_br_bl_tr,br_tl_tr_bl},
                                                            {tr_bl_br_tl,tl_br_bl_tr,br_tl_bl_tr},
                                                            {tr_bl_br_tl,tl_br_bl_tr,bl_tr_tl_br},
                                                            {tr_bl_br_tl,tl_br_bl_tr,bl_tr_br_tl},
                                                            {tr_bl_br_tl,tl_br_tr_bl,tr_bl_br_tl},
                                                            {tr_bl_br_tl,tl_br_tr_bl,tr_bl_tl_br},
                                                            {tr_bl_br_tl,tl_br_tr_bl,tl_br_bl_tr},
                                                            {tr_bl_br_tl,tl_br_tr_bl,tl_br_tr_bl},
                                                            {tr_bl_br_tl,tl_br_tr_bl,br_tl_tr_bl},
                                                            {tr_bl_br_tl,tl_br_tr_bl,br_tl_bl_tr},
                                                            {tr_bl_br_tl,tl_br_tr_bl,bl_tr_tl_br},
                                                            {tr_bl_br_tl,tl_br_tr_bl,bl_tr_br_tl},
                                                            {tr_bl_br_tl,br_tl_tr_bl,tr_bl_br_tl},
                                                            {tr_bl_br_tl,br_tl_tr_bl,tr_bl_tl_br},
                                                            {tr_bl_br_tl,br_tl_tr_bl,tl_br_bl_tr},
                                                            {tr_bl_br_tl,br_tl_tr_bl,tl_br_tr_bl},
                                                            {tr_bl_br_tl,br_tl_tr_bl,br_tl_tr_bl},
                                                            {tr_bl_br_tl,br_tl_tr_bl,br_tl_bl_tr},
                                                            {tr_bl_br_tl,br_tl_tr_bl,bl_tr_tl_br},
                                                            {tr_bl_br_tl,br_tl_tr_bl,bl_tr_br_tl},
                                                            {tr_bl_br_tl,br_tl_bl_tr,tr_bl_br_tl},
                                                            {tr_bl_br_tl,br_tl_bl_tr,tr_bl_tl_br},
                                                            {tr_bl_br_tl,br_tl_bl_tr,tl_br_bl_tr},
                                                            {tr_bl_br_tl,br_tl_bl_tr,tl_br_tr_bl},
                                                            {tr_bl_br_tl,br_tl_bl_tr,br_tl_tr_bl},
                                                            {tr_bl_br_tl,br_tl_bl_tr,br_tl_bl_tr},
                                                            {tr_bl_br_tl,br_tl_bl_tr,bl_tr_tl_br},
                                                            {tr_bl_br_tl,br_tl_bl_tr,bl_tr_br_tl},
                                                            {tr_bl_br_tl,bl_tr_tl_br,tr_bl_br_tl},
                                                            {tr_bl_br_tl,bl_tr_tl_br,tr_bl_tl_br},
                                                            {tr_bl_br_tl,bl_tr_tl_br,tl_br_bl_tr},
                                                            {tr_bl_br_tl,bl_tr_tl_br,tl_br_tr_bl},
                                                            {tr_bl_br_tl,bl_tr_tl_br,br_tl_tr_bl},
                                                            {tr_bl_br_tl,bl_tr_tl_br,br_tl_bl_tr},
                                                            {tr_bl_br_tl,bl_tr_tl_br,bl_tr_tl_br},
                                                            {tr_bl_br_tl,bl_tr_tl_br,bl_tr_br_tl},
                                                            {tr_bl_br_tl,bl_tr_br_tl,tr_bl_br_tl},
                                                            {tr_bl_br_tl,bl_tr_br_tl,tr_bl_tl_br},
                                                            {tr_bl_br_tl,bl_tr_br_tl,tl_br_bl_tr},
                                                            {tr_bl_br_tl,bl_tr_br_tl,tl_br_tr_bl},
                                                            {tr_bl_br_tl,bl_tr_br_tl,br_tl_tr_bl},
                                                            {tr_bl_br_tl,bl_tr_br_tl,br_tl_bl_tr},
                                                            {tr_bl_br_tl,bl_tr_br_tl,bl_tr_tl_br},
                                                            {tr_bl_br_tl,bl_tr_br_tl,bl_tr_br_tl},
                                                            {tr_bl_tl_br,tr_bl_br_tl,tr_bl_br_tl},
                                                            {tr_bl_tl_br,tr_bl_br_tl,tr_bl_tl_br},
                                                            {tr_bl_tl_br,tr_bl_br_tl,tl_br_bl_tr},
                                                            {tr_bl_tl_br,tr_bl_br_tl,tl_br_tr_bl},
                                                            {tr_bl_tl_br,tr_bl_br_tl,br_tl_tr_bl},
                                                            {tr_bl_tl_br,tr_bl_br_tl,br_tl_bl_tr},
                                                            {tr_bl_tl_br,tr_bl_br_tl,bl_tr_tl_br},
                                                            {tr_bl_tl_br,tr_bl_br_tl,bl_tr_br_tl},
                                                            {tr_bl_tl_br,tr_bl_tl_br,tr_bl_br_tl},
                                                            {tr_bl_tl_br,tr_bl_tl_br,tr_bl_tl_br},
                                                            {tr_bl_tl_br,tr_bl_tl_br,tl_br_bl_tr},
                                                            {tr_bl_tl_br,tr_bl_tl_br,tl_br_tr_bl},
                                                            {tr_bl_tl_br,tr_bl_tl_br,br_tl_tr_bl},
                                                            {tr_bl_tl_br,tr_bl_tl_br,br_tl_bl_tr},
                                                            {tr_bl_tl_br,tr_bl_tl_br,bl_tr_tl_br},
                                                            {tr_bl_tl_br,tr_bl_tl_br,bl_tr_br_tl},
                                                            {tr_bl_tl_br,tl_br_bl_tr,tr_bl_br_tl},
                                                            {tr_bl_tl_br,tl_br_bl_tr,tr_bl_tl_br},
                                                            {tr_bl_tl_br,tl_br_bl_tr,tl_br_bl_tr},
                                                            {tr_bl_tl_br,tl_br_bl_tr,tl_br_tr_bl},
                                                            {tr_bl_tl_br,tl_br_bl_tr,br_tl_tr_bl},
                                                            {tr_bl_tl_br,tl_br_bl_tr,br_tl_bl_tr},
                                                            {tr_bl_tl_br,tl_br_bl_tr,bl_tr_tl_br},
                                                            {tr_bl_tl_br,tl_br_bl_tr,bl_tr_br_tl},
                                                            {tr_bl_tl_br,tl_br_tr_bl,tr_bl_br_tl},
                                                            {tr_bl_tl_br,tl_br_tr_bl,tr_bl_tl_br},
                                                            {tr_bl_tl_br,tl_br_tr_bl,tl_br_bl_tr},
                                                            {tr_bl_tl_br,tl_br_tr_bl,tl_br_tr_bl},
                                                            {tr_bl_tl_br,tl_br_tr_bl,br_tl_tr_bl},
                                                            {tr_bl_tl_br,tl_br_tr_bl,br_tl_bl_tr},
                                                            {tr_bl_tl_br,tl_br_tr_bl,bl_tr_tl_br},
                                                            {tr_bl_tl_br,tl_br_tr_bl,bl_tr_br_tl},
                                                            {tr_bl_tl_br,br_tl_tr_bl,tr_bl_br_tl},
                                                            {tr_bl_tl_br,br_tl_tr_bl,tr_bl_tl_br},
                                                            {tr_bl_tl_br,br_tl_tr_bl,tl_br_bl_tr},
                                                            {tr_bl_tl_br,br_tl_tr_bl,tl_br_tr_bl},
                                                            {tr_bl_tl_br,br_tl_tr_bl,br_tl_tr_bl},
                                                            {tr_bl_tl_br,br_tl_tr_bl,br_tl_bl_tr},
                                                            {tr_bl_tl_br,br_tl_tr_bl,bl_tr_tl_br},
                                                            {tr_bl_tl_br,br_tl_tr_bl,bl_tr_br_tl},
                                                            {tr_bl_tl_br,br_tl_bl_tr,tr_bl_br_tl},
                                                            {tr_bl_tl_br,br_tl_bl_tr,tr_bl_tl_br},
                                                            {tr_bl_tl_br,br_tl_bl_tr,tl_br_bl_tr},
                                                            {tr_bl_tl_br,br_tl_bl_tr,tl_br_tr_bl},
                                                            {tr_bl_tl_br,br_tl_bl_tr,br_tl_tr_bl},
                                                            {tr_bl_tl_br,br_tl_bl_tr,br_tl_bl_tr},
                                                            {tr_bl_tl_br,br_tl_bl_tr,bl_tr_tl_br},
                                                            {tr_bl_tl_br,br_tl_bl_tr,bl_tr_br_tl},
                                                            {tr_bl_tl_br,bl_tr_tl_br,tr_bl_br_tl},
                                                            {tr_bl_tl_br,bl_tr_tl_br,tr_bl_tl_br},
                                                            {tr_bl_tl_br,bl_tr_tl_br,tl_br_bl_tr},
                                                            {tr_bl_tl_br,bl_tr_tl_br,tl_br_tr_bl},
                                                            {tr_bl_tl_br,bl_tr_tl_br,br_tl_tr_bl},
                                                            {tr_bl_tl_br,bl_tr_tl_br,br_tl_bl_tr},
                                                            {tr_bl_tl_br,bl_tr_tl_br,bl_tr_tl_br},
                                                            {tr_bl_tl_br,bl_tr_tl_br,bl_tr_br_tl},
                                                            {tr_bl_tl_br,bl_tr_br_tl,tr_bl_br_tl},
                                                            {tr_bl_tl_br,bl_tr_br_tl,tr_bl_tl_br},
                                                            {tr_bl_tl_br,bl_tr_br_tl,tl_br_bl_tr},
                                                            {tr_bl_tl_br,bl_tr_br_tl,tl_br_tr_bl},
                                                            {tr_bl_tl_br,bl_tr_br_tl,br_tl_tr_bl},
                                                            {tr_bl_tl_br,bl_tr_br_tl,br_tl_bl_tr},
                                                            {tr_bl_tl_br,bl_tr_br_tl,bl_tr_tl_br},
                                                            {tr_bl_tl_br,bl_tr_br_tl,bl_tr_br_tl},
                                                            {tl_br_bl_tr,tr_bl_br_tl,tr_bl_br_tl},
                                                            {tl_br_bl_tr,tr_bl_br_tl,tr_bl_tl_br},
                                                            {tl_br_bl_tr,tr_bl_br_tl,tl_br_bl_tr},
                                                            {tl_br_bl_tr,tr_bl_br_tl,tl_br_tr_bl},
                                                            {tl_br_bl_tr,tr_bl_br_tl,br_tl_tr_bl},
                                                            {tl_br_bl_tr,tr_bl_br_tl,br_tl_bl_tr},
                                                            {tl_br_bl_tr,tr_bl_br_tl,bl_tr_tl_br},
                                                            {tl_br_bl_tr,tr_bl_br_tl,bl_tr_br_tl},
                                                            {tl_br_bl_tr,tr_bl_tl_br,tr_bl_br_tl},
                                                            {tl_br_bl_tr,tr_bl_tl_br,tr_bl_tl_br},
                                                            {tl_br_bl_tr,tr_bl_tl_br,tl_br_bl_tr},
                                                            {tl_br_bl_tr,tr_bl_tl_br,tl_br_tr_bl},
                                                            {tl_br_bl_tr,tr_bl_tl_br,br_tl_tr_bl},
                                                            {tl_br_bl_tr,tr_bl_tl_br,br_tl_bl_tr},
                                                            {tl_br_bl_tr,tr_bl_tl_br,bl_tr_tl_br},
                                                            {tl_br_bl_tr,tr_bl_tl_br,bl_tr_br_tl},
                                                            {tl_br_bl_tr,tl_br_bl_tr,tr_bl_br_tl},
                                                            {tl_br_bl_tr,tl_br_bl_tr,tr_bl_tl_br},
                                                            {tl_br_bl_tr,tl_br_bl_tr,tl_br_bl_tr},
                                                            {tl_br_bl_tr,tl_br_bl_tr,tl_br_tr_bl},
                                                            {tl_br_bl_tr,tl_br_bl_tr,br_tl_tr_bl},
                                                            {tl_br_bl_tr,tl_br_bl_tr,br_tl_bl_tr},
                                                            {tl_br_bl_tr,tl_br_bl_tr,bl_tr_tl_br},
                                                            {tl_br_bl_tr,tl_br_bl_tr,bl_tr_br_tl},
                                                            {tl_br_bl_tr,tl_br_tr_bl,tr_bl_br_tl},
                                                            {tl_br_bl_tr,tl_br_tr_bl,tr_bl_tl_br},
                                                            {tl_br_bl_tr,tl_br_tr_bl,tl_br_bl_tr},
                                                            {tl_br_bl_tr,tl_br_tr_bl,tl_br_tr_bl},
                                                            {tl_br_bl_tr,tl_br_tr_bl,br_tl_tr_bl},
                                                            {tl_br_bl_tr,tl_br_tr_bl,br_tl_bl_tr},
                                                            {tl_br_bl_tr,tl_br_tr_bl,bl_tr_tl_br},
                                                            {tl_br_bl_tr,tl_br_tr_bl,bl_tr_br_tl},
                                                            {tl_br_bl_tr,br_tl_tr_bl,tr_bl_br_tl},
                                                            {tl_br_bl_tr,br_tl_tr_bl,tr_bl_tl_br},
                                                            {tl_br_bl_tr,br_tl_tr_bl,tl_br_bl_tr},
                                                            {tl_br_bl_tr,br_tl_tr_bl,tl_br_tr_bl},
                                                            {tl_br_bl_tr,br_tl_tr_bl,br_tl_tr_bl},
                                                            {tl_br_bl_tr,br_tl_tr_bl,br_tl_bl_tr},
                                                            {tl_br_bl_tr,br_tl_tr_bl,bl_tr_tl_br},
                                                            {tl_br_bl_tr,br_tl_tr_bl,bl_tr_br_tl},
                                                            {tl_br_bl_tr,br_tl_bl_tr,tr_bl_br_tl},
                                                            {tl_br_bl_tr,br_tl_bl_tr,tr_bl_tl_br},
                                                            {tl_br_bl_tr,br_tl_bl_tr,tl_br_bl_tr},
                                                            {tl_br_bl_tr,br_tl_bl_tr,tl_br_tr_bl},
                                                            {tl_br_bl_tr,br_tl_bl_tr,br_tl_tr_bl},
                                                            {tl_br_bl_tr,br_tl_bl_tr,br_tl_bl_tr},
                                                            {tl_br_bl_tr,br_tl_bl_tr,bl_tr_tl_br},
                                                            {tl_br_bl_tr,br_tl_bl_tr,bl_tr_br_tl},
                                                            {tl_br_bl_tr,bl_tr_tl_br,tr_bl_br_tl},
                                                            {tl_br_bl_tr,bl_tr_tl_br,tr_bl_tl_br},
                                                            {tl_br_bl_tr,bl_tr_tl_br,tl_br_bl_tr},
                                                            {tl_br_bl_tr,bl_tr_tl_br,tl_br_tr_bl},
                                                            {tl_br_bl_tr,bl_tr_tl_br,br_tl_tr_bl},
                                                            {tl_br_bl_tr,bl_tr_tl_br,br_tl_bl_tr},
                                                            {tl_br_bl_tr,bl_tr_tl_br,bl_tr_tl_br},
                                                            {tl_br_bl_tr,bl_tr_tl_br,bl_tr_br_tl},
                                                            {tl_br_bl_tr,bl_tr_br_tl,tr_bl_br_tl},
                                                            {tl_br_bl_tr,bl_tr_br_tl,tr_bl_tl_br},
                                                            {tl_br_bl_tr,bl_tr_br_tl,tl_br_bl_tr},
                                                            {tl_br_bl_tr,bl_tr_br_tl,tl_br_tr_bl},
                                                            {tl_br_bl_tr,bl_tr_br_tl,br_tl_tr_bl},
                                                            {tl_br_bl_tr,bl_tr_br_tl,br_tl_bl_tr},
                                                            {tl_br_bl_tr,bl_tr_br_tl,bl_tr_tl_br},
                                                            {tl_br_bl_tr,bl_tr_br_tl,bl_tr_br_tl},
                                                            {tl_br_tr_bl,tr_bl_br_tl,tr_bl_br_tl},
                                                            {tl_br_tr_bl,tr_bl_br_tl,tr_bl_tl_br},
                                                            {tl_br_tr_bl,tr_bl_br_tl,tl_br_bl_tr},
                                                            {tl_br_tr_bl,tr_bl_br_tl,tl_br_tr_bl},
                                                            {tl_br_tr_bl,tr_bl_br_tl,br_tl_tr_bl},
                                                            {tl_br_tr_bl,tr_bl_br_tl,br_tl_bl_tr},
                                                            {tl_br_tr_bl,tr_bl_br_tl,bl_tr_tl_br},
                                                            {tl_br_tr_bl,tr_bl_br_tl,bl_tr_br_tl},
                                                            {tl_br_tr_bl,tr_bl_tl_br,tr_bl_br_tl},
                                                            {tl_br_tr_bl,tr_bl_tl_br,tr_bl_tl_br},
                                                            {tl_br_tr_bl,tr_bl_tl_br,tl_br_bl_tr},
                                                            {tl_br_tr_bl,tr_bl_tl_br,tl_br_tr_bl},
                                                            {tl_br_tr_bl,tr_bl_tl_br,br_tl_tr_bl},
                                                            {tl_br_tr_bl,tr_bl_tl_br,br_tl_bl_tr},
                                                            {tl_br_tr_bl,tr_bl_tl_br,bl_tr_tl_br},
                                                            {tl_br_tr_bl,tr_bl_tl_br,bl_tr_br_tl},
                                                            {tl_br_tr_bl,tl_br_bl_tr,tr_bl_br_tl},
                                                            {tl_br_tr_bl,tl_br_bl_tr,tr_bl_tl_br},
                                                            {tl_br_tr_bl,tl_br_bl_tr,tl_br_bl_tr},
                                                            {tl_br_tr_bl,tl_br_bl_tr,tl_br_tr_bl},
                                                            {tl_br_tr_bl,tl_br_bl_tr,br_tl_tr_bl},
                                                            {tl_br_tr_bl,tl_br_bl_tr,br_tl_bl_tr},
                                                            {tl_br_tr_bl,tl_br_bl_tr,bl_tr_tl_br},
                                                            {tl_br_tr_bl,tl_br_bl_tr,bl_tr_br_tl},
                                                            {tl_br_tr_bl,tl_br_tr_bl,tr_bl_br_tl},
                                                            {tl_br_tr_bl,tl_br_tr_bl,tr_bl_tl_br},
                                                            {tl_br_tr_bl,tl_br_tr_bl,tl_br_bl_tr},
                                                            {tl_br_tr_bl,tl_br_tr_bl,tl_br_tr_bl},
                                                            {tl_br_tr_bl,tl_br_tr_bl,br_tl_tr_bl},
                                                            {tl_br_tr_bl,tl_br_tr_bl,br_tl_bl_tr},
                                                            {tl_br_tr_bl,tl_br_tr_bl,bl_tr_tl_br},
                                                            {tl_br_tr_bl,tl_br_tr_bl,bl_tr_br_tl},
                                                            {tl_br_tr_bl,br_tl_tr_bl,tr_bl_br_tl},
                                                            {tl_br_tr_bl,br_tl_tr_bl,tr_bl_tl_br},
                                                            {tl_br_tr_bl,br_tl_tr_bl,tl_br_bl_tr},
                                                            {tl_br_tr_bl,br_tl_tr_bl,tl_br_tr_bl},
                                                            {tl_br_tr_bl,br_tl_tr_bl,br_tl_tr_bl},
                                                            {tl_br_tr_bl,br_tl_tr_bl,br_tl_bl_tr},
                                                            {tl_br_tr_bl,br_tl_tr_bl,bl_tr_tl_br},
                                                            {tl_br_tr_bl,br_tl_tr_bl,bl_tr_br_tl},
                                                            {tl_br_tr_bl,br_tl_bl_tr,tr_bl_br_tl},
                                                            {tl_br_tr_bl,br_tl_bl_tr,tr_bl_tl_br},
                                                            {tl_br_tr_bl,br_tl_bl_tr,tl_br_bl_tr},
                                                            {tl_br_tr_bl,br_tl_bl_tr,tl_br_tr_bl},
                                                            {tl_br_tr_bl,br_tl_bl_tr,br_tl_tr_bl},
                                                            {tl_br_tr_bl,br_tl_bl_tr,br_tl_bl_tr},
                                                            {tl_br_tr_bl,br_tl_bl_tr,bl_tr_tl_br},
                                                            {tl_br_tr_bl,br_tl_bl_tr,bl_tr_br_tl},
                                                            {tl_br_tr_bl,bl_tr_tl_br,tr_bl_br_tl},
                                                            {tl_br_tr_bl,bl_tr_tl_br,tr_bl_tl_br},
                                                            {tl_br_tr_bl,bl_tr_tl_br,tl_br_bl_tr},
                                                            {tl_br_tr_bl,bl_tr_tl_br,tl_br_tr_bl},
                                                            {tl_br_tr_bl,bl_tr_tl_br,br_tl_tr_bl},
                                                            {tl_br_tr_bl,bl_tr_tl_br,br_tl_bl_tr},
                                                            {tl_br_tr_bl,bl_tr_tl_br,bl_tr_tl_br},
                                                            {tl_br_tr_bl,bl_tr_tl_br,bl_tr_br_tl},
                                                            {tl_br_tr_bl,bl_tr_br_tl,tr_bl_br_tl},
                                                            {tl_br_tr_bl,bl_tr_br_tl,tr_bl_tl_br},
                                                            {tl_br_tr_bl,bl_tr_br_tl,tl_br_bl_tr},
                                                            {tl_br_tr_bl,bl_tr_br_tl,tl_br_tr_bl},
                                                            {tl_br_tr_bl,bl_tr_br_tl,br_tl_tr_bl},
                                                            {tl_br_tr_bl,bl_tr_br_tl,br_tl_bl_tr},
                                                            {tl_br_tr_bl,bl_tr_br_tl,bl_tr_tl_br},
                                                            {tl_br_tr_bl,bl_tr_br_tl,bl_tr_br_tl},
                                                            {br_tl_tr_bl,tr_bl_br_tl,tr_bl_br_tl},
                                                            {br_tl_tr_bl,tr_bl_br_tl,tr_bl_tl_br},
                                                            {br_tl_tr_bl,tr_bl_br_tl,tl_br_bl_tr},
                                                            {br_tl_tr_bl,tr_bl_br_tl,tl_br_tr_bl},
                                                            {br_tl_tr_bl,tr_bl_br_tl,br_tl_tr_bl},
                                                            {br_tl_tr_bl,tr_bl_br_tl,br_tl_bl_tr},
                                                            {br_tl_tr_bl,tr_bl_br_tl,bl_tr_tl_br},
                                                            {br_tl_tr_bl,tr_bl_br_tl,bl_tr_br_tl},
                                                            {br_tl_tr_bl,tr_bl_tl_br,tr_bl_br_tl},
                                                            {br_tl_tr_bl,tr_bl_tl_br,tr_bl_tl_br},
                                                            {br_tl_tr_bl,tr_bl_tl_br,tl_br_bl_tr},
                                                            {br_tl_tr_bl,tr_bl_tl_br,tl_br_tr_bl},
                                                            {br_tl_tr_bl,tr_bl_tl_br,br_tl_tr_bl},
                                                            {br_tl_tr_bl,tr_bl_tl_br,br_tl_bl_tr},
                                                            {br_tl_tr_bl,tr_bl_tl_br,bl_tr_tl_br},
                                                            {br_tl_tr_bl,tr_bl_tl_br,bl_tr_br_tl},
                                                            {br_tl_tr_bl,tl_br_bl_tr,tr_bl_br_tl},
                                                            {br_tl_tr_bl,tl_br_bl_tr,tr_bl_tl_br},
                                                            {br_tl_tr_bl,tl_br_bl_tr,tl_br_bl_tr},
                                                            {br_tl_tr_bl,tl_br_bl_tr,tl_br_tr_bl},
                                                            {br_tl_tr_bl,tl_br_bl_tr,br_tl_tr_bl},
                                                            {br_tl_tr_bl,tl_br_bl_tr,br_tl_bl_tr},
                                                            {br_tl_tr_bl,tl_br_bl_tr,bl_tr_tl_br},
                                                            {br_tl_tr_bl,tl_br_bl_tr,bl_tr_br_tl},
                                                            {br_tl_tr_bl,tl_br_tr_bl,tr_bl_br_tl},
                                                            {br_tl_tr_bl,tl_br_tr_bl,tr_bl_tl_br},
                                                            {br_tl_tr_bl,tl_br_tr_bl,tl_br_bl_tr},
                                                            {br_tl_tr_bl,tl_br_tr_bl,tl_br_tr_bl},
                                                            {br_tl_tr_bl,tl_br_tr_bl,br_tl_tr_bl},
                                                            {br_tl_tr_bl,tl_br_tr_bl,br_tl_bl_tr},
                                                            {br_tl_tr_bl,tl_br_tr_bl,bl_tr_tl_br},
                                                            {br_tl_tr_bl,tl_br_tr_bl,bl_tr_br_tl},
                                                            {br_tl_tr_bl,br_tl_tr_bl,tr_bl_br_tl},
                                                            {br_tl_tr_bl,br_tl_tr_bl,tr_bl_tl_br},
                                                            {br_tl_tr_bl,br_tl_tr_bl,tl_br_bl_tr},
                                                            {br_tl_tr_bl,br_tl_tr_bl,tl_br_tr_bl},
                                                            {br_tl_tr_bl,br_tl_tr_bl,br_tl_tr_bl},
                                                            {br_tl_tr_bl,br_tl_tr_bl,br_tl_bl_tr},
                                                            {br_tl_tr_bl,br_tl_tr_bl,bl_tr_tl_br},
                                                            {br_tl_tr_bl,br_tl_tr_bl,bl_tr_br_tl},
                                                            {br_tl_tr_bl,br_tl_bl_tr,tr_bl_br_tl},
                                                            {br_tl_tr_bl,br_tl_bl_tr,tr_bl_tl_br},
                                                            {br_tl_tr_bl,br_tl_bl_tr,tl_br_bl_tr},
                                                            {br_tl_tr_bl,br_tl_bl_tr,tl_br_tr_bl},
                                                            {br_tl_tr_bl,br_tl_bl_tr,br_tl_tr_bl},
                                                            {br_tl_tr_bl,br_tl_bl_tr,br_tl_bl_tr},
                                                            {br_tl_tr_bl,br_tl_bl_tr,bl_tr_tl_br},
                                                            {br_tl_tr_bl,br_tl_bl_tr,bl_tr_br_tl},
                                                            {br_tl_tr_bl,bl_tr_tl_br,tr_bl_br_tl},
                                                            {br_tl_tr_bl,bl_tr_tl_br,tr_bl_tl_br},
                                                            {br_tl_tr_bl,bl_tr_tl_br,tl_br_bl_tr},
                                                            {br_tl_tr_bl,bl_tr_tl_br,tl_br_tr_bl},
                                                            {br_tl_tr_bl,bl_tr_tl_br,br_tl_tr_bl},
                                                            {br_tl_tr_bl,bl_tr_tl_br,br_tl_bl_tr},
                                                            {br_tl_tr_bl,bl_tr_tl_br,bl_tr_tl_br},
                                                            {br_tl_tr_bl,bl_tr_tl_br,bl_tr_br_tl},
                                                            {br_tl_tr_bl,bl_tr_br_tl,tr_bl_br_tl},
                                                            {br_tl_tr_bl,bl_tr_br_tl,tr_bl_tl_br},
                                                            {br_tl_tr_bl,bl_tr_br_tl,tl_br_bl_tr},
                                                            {br_tl_tr_bl,bl_tr_br_tl,tl_br_tr_bl},
                                                            {br_tl_tr_bl,bl_tr_br_tl,br_tl_tr_bl},
                                                            {br_tl_tr_bl,bl_tr_br_tl,br_tl_bl_tr},
                                                            {br_tl_tr_bl,bl_tr_br_tl,bl_tr_tl_br},
                                                            {br_tl_tr_bl,bl_tr_br_tl,bl_tr_br_tl},
                                                            {br_tl_bl_tr,tr_bl_br_tl,tr_bl_br_tl},
                                                            {br_tl_bl_tr,tr_bl_br_tl,tr_bl_tl_br},
                                                            {br_tl_bl_tr,tr_bl_br_tl,tl_br_bl_tr},
                                                            {br_tl_bl_tr,tr_bl_br_tl,tl_br_tr_bl},
                                                            {br_tl_bl_tr,tr_bl_br_tl,br_tl_tr_bl},
                                                            {br_tl_bl_tr,tr_bl_br_tl,br_tl_bl_tr},
                                                            {br_tl_bl_tr,tr_bl_br_tl,bl_tr_tl_br},
                                                            {br_tl_bl_tr,tr_bl_br_tl,bl_tr_br_tl},
                                                            {br_tl_bl_tr,tr_bl_tl_br,tr_bl_br_tl},
                                                            {br_tl_bl_tr,tr_bl_tl_br,tr_bl_tl_br},
                                                            {br_tl_bl_tr,tr_bl_tl_br,tl_br_bl_tr},
                                                            {br_tl_bl_tr,tr_bl_tl_br,tl_br_tr_bl},
                                                            {br_tl_bl_tr,tr_bl_tl_br,br_tl_tr_bl},
                                                            {br_tl_bl_tr,tr_bl_tl_br,br_tl_bl_tr},
                                                            {br_tl_bl_tr,tr_bl_tl_br,bl_tr_tl_br},
                                                            {br_tl_bl_tr,tr_bl_tl_br,bl_tr_br_tl},
                                                            {br_tl_bl_tr,tl_br_bl_tr,tr_bl_br_tl},
                                                            {br_tl_bl_tr,tl_br_bl_tr,tr_bl_tl_br},
                                                            {br_tl_bl_tr,tl_br_bl_tr,tl_br_bl_tr},
                                                            {br_tl_bl_tr,tl_br_bl_tr,tl_br_tr_bl},
                                                            {br_tl_bl_tr,tl_br_bl_tr,br_tl_tr_bl},
                                                            {br_tl_bl_tr,tl_br_bl_tr,br_tl_bl_tr},
                                                            {br_tl_bl_tr,tl_br_bl_tr,bl_tr_tl_br},
                                                            {br_tl_bl_tr,tl_br_bl_tr,bl_tr_br_tl},
                                                            {br_tl_bl_tr,tl_br_tr_bl,tr_bl_br_tl},
                                                            {br_tl_bl_tr,tl_br_tr_bl,tr_bl_tl_br},
                                                            {br_tl_bl_tr,tl_br_tr_bl,tl_br_bl_tr},
                                                            {br_tl_bl_tr,tl_br_tr_bl,tl_br_tr_bl},
                                                            {br_tl_bl_tr,tl_br_tr_bl,br_tl_tr_bl},
                                                            {br_tl_bl_tr,tl_br_tr_bl,br_tl_bl_tr},
                                                            {br_tl_bl_tr,tl_br_tr_bl,bl_tr_tl_br},
                                                            {br_tl_bl_tr,tl_br_tr_bl,bl_tr_br_tl},
                                                            {br_tl_bl_tr,br_tl_tr_bl,tr_bl_br_tl},
                                                            {br_tl_bl_tr,br_tl_tr_bl,tr_bl_tl_br},
                                                            {br_tl_bl_tr,br_tl_tr_bl,tl_br_bl_tr},
                                                            {br_tl_bl_tr,br_tl_tr_bl,tl_br_tr_bl},
                                                            {br_tl_bl_tr,br_tl_tr_bl,br_tl_tr_bl},
                                                            {br_tl_bl_tr,br_tl_tr_bl,br_tl_bl_tr},
                                                            {br_tl_bl_tr,br_tl_tr_bl,bl_tr_tl_br},
                                                            {br_tl_bl_tr,br_tl_tr_bl,bl_tr_br_tl},
                                                            {br_tl_bl_tr,br_tl_bl_tr,tr_bl_br_tl},
                                                            {br_tl_bl_tr,br_tl_bl_tr,tr_bl_tl_br},
                                                            {br_tl_bl_tr,br_tl_bl_tr,tl_br_bl_tr},
                                                            {br_tl_bl_tr,br_tl_bl_tr,tl_br_tr_bl},
                                                            {br_tl_bl_tr,br_tl_bl_tr,br_tl_tr_bl},
                                                            {br_tl_bl_tr,br_tl_bl_tr,br_tl_bl_tr},
                                                            {br_tl_bl_tr,br_tl_bl_tr,bl_tr_tl_br},
                                                            {br_tl_bl_tr,br_tl_bl_tr,bl_tr_br_tl},
                                                            {br_tl_bl_tr,bl_tr_tl_br,tr_bl_br_tl},
                                                            {br_tl_bl_tr,bl_tr_tl_br,tr_bl_tl_br},
                                                            {br_tl_bl_tr,bl_tr_tl_br,tl_br_bl_tr},
                                                            {br_tl_bl_tr,bl_tr_tl_br,tl_br_tr_bl},
                                                            {br_tl_bl_tr,bl_tr_tl_br,br_tl_tr_bl},
                                                            {br_tl_bl_tr,bl_tr_tl_br,br_tl_bl_tr},
                                                            {br_tl_bl_tr,bl_tr_tl_br,bl_tr_tl_br},
                                                            {br_tl_bl_tr,bl_tr_tl_br,bl_tr_br_tl},
                                                            {br_tl_bl_tr,bl_tr_br_tl,tr_bl_br_tl},
                                                            {br_tl_bl_tr,bl_tr_br_tl,tr_bl_tl_br},
                                                            {br_tl_bl_tr,bl_tr_br_tl,tl_br_bl_tr},
                                                            {br_tl_bl_tr,bl_tr_br_tl,tl_br_tr_bl},
                                                            {br_tl_bl_tr,bl_tr_br_tl,br_tl_tr_bl},
                                                            {br_tl_bl_tr,bl_tr_br_tl,br_tl_bl_tr},
                                                            {br_tl_bl_tr,bl_tr_br_tl,bl_tr_tl_br},
                                                            {br_tl_bl_tr,bl_tr_br_tl,bl_tr_br_tl},
                                                            {bl_tr_tl_br,tr_bl_br_tl,tr_bl_br_tl},
                                                            {bl_tr_tl_br,tr_bl_br_tl,tr_bl_tl_br},
                                                            {bl_tr_tl_br,tr_bl_br_tl,tl_br_bl_tr},
                                                            {bl_tr_tl_br,tr_bl_br_tl,tl_br_tr_bl},
                                                            {bl_tr_tl_br,tr_bl_br_tl,br_tl_tr_bl},
                                                            {bl_tr_tl_br,tr_bl_br_tl,br_tl_bl_tr},
                                                            {bl_tr_tl_br,tr_bl_br_tl,bl_tr_tl_br},
                                                            {bl_tr_tl_br,tr_bl_br_tl,bl_tr_br_tl},
                                                            {bl_tr_tl_br,tr_bl_tl_br,tr_bl_br_tl},
                                                            {bl_tr_tl_br,tr_bl_tl_br,tr_bl_tl_br},
                                                            {bl_tr_tl_br,tr_bl_tl_br,tl_br_bl_tr},
                                                            {bl_tr_tl_br,tr_bl_tl_br,tl_br_tr_bl},
                                                            {bl_tr_tl_br,tr_bl_tl_br,br_tl_tr_bl},
                                                            {bl_tr_tl_br,tr_bl_tl_br,br_tl_bl_tr},
                                                            {bl_tr_tl_br,tr_bl_tl_br,bl_tr_tl_br},
                                                            {bl_tr_tl_br,tr_bl_tl_br,bl_tr_br_tl},
                                                            {bl_tr_tl_br,tl_br_bl_tr,tr_bl_br_tl},
                                                            {bl_tr_tl_br,tl_br_bl_tr,tr_bl_tl_br},
                                                            {bl_tr_tl_br,tl_br_bl_tr,tl_br_bl_tr},
                                                            {bl_tr_tl_br,tl_br_bl_tr,tl_br_tr_bl},
                                                            {bl_tr_tl_br,tl_br_bl_tr,br_tl_tr_bl},
                                                            {bl_tr_tl_br,tl_br_bl_tr,br_tl_bl_tr},
                                                            {bl_tr_tl_br,tl_br_bl_tr,bl_tr_tl_br},
                                                            {bl_tr_tl_br,tl_br_bl_tr,bl_tr_br_tl},
                                                            {bl_tr_tl_br,tl_br_tr_bl,tr_bl_br_tl},
                                                            {bl_tr_tl_br,tl_br_tr_bl,tr_bl_tl_br},
                                                            {bl_tr_tl_br,tl_br_tr_bl,tl_br_bl_tr},
                                                            {bl_tr_tl_br,tl_br_tr_bl,tl_br_tr_bl},
                                                            {bl_tr_tl_br,tl_br_tr_bl,br_tl_tr_bl},
                                                            {bl_tr_tl_br,tl_br_tr_bl,br_tl_bl_tr},
                                                            {bl_tr_tl_br,tl_br_tr_bl,bl_tr_tl_br},
                                                            {bl_tr_tl_br,tl_br_tr_bl,bl_tr_br_tl},
                                                            {bl_tr_tl_br,br_tl_tr_bl,tr_bl_br_tl},
                                                            {bl_tr_tl_br,br_tl_tr_bl,tr_bl_tl_br},
                                                            {bl_tr_tl_br,br_tl_tr_bl,tl_br_bl_tr},
                                                            {bl_tr_tl_br,br_tl_tr_bl,tl_br_tr_bl},
                                                            {bl_tr_tl_br,br_tl_tr_bl,br_tl_tr_bl},
                                                            {bl_tr_tl_br,br_tl_tr_bl,br_tl_bl_tr},
                                                            {bl_tr_tl_br,br_tl_tr_bl,bl_tr_tl_br},
                                                            {bl_tr_tl_br,br_tl_tr_bl,bl_tr_br_tl},
                                                            {bl_tr_tl_br,br_tl_bl_tr,tr_bl_br_tl},
                                                            {bl_tr_tl_br,br_tl_bl_tr,tr_bl_tl_br},
                                                            {bl_tr_tl_br,br_tl_bl_tr,tl_br_bl_tr},
                                                            {bl_tr_tl_br,br_tl_bl_tr,tl_br_tr_bl},
                                                            {bl_tr_tl_br,br_tl_bl_tr,br_tl_tr_bl},
                                                            {bl_tr_tl_br,br_tl_bl_tr,br_tl_bl_tr},
                                                            {bl_tr_tl_br,br_tl_bl_tr,bl_tr_tl_br},
                                                            {bl_tr_tl_br,br_tl_bl_tr,bl_tr_br_tl},
                                                            {bl_tr_tl_br,bl_tr_tl_br,tr_bl_br_tl},
                                                            {bl_tr_tl_br,bl_tr_tl_br,tr_bl_tl_br},
                                                            {bl_tr_tl_br,bl_tr_tl_br,tl_br_bl_tr},
                                                            {bl_tr_tl_br,bl_tr_tl_br,tl_br_tr_bl},
                                                            {bl_tr_tl_br,bl_tr_tl_br,br_tl_tr_bl},
                                                            {bl_tr_tl_br,bl_tr_tl_br,br_tl_bl_tr},
                                                            {bl_tr_tl_br,bl_tr_tl_br,bl_tr_tl_br},
                                                            {bl_tr_tl_br,bl_tr_tl_br,bl_tr_br_tl},
                                                            {bl_tr_tl_br,bl_tr_br_tl,tr_bl_br_tl},
                                                            {bl_tr_tl_br,bl_tr_br_tl,tr_bl_tl_br},
                                                            {bl_tr_tl_br,bl_tr_br_tl,tl_br_bl_tr},
                                                            {bl_tr_tl_br,bl_tr_br_tl,tl_br_tr_bl},
                                                            {bl_tr_tl_br,bl_tr_br_tl,br_tl_tr_bl},
                                                            {bl_tr_tl_br,bl_tr_br_tl,br_tl_bl_tr},
                                                            {bl_tr_tl_br,bl_tr_br_tl,bl_tr_tl_br},
                                                            {bl_tr_tl_br,bl_tr_br_tl,bl_tr_br_tl},
                                                            {bl_tr_br_tl,tr_bl_br_tl,tr_bl_br_tl},
                                                            {bl_tr_br_tl,tr_bl_br_tl,tr_bl_tl_br},
                                                            {bl_tr_br_tl,tr_bl_br_tl,tl_br_bl_tr},
                                                            {bl_tr_br_tl,tr_bl_br_tl,tl_br_tr_bl},
                                                            {bl_tr_br_tl,tr_bl_br_tl,br_tl_tr_bl},
                                                            {bl_tr_br_tl,tr_bl_br_tl,br_tl_bl_tr},
                                                            {bl_tr_br_tl,tr_bl_br_tl,bl_tr_tl_br},
                                                            {bl_tr_br_tl,tr_bl_br_tl,bl_tr_br_tl},
                                                            {bl_tr_br_tl,tr_bl_tl_br,tr_bl_br_tl},
                                                            {bl_tr_br_tl,tr_bl_tl_br,tr_bl_tl_br},
                                                            {bl_tr_br_tl,tr_bl_tl_br,tl_br_bl_tr},
                                                            {bl_tr_br_tl,tr_bl_tl_br,tl_br_tr_bl},
                                                            {bl_tr_br_tl,tr_bl_tl_br,br_tl_tr_bl},
                                                            {bl_tr_br_tl,tr_bl_tl_br,br_tl_bl_tr},
                                                            {bl_tr_br_tl,tr_bl_tl_br,bl_tr_tl_br},
                                                            {bl_tr_br_tl,tr_bl_tl_br,bl_tr_br_tl},
                                                            {bl_tr_br_tl,tl_br_bl_tr,tr_bl_br_tl},
                                                            {bl_tr_br_tl,tl_br_bl_tr,tr_bl_tl_br},
                                                            {bl_tr_br_tl,tl_br_bl_tr,tl_br_bl_tr},
                                                            {bl_tr_br_tl,tl_br_bl_tr,tl_br_tr_bl},
                                                            {bl_tr_br_tl,tl_br_bl_tr,br_tl_tr_bl},
                                                            {bl_tr_br_tl,tl_br_bl_tr,br_tl_bl_tr},
                                                            {bl_tr_br_tl,tl_br_bl_tr,bl_tr_tl_br},
                                                            {bl_tr_br_tl,tl_br_bl_tr,bl_tr_br_tl},
                                                            {bl_tr_br_tl,tl_br_tr_bl,tr_bl_br_tl},
                                                            {bl_tr_br_tl,tl_br_tr_bl,tr_bl_tl_br},
                                                            {bl_tr_br_tl,tl_br_tr_bl,tl_br_bl_tr},
                                                            {bl_tr_br_tl,tl_br_tr_bl,tl_br_tr_bl},
                                                            {bl_tr_br_tl,tl_br_tr_bl,br_tl_tr_bl},
                                                            {bl_tr_br_tl,tl_br_tr_bl,br_tl_bl_tr},
                                                            {bl_tr_br_tl,tl_br_tr_bl,bl_tr_tl_br},
                                                            {bl_tr_br_tl,tl_br_tr_bl,bl_tr_br_tl},
                                                            {bl_tr_br_tl,br_tl_tr_bl,tr_bl_br_tl},
                                                            {bl_tr_br_tl,br_tl_tr_bl,tr_bl_tl_br},
                                                            {bl_tr_br_tl,br_tl_tr_bl,tl_br_bl_tr},
                                                            {bl_tr_br_tl,br_tl_tr_bl,tl_br_tr_bl},
                                                            {bl_tr_br_tl,br_tl_tr_bl,br_tl_tr_bl},
                                                            {bl_tr_br_tl,br_tl_tr_bl,br_tl_bl_tr},
                                                            {bl_tr_br_tl,br_tl_tr_bl,bl_tr_tl_br},
                                                            {bl_tr_br_tl,br_tl_tr_bl,bl_tr_br_tl},
                                                            {bl_tr_br_tl,br_tl_bl_tr,tr_bl_br_tl},
                                                            {bl_tr_br_tl,br_tl_bl_tr,tr_bl_tl_br},
                                                            {bl_tr_br_tl,br_tl_bl_tr,tl_br_bl_tr},
                                                            {bl_tr_br_tl,br_tl_bl_tr,tl_br_tr_bl},
                                                            {bl_tr_br_tl,br_tl_bl_tr,br_tl_tr_bl},
                                                            {bl_tr_br_tl,br_tl_bl_tr,br_tl_bl_tr},
                                                            {bl_tr_br_tl,br_tl_bl_tr,bl_tr_tl_br},
                                                            {bl_tr_br_tl,br_tl_bl_tr,bl_tr_br_tl},
                                                            {bl_tr_br_tl,bl_tr_tl_br,tr_bl_br_tl},
                                                            {bl_tr_br_tl,bl_tr_tl_br,tr_bl_tl_br},
                                                            {bl_tr_br_tl,bl_tr_tl_br,tl_br_bl_tr},
                                                            {bl_tr_br_tl,bl_tr_tl_br,tl_br_tr_bl},
                                                            {bl_tr_br_tl,bl_tr_tl_br,br_tl_tr_bl},
                                                            {bl_tr_br_tl,bl_tr_tl_br,br_tl_bl_tr},
                                                            {bl_tr_br_tl,bl_tr_tl_br,bl_tr_tl_br},
                                                            {bl_tr_br_tl,bl_tr_tl_br,bl_tr_br_tl},
                                                            {bl_tr_br_tl,bl_tr_br_tl,tr_bl_br_tl},
                                                            {bl_tr_br_tl,bl_tr_br_tl,tr_bl_tl_br},
                                                            {bl_tr_br_tl,bl_tr_br_tl,tl_br_bl_tr},
                                                            {bl_tr_br_tl,bl_tr_br_tl,tl_br_tr_bl},
                                                            {bl_tr_br_tl,bl_tr_br_tl,br_tl_tr_bl},
                                                            {bl_tr_br_tl,bl_tr_br_tl,br_tl_bl_tr},
                                                            {bl_tr_br_tl,bl_tr_br_tl,bl_tr_tl_br},
                                                            {bl_tr_br_tl,bl_tr_br_tl,bl_tr_br_tl}

                                                           };


    int fp_x = -1;
    int fp_y = -1;

    int n_points = path.size();

    for(int i = 0; i < n_points; i++){

        vector<Point> sub_points;

        if (i + depth <= n_points){
            for(int k = 0; k < depth; k++)
                sub_points.push_back(path[i+k]);
        } else {
            for(int k = 0; k < (n_points - i); k++)
                sub_points.push_back(path[i+k]);
        }


        int n_sub_points = sub_points.size();
        double min_d = LOCAL_INFINITY;
        vector<int> min_stitch_x(4, 0);
        vector<int> min_stitch_y(4, 0);

        for(auto product : stitich_functions){

            vector<Point> pb_points(n_sub_points);
            vector<Point> pf_points(n_sub_points);

            vector<int> first_stitch_x(4, 0);
            vector<int> first_stitch_y(4, 0);

            vector<int> stitch_x(4, 0);
            vector<int> stitch_y(4, 0);

            for(int j = 0; j < n_sub_points; j++){
                product[j](stitch_x, stitch_y, sub_points[j].get_x(), sub_points[j].get_y());

                if (j == 0){
                    first_stitch_x = stitch_x;
                    first_stitch_y = stitch_y;
                }

                pb_points[j] = Point(stitch_x[0], stitch_y[0]);
                pf_points[j] = Point(stitch_x[3], stitch_y[3]);
            }

            // Check if the end point of the previous stitch
            // does not match the starting point of the first stitch.
            if ((fp_x == first_stitch_x[0]) && (fp_y == first_stitch_y[0]))
                continue;

            bool is_order_valid = check_if_stitch_order_is_valid(pb_points, pf_points);
            if (is_order_valid == false)
                continue;

            // Take in to account the previous point in the path.
            // For i == 0 there are no previous points.
            if (i != 0){
                // The first pb point is not used in the calcualtion of stitch_d
                // so we use a default constructor.
                // deque would be more appropriate (nned to change in the future)
                pb_points.insert(pb_points.begin(), Point());
                pf_points.insert(pf_points.begin(), Point(fp_x, fp_y));
            }


            double stitch_d = calculate_stitch_d_distance(pb_points, pf_points);

            if (stitch_d < min_d){
                min_d = stitch_d;
                min_stitch_x = first_stitch_x;
                min_stitch_y = first_stitch_y;
            }


        }

        ret.push_back(to_string(min_stitch_x[0]) + " " + to_string(min_stitch_y[0]));
        ret.push_back(to_string(min_stitch_x[1]) + " " + to_string(min_stitch_y[1]));
        ret.push_back(to_string(min_stitch_x[2]) + " " + to_string(min_stitch_y[2]));
        ret.push_back(to_string(min_stitch_x[3]) + " " + to_string(min_stitch_y[3]));
        
        fp_x = min_stitch_x[3];
        fp_y = min_stitch_y[3];

    }

}


vector<string> PointCollection::make_stitchs(int depth){

    vector<string> ret;

    for(int letter_index = 0; letter_index < NUMBER_OF_LETTERS; letter_index++){
        int n_points = m_point_collection[letter_index].size();
        if (n_points == 0)
            continue;

        char letter = revert_get_letter_index(letter_index);
        ret.push_back( string(1, letter) );
        
        product_stitch(ret, m_point_collection[letter_index], depth);
    }

    return ret;
}


void PointCollection::rotate_collection(){

    for(int letter_index = 0; letter_index < NUMBER_OF_LETTERS; letter_index++){

        int n_points = m_point_collection[letter_index].size();

        if (n_points == 0)
            continue;

        int shift_index = 0;
        int minimal_shift = 0;
        double min_d = get_path_length_in_collection(m_point_collection[letter_index], letter_index);
        double new_d = min_d;
        for(int i = 0; i < n_points; i++){

            // distance first - second
            int b1 = m_point_collection[letter_index][0].get_b();
            int b2 = m_point_collection[letter_index][1].get_b();
            double dfs = m_distance_matrix[letter_index][b1][b2];

            std::rotate(m_point_collection[letter_index].begin(), m_point_collection[letter_index].begin() + 1, m_point_collection[letter_index].end());
            shift_index = shift_index + 1;

            // distance last - last but one
            b1 = m_point_collection[letter_index][n_points - 1].get_b();
            b2 = m_point_collection[letter_index][n_points - 2].get_b();
            double dll = m_distance_matrix[letter_index][b1][b2];

            new_d = new_d - dfs + dll;

            if (new_d < min_d){
                minimal_shift = shift_index;
                min_d = new_d;
            }
        }

        std::rotate(m_point_collection[letter_index].begin(), m_point_collection[letter_index].begin() + minimal_shift, m_point_collection[letter_index].end());
    }
}


void x_y_to_corners(vector<int> &corners, int x, int y){

    int tr_x = x + 1;
    int tr_y = y;

    int bl_x = x;
    int bl_y = y + 1;

    int br_x = x + 1;
    int br_y = y + 1;

    int tl_x = x;
    int tl_y = y;

    corners[0] = tl_x;
    corners[1] = tl_y;
    corners[2] = tr_x;
    corners[3] = tr_y;
    corners[4] = bl_x;
    corners[5] = bl_y;
    corners[6] = br_x;
    corners[7] = br_y;

}


void tr_bl_br_tl(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y){

    stitch_x[0] = x + 1; // tr_x = x + 1
    stitch_x[1] = x;     // bl_x = x
    stitch_x[2] = x + 1; // br_x = x + 1
    stitch_x[3] = x;     // tl_x = x
    
    stitch_y[0] = y;     // tr_y = y
    stitch_y[1] = y + 1; // bl_y = y + 1
    stitch_y[2] = y + 1; // br_y = y + 1
    stitch_y[3] = y;     // tl_y = y
}


void tr_bl_tl_br(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y){

    stitch_x[0] = x + 1; // tr_x = x + 1
    stitch_x[1] = x;     // bl_x = x
    stitch_x[2] = x;     // tl_x = x
    stitch_x[3] = x + 1; // br_x = x + 1
    
    stitch_y[0] = y;     // tr_y = y
    stitch_y[1] = y + 1; // bl_y = y + 1
    stitch_y[2] = y;     // tl_y = y
    stitch_y[3] = y + 1; // br_y = y + 1

}

void tl_br_bl_tr(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y){

    stitch_x[0] = x;     // tl_x = x
    stitch_x[1] = x + 1; // br_x = x + 1
    stitch_x[2] = x;     // bl_x = x
    stitch_x[3] = x + 1; // tr_x = x + 1
    
    stitch_y[0] = y;     // tl_y = y
    stitch_y[1] = y + 1; // br_y = y + 1
    stitch_y[2] = y + 1; // bl_y = y + 1
    stitch_y[3] = y;     // tr_y = y

}

void tl_br_tr_bl(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y){

    stitch_x[0] = x;     // tl_x = x
    stitch_x[1] = x + 1; // br_x = x + 1
    stitch_x[2] = x + 1; // tr_x = x + 1
    stitch_x[3] = x;     // bl_x = x
  
    stitch_y[0] = y;     // tl_y = y
    stitch_y[1] = y + 1; // br_y = y + 1
    stitch_y[2] = y;     // tr_y = y
    stitch_y[3] = y + 1; // bl_y = y + 1

}


void br_tl_tr_bl(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y){

    stitch_x[0] = x + 1; // br_x = x + 1
    stitch_x[1] = x;     // tl_x = x
    stitch_x[2] = x + 1; // tr_x = x + 1
    stitch_x[3] = x;     // bl_x = x
  
    stitch_y[0] = y + 1; // br_y = y + 1
    stitch_y[1] = y;     // tl_y = y
    stitch_y[2] = y;     // tr_y = y
    stitch_y[3] = y + 1; // bl_y = y + 1

}

void br_tl_bl_tr(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y){

    stitch_x[0] = x + 1; // br_x = x + 1
    stitch_x[1] = x;     // tl_x = x
    stitch_x[2] = x;     // bl_x = x
    stitch_x[3] = x + 1; // tr_x = x + 1
  
    stitch_y[0] = y + 1; // br_y = y + 1
    stitch_y[1] = y;     // tl_y = y
    stitch_y[2] = y + 1; // bl_y = y + 1
    stitch_y[3] = y;     // tr_y = y

}


void bl_tr_tl_br(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y){

    stitch_x[0] = x;     // bl_x = x
    stitch_x[1] = x + 1; // tr_x = x + 1
    stitch_x[2] = x;     // tl_x = x
    stitch_x[3] = x + 1; // br_x = x + 1
  
    stitch_y[0] = y + 1; // bl_y = y + 1
    stitch_y[1] = y;     // tr_y = y
    stitch_y[2] = y;     // tl_y = y
    stitch_y[3] = y + 1; // br_y = y + 1

}

void bl_tr_br_tl(vector<int> &stitch_x, vector<int> &stitch_y, int x, int y){

    stitch_x[0] = x;     // bl_x = x
    stitch_x[1] = x + 1; // tr_x = x + 1
    stitch_x[2] = x + 1; // br_x = x + 1
    stitch_x[3] = x;     // tl_x = x
  
    stitch_y[0] = y + 1; // bl_y = y + 1
    stitch_y[1] = y;     // tr_y = y
    stitch_y[2] = y + 1; // br_y = y + 1
    stitch_y[3] = y;     // tl_y = y

}


class CrossStitch {
public:
    vector<string> embroider(vector<string> pattern) {

        function_pointer f_ptr;

        f_ptr = &bl_tr_br_tl;


        Point p1;
        Point p2(-1, -2, -2);

        //cerr << "P1 == p2: " << (p1 == p2) << endl;

        std::chrono::time_point<std::chrono::system_clock> start, end;


        start = std::chrono::system_clock::now();


        PointCollection point_collection(pattern);
        // point_collection.eprint_number_of_points_in_each_path();
        point_collection.eprint_average_path_length_of_collection();


        point_collection.nearest_neighbour_all_path_optimize();
        point_collection.eprint_average_path_length_of_collection();


        point_collection.two_opt_all_path_optimize();
        point_collection.eprint_average_path_length_of_collection();


        //point_collection.markov_monte_carlo_like_all_path_optimize();
        //point_collection.eprint_average_path_length_of_collection();

        point_collection.rotate_collection();
        point_collection.eprint_average_path_length_of_collection();


        double s = sqrt(2.0);
        vector<vector<vector<double>>> distance_matrix = {{
                                                   {LOCAL_INFINITY,   2, s, 2*s, 2},
                                                   {2,   LOCAL_INFINITY, s,   2, 2*s},
                                                   {s,   s, LOCAL_INFINITY,   s, s},
                                                   {2*s, 2, s,   LOCAL_INFINITY, 2},
                                                   {2, 2*s, s,   2, LOCAL_INFINITY},
                                                  }};

   
    /*
        vector<vector<vector<double>>> distance_matrix = {{
                                                   {0, 2, 0, 6, 0},
                                                   {2, 0, 3, 8, 5},
                                                   {0, 3, 0, 0, 7},
                                                   {6, 8, 0, 0, 9},
                                                   {0, 5, 7, 9, 0},
                                                  }};

   */

        /*
        vector<vector<vector<double>>> distance_matrix = {{
                                                   {0,  4, 0, 0, 0, 0, 0, 7, 0},//0
                                                   {4,  0, 8, 0, 0, 0, 0,11, 0},//1
                                                   {0,  8, 0, 7, 0, 4, 0, 0, 2},//2
                                                   {0,  0, 7, 0, 9,14, 0, 0, 0},//3
                                                   {0,  0, 0, 9, 0,10, 0, 0, 0},//4
                                                   {0,  0, 4,14,10, 0, 2, 0, 0},//5
                                                   {0,  0, 0, 0, 0, 2, 0, 1, 6},//6
                                                   {7, 11, 0, 0, 0, 0, 1, 0, 7},//7
                                                   {0,  0, 2, 0, 0, 0, 6, 7, 0},//8
                                                   }};

        */

        int li = 0;
    
        vector<Point> path = {Point(), Point(), Point(), Point(), Point()};
    
        point_collection.minimal_spannig_tree(path, li, distance_matrix);

        //point_collection.eprint_collection();
        //point_collection.eprint_distance_matrixes(1);
        vector<string> ret;
        bool flag = true;

        if (flag == true) {
            ret = point_collection.make_stitchs(2);

        } else {       
            int S = pattern.size();
            // for each color, for each cell (r, c) do two stitches (r+1, c)-(r, c+1)-(r+1, c+1)-(r, c)
            for (char col = 'a'; col <= 'z'; ++col) {
                bool first = true;
                for (int r = 0; r < S; ++r)
                for (int c = 0; c < S; ++c)
                    if (pattern[r][c] == col) {
                        if (first) {
                            ret.push_back(string(1, col));
                            first = false;
                        }
                        ret.push_back(to_string(r+1) + " " + to_string(c));
                        ret.push_back(to_string(r) + " " + to_string(c+1));
                        ret.push_back(to_string(r+1) + " " + to_string(c+1));
                        ret.push_back(to_string(r) + " " + to_string(c));
                    }
            }
        }


        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_in_main;
        elapsed_seconds_in_main = end-start;
        cerr << "elapsed time: " << elapsed_seconds_in_main.count() << endl;

        return ret;
    }
};


// -------8<------- end of solution submitted to the website -------8<-------

template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}

int main() {
    int S;
    cin >> S;
    vector<string> pattern(S);
    getVector(pattern);
    
    CrossStitch cs;
    vector<string> ret = cs.embroider(pattern);
    cout << ret.size() << endl;
    for (int i = 0; i < (int)ret.size(); ++i)
        cout << ret[i] << endl;
    cout.flush();
}

