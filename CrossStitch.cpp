#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>
#include <chrono>
#include <random>

#define NUMBER_OF_LETTERS 26
#define LOCAL_INFINITY 2147483647.0 // 2^31 -1

using namespace std;

inline int get_letter_index(char &letter) {
    return letter - 97;
}

inline char revert_get_letter_index(int &letter_index){
    return letter_index + 97;
}


class Point {
    private:
        int m_x;
        int m_y;
        int m_b; // # badge

    public:
        Point() { m_x = -1.0; m_y = -1.0; m_b = -1; }
        Point(int x, int y): m_x(x), m_y(y) { m_b = -1; }        
        Point(int x, int y, int b): m_x(x), m_y(y), m_b(b) {}

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


//-------------------------------------------------------
//--------------------Point------------------------------
//-------------------------------------------------------


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
        void eprint_distance_matrixes(int letter_index);
        void eprint_average_path_length_of_collection();

        double get_path_length_in_collection(vector<Point> &path);

        vector<Point> nearest_neighbour_single_path_optimize(vector<Point> &path, int &letter_index);
        void nearest_neighbour_all_path_minimize();

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

        apl = apl + get_path_length_in_collection(m_point_collection[letter_index]);
        nvp = nvp + 1;
    }
    cerr << "\nAverage path length: " << apl/nvp << endl << endl;
}


double PointCollection::get_path_length_in_collection(vector<Point> &path){

    int n_points = path.size();

    double pl = 0.0; // path length
    for(int i = 1; i < n_points; i++){
        double d = Point::euclidean_distance(path[i], path[i - 1]);
        pl = pl + d;
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


void PointCollection::nearest_neighbour_all_path_minimize(){

    for(int letter_index = 0; letter_index < NUMBER_OF_LETTERS; letter_index++){
        if (m_point_collection[letter_index].size() == 0)
            continue;

        m_point_collection[letter_index] = nearest_neighbour_single_path_optimize(m_point_collection[letter_index], letter_index);
    }

}


class CrossStitch {
public:
    vector<string> embroider(vector<string> pattern) {

        Point p1;
        Point p2(-1, -2, -2);

        //cerr << "P1 == p2: " << (p1 == p2) << endl;

        std::chrono::time_point<std::chrono::system_clock> start, end;


        start = std::chrono::system_clock::now();


        PointCollection point_collection(pattern);
        point_collection.eprint_average_path_length_of_collection();


        point_collection.nearest_neighbour_all_path_minimize();
        point_collection.eprint_average_path_length_of_collection();


        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_in_main;
        elapsed_seconds_in_main = end-start;
        cerr << "elapsed time: " << elapsed_seconds_in_main.count() << endl;

        //point_collection.eprint_collection();
        //point_collection.eprint_distance_matrixes(1);

        vector<string> ret;
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

