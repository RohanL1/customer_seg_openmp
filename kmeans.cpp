#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cstring>
#include </opt/homebrew/Cellar/libomp/15.0.4/include/omp.h>

using namespace std;

/// Global variable
int nthreads = 0 ;
string output_file_name = "output.csv";

/// Structure to store features of a point
struct Point {

    double x, y;     /// coordinates
    int cluster;     /// cluster ID
    double minDist;  /// Minimun distance to the nearest cluster

    /// Setting the default values
    Point() : x(0.0), y(0.0), cluster(-1), minDist(__DBL_MAX__) {}
    Point(double x, double y) : x(x), y(y), cluster(-1), minDist(__DBL_MAX__) {}

    
    /// Compute the Distance using Euclidean Distance formula between 2 points
    double distance(Point p) 
    {
        return ((p.x - x) * (p.x - x) + (p.y - y) * (p.y - y));
    }
};

/// Reading the input ".csv" file and converting it into vector of points
vector<Point> readcsv(string input_file_name) 
{
    vector<Point> points;
    string line;
    ifstream file(input_file_name);
    while (getline(file, line)) 
    {
        stringstream lineStream(line);
        string bit;
        double x, y;
        getline(lineStream, bit, ',');
        x = stof(bit);
        getline(lineStream, bit, '\n');
        y = stof(bit);
        points.push_back(Point(x, y));
    }
    return points;
}


/// Perform k-means clustering serially
/// @param points - pointer to vector of points
/// @param iterations - number of k means iterations
/// @param ncluster - the number of initial centroids

void kMeansClusteringSerial(vector<Point>* points, int iterations, int ncluster) {
    int n = points->size();

    /// Randomly initialise centroids
    /// The index of the centroid represents the cluster ID
    vector<Point> centroids;
    srand(time(0));
    for (int i = 0; i < ncluster; ++i) {
        centroids.push_back(points->at(rand() % n));
    }

    for (int i = 0; i < iterations; ++i) 
    {
        /// For each centroid, compute distance between each point and centroid
        /// and update min Distance and Cluster ID, if changed
        for (vector<Point>::iterator c = begin(centroids); c != end(centroids); ++c) 
        {
            int clusterId = c - begin(centroids);

            for (vector<Point>::iterator it = points->begin(); it != points->end(); ++it) 
            {
                Point p = *it;
                double dist = c->distance(p);
                if (dist < p.minDist) 
                {
                    p.minDist = dist;
                    p.cluster = clusterId;
                }
                *it = p;
            }
        }

        /// Vectors to keep track of data needed to compute new centroids
        vector<int> nPoints;
        vector<double> sumX, sumY;
        for (int j = 0; j < ncluster; ++j) 
        {
            nPoints.push_back(0);
            sumX.push_back(0.0);
            sumY.push_back(0.0);
        }
        for (vector<Point>::iterator it = points->begin(); it != points->end(); ++it) 
        {
            int clusterId = it->cluster;
            nPoints[clusterId] += 1;
            sumX[clusterId] += it->x;
            sumY[clusterId] += it->y;

            it->minDist = __DBL_MAX__;  /// reset min distance 
        }

        /// Compute the new centroids
        for (vector<Point>::iterator c = begin(centroids); c != end(centroids); ++c) 
        {
            int clusterId = c - begin(centroids);
            c->x = sumX[clusterId] / nPoints[clusterId];
            c->y = sumY[clusterId] / nPoints[clusterId];
        }
    }

    /// Writing output to csv
    ofstream output_file;
    output_file.open(output_file_name + "_serial.csv");
    output_file << "x,y,cluster_id" << endl;

    for (vector<Point>::iterator it = points->begin(); it != points->end(); ++it) 
    {
        output_file << it->x << "," << it->y << "," << it->cluster << endl;
    }
    output_file.close();
}

/// Perform k-means clustering Parallely
/// @param points - pointer to vector of points
/// @param iterations - number of k means iterations
/// @param ncluster - the number of initial centroids

void kMeansClusteringParallel (vector<Point>* points, int iterations, int ncluster) 
{
    int n = points->size(); /// Size of the data

    /// Vectors to keep track of data needed to compute new centroids
    vector<int> nPoints;
    vector<double> sumX, sumY;

    /// Randomly initialise centroids
    /// The index of the centroid represents the cluster ID
    vector<Point> centroids;
    srand(time(0));
    for (int i = 0; i < ncluster; ++i) 
    {
        centroids.push_back(points->at(rand() % n));
    }

    for (int i = 0; i < iterations; ++i) 
    {
        /// For each centroid, compute distance between each point and centroid
        /// and update min Distance and Cluster ID, if changed
        for (unsigned int j=0; j < centroids.size(); j++) 
        {
            int clusterId = j;
                #pragma omp parallel for
                for(unsigned int i = 0; i < (*points).size(); i++)
                {
                    double dist = centroids[j].distance((*points)[i]);
                    if (dist < (*points)[i].minDist) 
                    {
                        (*points)[i].minDist = dist;
                        (*points)[i].cluster = clusterId;
                    }
                }
        }
    
        #pragma omp parallel
        {
            /// cout << "\n Thread Count : " << omp_get_num_threads();
            #pragma omp single 
            {
                for (int j = 0; j < ncluster; ++j) 
                {
                    nPoints.push_back(0);
                    sumX.push_back(0.0);
                    sumY.push_back(0.0);
                }
            }
            #pragma omp barrier

            #pragma omp for
            for (unsigned int j=0; j < (*points).size(); j++) 
            {
                int clusterId = (*points)[j].cluster;
                nPoints[clusterId] += 1;
                sumX[clusterId] += (*points)[j].x;
                sumY[clusterId] += (*points)[j].y;

                (*points)[j].minDist = __DBL_MAX__;  /// reset min distance
            }
            #pragma omp barrier

            /// Compute the new centroids
            #pragma omp for
            for (unsigned int j=0; j < centroids.size(); j++) 
            {
                int clusterId =j;
                centroids[j].x = sumX[clusterId] / nPoints[clusterId];
                centroids[j].y = sumY[clusterId] / nPoints[clusterId];
            }
        }
    }

    /// Writing output to csv
    ofstream output_file;        
    output_file.open(output_file_name + "_parallel.csv");
    output_file << "x,y,cluster_id" << endl;

    for (vector<Point>::iterator it = points->begin(); it != points->end(); ++it) 
    {
        output_file << it->x << "," << it->y << "," << it->cluster << endl;
    }
    output_file.close();
}

int main(int argc, char *argv[]) 
{
    int no_iterations = 0;
    int no_clusters = 0;
    string input_file_name = "input.csv";

    /// Setting command line arguments
    if(argc < 4)
    {
        cerr << "Invalid options." << endl << 
        "<program> <Input Data Filename> <Output Filename> <Number of Iterations> <Number of Clusters> [-t <num_threads>]" << endl;
        exit(1);
    }

    input_file_name = argv[1];
    output_file_name = argv[2];
    no_iterations = atoi(argv[3]);
    no_clusters = atoi(argv[4]);

    if(argc == 7 && strcasecmp(argv[5], "-t") == 0){
        nthreads = atoi(argv[6]);
        omp_set_num_threads(nthreads);
    }

    /// Reading data file and populating points in a vector
    vector<Point> points = readcsv(input_file_name);
    
    /// Running K-Means algorithm for Customer Segmentation
    /// Serial Code
    double t1,t2;
    t1 = omp_get_wtime();
    kMeansClusteringSerial(&points, no_iterations, no_clusters);
    t2 = omp_get_wtime();

    cout << "\nExecution time for serial: " << (t2-t1) << endl;

    /// Parallel Code
    t1 = omp_get_wtime();
    kMeansClusteringParallel(&points, no_iterations, no_clusters);
    t2 = omp_get_wtime();

    cout << "\nExecution time for parallel: " << (t2-t1) << endl;

}